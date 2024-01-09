#include "../include/solver.h"

int counter = 0;
std::vector<std::vector<int>> grid_to_net_1;
std::vector<std::vector<int>> grid_to_net_2;
std::vector<std::vector<int>> grid_to_net_3;

Solver::Solver(int argc, char *argv[], int load_map)
{

    // Parse dos argumentos de entrada
    //parse_input_arguments(argc, argv);
    parse_input_arguments_v2(argc, argv, load_map);
    //print_parameters();   // DEBUG

    // Criar os diretorios para armazenar os arquivos das solucoes
    create_folders();

    // Ler e alocar memoria para geometria do grid vinda do arquivo de entrada
    read_and_allocate_tissue_geometry(this->tissue_filename);
    //print_cells_tissue();         // DEBUG
    
    // Numero de passos de tempo
    this->M = nearbyint(tmax / dt);

    // Cria o grafo a partir do arquivo de malha
    this->g1 = new Graph(network_filename_1, dx);
    this->g2 = new Graph(network_filename_2, dx);
    this->g3 = new Graph(network_filename_3, dx);

    // Alocando memoria para as celulas das arvores que serao resolvidas pelo metodo numerico
    allocate_memory_for_networks();

    // Inicilizar os tipos de celulas para as arvores
    set_type_cells_networks();
    //print_cells_network(1);       // DEBUG

    // Acoplar os terminais das arvores no tecido.
    // Ajustando o tipo de celula na estrutura do tecido
    set_type_cells_tissue();
    //write_tissue_vtk("outputs/tissue.vtk", this->n_tissue, this->n_tissue); // DEBUG
}

Solver::~Solver()
{
    // Free memory for the networks
    if (this->cell_types_network_1)
        free(this->cell_types_network_1);
    if (this->roi_network_1)
        free(this->roi_network_1);
    if (this->concentration_old_network_1)
        free(this->concentration_old_network_1);
    if (this->concentration_new_network_1)
        free(this->concentration_new_network_1);

    if (this->cell_types_network_2)
        free(this->cell_types_network_2);
    if (this->roi_network_2)
        free(this->roi_network_2);
    if (this->concentration_old_network_2)
        free(this->concentration_old_network_2);
    if (this->concentration_new_network_2)
        free(this->concentration_new_network_2);

    if (this->cell_types_network_3)
        free(this->cell_types_network_3);
    if (this->roi_network_3)
        free(this->roi_network_3);
    if (this->concentration_old_network_3)
        free(this->concentration_old_network_3);
    if (this->concentration_new_network_3)
        free(this->concentration_new_network_3);

    // Free memory for the tissue
    if (this->grid)
        free(this->grid);
    if (this->cell_types_tissue)
        free(this->cell_types_tissue);
    if (this->concentration_old_tissue)
        free(this->concentration_old_tissue);
    if (this->concentration_new_tissue)
        free(this->concentration_new_tissue);
    if (this->concentration_old_fibrosis)
        free(this->concentration_old_fibrosis);
    if (this->concentration_new_fibrosis)
        free(this->concentration_new_fibrosis);
    
    // Free memory for auxiliary vectors
    if (this->C_1D_old)
        free(this->C_1D_old);
    if (this->C_1D_new)
        free(this->C_1D_new);

    // Free memory for the Graph networks
    delete g1;
    delete g2;
    delete g3;
}

void Solver::parse_input_arguments(int argc, char *argv[])
{
    this->dt = atof(argv[1]);
    this->tmax = atof(argv[2]);
    this->network_filename_1 = argv[3];
    this->network_filename_2 = argv[4];
    this->network_filename_3 = argv[5];
    this->n_tissue = atoi(argv[6]);
    this->tissue_filename = argv[7];
    this->n_threads = atoi(argv[8]);
}

void Solver::parse_input_arguments_v2 (int argc, char *argv[], int load_map)
{
    this->dt = atof(argv[1]);
    this->tmax = atof(argv[2]);
    this->network_filename_1 = argv[3];
    this->network_filename_2 = argv[4];
    this->network_filename_3 = argv[5];
    this->tissue_filename = argv[6];
    this->output_dir = argv[7];
}

void Solver::print_parameters()
{
    printf("%s\n", LINE_2.c_str());
    printf("INPUT PARAMETERS\n");
    printf("%s\n", LINE_2.c_str());
    printf("dt = %g\n", this->dt);
    printf("tmax = %g\n", this->tmax);
    printf("network_filename_1 = %s\n", this->network_filename_1.c_str());
    printf("network_filename_2 = %s\n", this->network_filename_2.c_str());
    printf("network_filename_3 = %s\n", this->network_filename_3.c_str());
    printf("tissue_filename = %s\n", this->tissue_filename.c_str());
    printf("output_dir = %s\n", this->output_dir.c_str());
    printf("load_graph_map_from_file = %d\n", (int)this->load_map_from_file);
    printf("%s\n", LINE_2.c_str());
}

void Solver::print_cells_network(int net_id)
{
    int np;
    int *cell_type;
    int *roi;
    double *c_old;
    double *c_new;

    if (net_id == 1)
    {
        np = this->g1->getTotalNodes();
        cell_type = this->cell_types_network_1;
        roi = this->roi_network_1;
        c_old = this->concentration_new_network_1;
        c_new = this->concentration_new_network_1;

        printf("[Network 1]\n");
        for (int i = 0; i < np; i++)
        {
            printf("Node %d = [cell_type=%d, roi=%d, c_old=%g, c_new=%g]\n", i, cell_type[i], roi[i], c_old[i], c_new[i]);
        }
    }
    else if (net_id == 2)
    {
        np = this->g2->getTotalNodes();
        cell_type = this->cell_types_network_2;
        roi = this->roi_network_2;
        c_old = this->concentration_new_network_2;
        c_new = this->concentration_new_network_2;

        printf("[Network 2]\n");
        for (int i = 0; i < np; i++)
        {
            printf("Node %d = [cell_type=%d, roi=%d, c_old=%g, c_new=%g]\n", i, cell_type[i], roi[i], c_old[i], c_new[i]);
        }
    }
    else
    {
        np = this->g3->getTotalNodes();
        cell_type = this->cell_types_network_3;
        roi = this->roi_network_3;
        c_old = this->concentration_new_network_3;
        c_new = this->concentration_new_network_3;

        printf("[Network 3]\n");
        for (int i = 0; i < np; i++)
        {
            printf("Node %d = [cell_type=%d, roi=%d, c_old=%g, c_new=%g]\n", i, cell_type[i], roi[i], c_old[i], c_new[i]);
        }
    }
}

void Solver::print_cells_tissue()
{
    for (int j = 0; j < this->n_tissue; j++)
    {
        for (int i = 0; i < this->n_tissue; i++)
        {
            printf("%d ", this->grid[GET_IDX_2D(j, i, this->n_tissue, this->n_tissue)]);
        }
        printf("\n");
    }
}

void Solver::set_type_cells_networks()
{
    set_type_cells_using_network(this->g1, this->cell_types_network_1, this->concentration_old_network_1, this->concentration_new_network_1);
    set_type_cells_using_network(this->g2, this->cell_types_network_2, this->concentration_old_network_2, this->concentration_new_network_2);
    set_type_cells_using_network(this->g3, this->cell_types_network_3, this->concentration_old_network_3, this->concentration_new_network_3);
}

void Solver::allocate_memory_for_networks()
{

    // Allocate memory for the networks
    allocate_memory_for_network_cells(1, g1->getTotalNodes()); // Network 1
    allocate_memory_for_network_cells(2, g2->getTotalNodes()); // Network 2
    allocate_memory_for_network_cells(3, g3->getTotalNodes()); // Network 3

    // Allocate memory for the tissue
    //allocate_memory_for_tissue(); // Tissue + Ischemia

    // Allocate memory for auxiliary vectors
    this->C_1D_old = new double[n_1D]();
    this->C_1D_new = new double[n_1D]();
}

void Solver::set_type_cells_tissue()
{
    int i, j, k;

    // Initialize all cell types as normal
    int *cell_types = this->cell_types_tissue;
    for (i = 0; i < n_tissue; i++)
    {
        for (j = 0; j < n_tissue; j++)
        {
            cell_types[GET_IDX_2D(i, j, n_tissue, n_tissue)] = 0;
        }
    }

    // Check the terminals for each network and coupled them to the tissue
    set_type_cells_tissue_using_network(1);
    set_type_cells_tissue_using_network(2);
    set_type_cells_tissue_using_network(3);
}

// TODO: Trocar esta funcao para revisao do paper ...
// Busca pelo mais proximo
void Solver::busca_proximal(int id, int p[], int graph_id)
{
    int i, j;
    double aux;
    double dist = 100000.0;
    Node *ptr;

    if (graph_id == 0)
        ptr = g1->searchNode(id);
    else if (graph_id == 1)
        ptr = g2->searchNode(id);
    else if (graph_id == 2)
        ptr = g3->searchNode(id);

    // Loop esperto para nao iterar todas as celulas, mas somente as mais proximas do ponto da arvore
    for (i = (ptr->y - h_gadolinio) / h_gadolinio; i <= (ptr->y + h_gadolinio) / h_gadolinio; i++)
    {
        for (j = (ptr->x - h_gadolinio) / h_gadolinio; j <= (ptr->x + h_gadolinio) / h_gadolinio; j++)
        {
            aux = calc_norm(ptr->x, ptr->y, ptr->z, j * h_gadolinio, i * h_gadolinio, 0.);
            if (aux < dist)
            {
                dist = aux;
                p[0] = i;
                p[1] = j;
            }
        }
    }
}

void Solver::search_region(int id, int p[], int graph_id, double *peso)
{
    int i, j;
    double aux, dist = 100000.0;

    Node *ptr;

    if (graph_id == 0)
        ptr = g1->searchNode(id);
    else if (graph_id == 1)
        ptr = g2->searchNode(id);
    else if (graph_id == 2)
        ptr = g3->searchNode(id);

    // Loop esperto para nao iterar todas as celulas, mas somente as mais proximas do ponto da arvore
    for (i = (ptr->y - h_gadolinio) / h_gadolinio; i <= (ptr->y + h_gadolinio) / h_gadolinio; i++) {
        for (j = (ptr->x - h_gadolinio) / h_gadolinio; j <= (ptr->x + h_gadolinio) / h_gadolinio; j++) {
            aux = calc_norm(ptr->x, ptr->y, ptr->z, j * h_gadolinio, i * h_gadolinio, 0.);
            if (aux < dist) {
                dist = aux;
                p[0] = i;
                p[1] = j;
            }
        }
    }

    // Calculate if the network cell is healthy or ischemic
    int *cell_types; 
    int *roi; 

    // Network 1
    if (graph_id == 0) {
        cell_types = this->cell_types_network_1;
        roi = this->roi_network_1;
        if ((peso[GET_IDX_2D(p[0], p[1], n_tissue, n_tissue)] >= 0.75 && p[1] < n_tissue / 3. & p[0] < 10. * n_tissue / 16. && p[0] > 6. * n_tissue / 16.) && \
            (cell_types[id] == 1))
                roi[id] = 0;
        else
                roi[id] = 1;
    }
    // Network 2
    else if (graph_id == 1) {
        cell_types = this->cell_types_network_2;
        roi = this->roi_network_2;
        if ((peso[GET_IDX_2D(p[0], p[1], n_tissue, n_tissue)] >= 0.75 && p[1] < n_tissue / 3. && p[0] < 10. * n_tissue / 16. && p[0] > 6. * n_tissue / 16.) && \
            (cell_types[id] == 1))
                roi[id] = 0;
        else
                roi[id] = 1;
    }
    // Network 3
    else if (graph_id == 2) {
        cell_types = this->cell_types_network_3;
        roi = this->roi_network_3;
        if ((peso[GET_IDX_2D(p[0], p[1], n_tissue, n_tissue)] >= 0.75 && p[1] < n_tissue / 3. && p[0] < 10. * n_tissue / 16. && p[0] > 6. * n_tissue / 16.) && \
            (cell_types[id] == 1))
            roi[id] = 0;
        else
            roi[id] = 1;
    }

    //     else if(grid[GET_IDX_2D(p[0],p[1],n_tissue,n_tissue)] == 50)
    //         cOld[id].roi = 1;
    //     else
    //         cOld[id].roi = 2;
}

int Solver::number_tissue()
{
    int i, j;
    int N = 0;

    for (i = 0; i < n_tissue; i++)
    {
        for (j = 0; j < n_tissue; j++)
        {
            if (grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 128 || grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 50)
                N++;
        }
    }
    return N;
}

double Solver::integral_extra(int N)
{
    int i, j;

    double V = 0.;

    for (i = 1; i < n_tissue - 1; i++)
    {
        for (j = 1; j < n_tissue - 1; j++)
        {
            if (grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 128 || grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 50) // grid sadio ou grid infartado
            {
                V += cNew[GET_IDX_2D(i, j, n_tissue, n_tissue)].concentration;
            }
        }
    }

    //     V = V/10.; // soma das concentrações dividida pelo número de terminais
    V = V / N;
    return V;
}

double Solver::integral_intra()
{
    int j;
    double V = 0.;

    int np;
    int n_terminais = 0;
    ;

    np = g1->getTotalNodes();

    for (int j = 0; j < np; j++)
    {
        if (cOld[j].type == 1) // terminais
        {
            V += cNew[j].concentration;
            n_terminais++;
        }
    }

    np = g2->getTotalNodes();

    for (int j = 0; j < np; j++)
    {
        if (cOld[j].type == 1) // terminais
        {
            V += cNew[j].concentration;
            n_terminais++;
        }
    }

    np = g3->getTotalNodes();

    for (int j = 0; j < np; j++)
    {
        if (cOld[j].type == 1) // terminais
        {
            V += cNew[j].concentration;
            n_terminais++;
        }
    }

    return V / n_terminais;
}

void Solver::solve1D_source (double *C_1D, double V, double adv_1D, double difusao_1D) {
    int i;
    double *new_C_1D = new double[n_1D];

    C_1D[0] = V; //+= ou +?
                 //    fou1D(C_1D, P_1D, vswap, C_face_1D, v_1D);
    double dif_1D_global = 0.0;

    for (i = 1; i < n_1D; i++)
    {
        if (i != n_1D - 1)
        {
            adv_1D = (v_1D * C_1D[i] - v_1D * C_1D[i - 1]) / h_1D;

            difusao_1D = (dif_1D_global * (C_1D[i + 1] - C_1D[i]) - dif_1D_global * (C_1D[i] - C_1D[i - 1])) / (h_1D * h_1D);
        }

        else
        {
            adv_1D = (v_1D * C_1D[i] - v_1D * C_1D[i - 1]) / h_1D;
            difusao_1D = 0.;
        }

        if (i == n_1D - 1)
            new_C_1D[i] = (-adv_1D - rins * C_1D[i] - (dif_1D_global * (C_1D[i] - C_1D[i - 1])) / (h_1D * h_1D)) * dt / (phi_1D) + C_1D[i];

        else
            new_C_1D[i] = (difusao_1D - adv_1D - rins * C_1D[i]) * dt / (phi_1D) + C_1D[i];
    }

    for (i = 1; i < n_1D; i++)
    {
        C_1D[i] = new_C_1D[i];
    }
}

double Solver::calcula_fonte(double *C_1D, double t) {
    double fonte_new = (1. / (sigma * sqrt(2 * M_PI))) * exp(-0.5 * pow((t - T_peak) / sigma, 2));
    double sumidouro_new;

    //if (C_1D[n_1D - 1] - fonte_new < 0.0)
    //    fonte_new += 0.0;
    //else
    //    fonte_new += -r3 * (C_1D[n_1D - 1] - fonte_new);
    if (this->C_1D_old[n_1D - 1] - fonte_new < 0.0)
        fonte_new += 0.0;
    else
        fonte_new += -r3 * (this->C_1D_old[n_1D - 1] - fonte_new);

    return fonte_new;
}

double Solver::calcula_fonte_arteria(double *C_1D, int i, double dt, double C)
{
    double fonte_new = (1. / (sigma * sqrt(2 * M_PI))) * exp(-0.5 * pow((i * dt - T_peak) / sigma, 2));

    if (C_1D[n_1D - 1] - fonte_new < 0.)
        C = 0.;
    else
        C = r3 * (C_1D[n_1D - 1] - fonte_new);

    fonte_new += C;

    return fonte_new;
}

void Solver::solve(int load_map)
{
    char filename[500];

    int i, j, p, q, np, u[2];
    double dif_tissue = 0.0, aux = 0.0;

    // Signal intensity variables
    double SI_health_tissue = 0.0;
    double SI = 0.0;
    double SI_ischemic = 0.0;
    double SI_ischemic_tissue = 0.0;
    double SI_fibrose = 0.0;
    double SI_health_intra = 0.0;

    Node *ptr = g1->getListNodes();
    int N = count_number_active_tissue_cells();

#ifdef OUTPUT
    printf("[!] Solving transient problem ... \n");
    printf("[!] Progress\n");
    fflush(stdout);
#endif

    // Inicializa os termos de comunicacao do grid
    double *reaction = (double *)malloc(sizeof(double) * n_tissue * n_tissue);
    double *peso = (double *)malloc(sizeof(double) * n_tissue * n_tissue);

    // Initialize the solution of the model
    set_initial_conditions_tissue(reaction, peso);

    // Parte 1D da recirculação
    double difusao_1D = 0.0, adv_1D = 0.0, V = 0.0;
    //double *C_1D = (double*)malloc(sizeof(double)*n_1D*2);
    double *C_1D = new double[n_1D*2];
    for (i = 0; i < 2*n_1D; i++) {
        C_1D[i] = 0.0;
    }

    // Calcular o vetor peso
    contorno_perm_het(peso, grid);
    peso = perm_het(peso, grid);

    // Atualizar o vetor 'roi' nas arvores (define se os terminais serao healthy ou ischemia)
    for (i = 0; i < 3; i++) {
        if (i == 0) {
            np = g1->getTotalNodes();
        }
        else if (i == 1) {
            np = g2->getTotalNodes();
        }
        else if (i == 2) {
            np = g3->getTotalNodes();
        }
        for (j = 0; j < np; j++) {
            search_region(j, u, i, peso);
        }
    }

    // Save the signal intensity for each print rate timestep
    double *SI_health_vec = (double*)malloc(sizeof(double)*(M / SAVEEACH));
    double *SI_ischemic_vec = (double*)malloc(sizeof(double)*(M / SAVEEACH));
    double *SI_fibrose_vec = (double*)malloc(sizeof(double)*(M / SAVEEACH));

    // ===================================================================================
    // Time loop starts
    for (int i = 0; i < M; i++) {
        double t = i * dt;

        // Reset this variable
        SI_health_intra = 0.0;
        
        // 1) RECIRCULACAO:
        //   1.1) Calcula a média da concentração no domínio
        //V = integral_intra();
        V = integral_intra_v2();

        //   1.2) Resolve a parte 1D (Recirculação), 
        //        com V sendo a condição de contorno prescrita no lado esquerdo
        solve1D_source(C_1D, V, adv_1D, difusao_1D);
        //solve1D_source_v2(C_1D, V, adv_1D, difusao_1D);

        //   1.3) Calcula o termo de fonte que será atribudo como fluxo no contorno
        double fonte_new = calcula_fonte(C_1D, t);
        // Print the progress of the solution
        #ifdef OUTPUT
        print_progress(i, M);
        #endif
        
        // Write the solution to .vtk file
        #ifdef VTK
        if (i % SAVEEACH == 0) {
            // Write the networks    
            write_vtk_file_networks(this->g1, this->cell_types_network_1, this->roi_network_1, this->concentration_new_network_1, 0, i);
            write_vtk_file_networks(this->g2, this->cell_types_network_2, this->roi_network_2, this->concentration_new_network_2, 1, i);
            write_vtk_file_networks(this->g3, this->cell_types_network_3, this->roi_network_3, this->concentration_new_network_3, 2, i);
            
            // Write the tissue
            sprintf(filename, "%s/tissue/tissue%d.vtk", this->output_dir.c_str(), i);
            write_vtk_file_tissue(filename, n_tissue, n_tissue);
        }
        #endif

        // 2) CONCENTRACAO
            // Resolver as equacoes para as arvores primeiro ...
        solve_concentration_for_network(g1, 0, this->concentration_old_network_1, \
                                            this->concentration_new_network_1, \
                                            this->cell_types_network_1, \
                                            this->roi_network_1, \
                                            reaction, \
                                            fonte_new, \
                                            SI_health_intra);
        solve_concentration_for_network(g2, 1, this->concentration_old_network_2, \
                                            this->concentration_new_network_2, \
                                            this->cell_types_network_2, \
                                            this->roi_network_2, \
                                            reaction, \
                                            fonte_new, \
                                            SI_health_intra);
        solve_concentration_for_network(g3, 2, this->concentration_old_network_3, \
                                            this->concentration_new_network_3, \
                                            this->cell_types_network_3, \
                                            this->roi_network_3, \
                                            reaction, \
                                            fonte_new, \
                                            SI_health_intra);

            // Resolver as equacoes para o tecido
        solve_concentration_for_tissue(reaction);

        // Cálculo da Intensidade do Sinal do Tecido (extravascular)
        calculate_signal_intensity(i, peso, SI_ischemic, SI_health_tissue, SI_fibrose, SI, SI_ischemic_vec, SI_fibrose_vec, SI_health_vec);
        
        // Pular para a proxima iteracao
        next_timestep();
    }

    post_processing(SI_health_vec, SI_ischemic_vec, SI_fibrose_vec);

    // Liberar memoria
    free(SI_health_vec);
    free(SI_ischemic_vec);
    free(SI_fibrose_vec);
    free(C_1D);
    free(reaction);
    free(peso);
}

void Solver::calcula_reaction(double *reaction, int k, int i, int j, int graph_id) {
    
    double *c_old_network = NULL;
    double *c_new_network = NULL;
    double *c_old_tissue = this->concentration_old_tissue;
    double *c_new_tissue = this->concentration_new_tissue;
    if (graph_id == 0) {
        c_old_network = this->concentration_old_network_1;
        c_new_network = this->concentration_new_network_1;
    }
    else if (graph_id == 1) {
        c_old_network = this->concentration_old_network_2;
        c_new_network = this->concentration_new_network_2;
    }
    else if (graph_id == 2) {
        c_old_network = this->concentration_old_network_3;
        c_new_network = this->concentration_new_network_3;
    }

    if (c_old_network[k] < c_old_tissue[GET_IDX_2D(i, j, n_tissue, n_tissue)]) {
        reaction[GET_IDX_2D(i, j, n_tissue, n_tissue)] = 0.;
    }
    else {
        if (!(grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 50)) {
            reaction[GET_IDX_2D(i, j, n_tissue, n_tissue)] = -r * (c_old_network[k] - c_old_tissue[GET_IDX_2D(i, j, n_tissue, n_tissue)]);
        }
        else {
            reaction[GET_IDX_2D(i, j, n_tissue, n_tissue)] = -r_problema * ( c_old_network[k] - c_old_tissue[GET_IDX_2D(i, j, n_tissue, n_tissue)]);
        }
    }
}

void Solver::contorno_perm_het(double *p, int *grid)
{
    int i, j;
    for (i = 1; i < n_tissue - 1; i++) {
        for (j = 1; j < n_tissue - 1; j++) {
            if (grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 128 || \
                grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 50) {
                if (grid[GET_IDX_2D(i + 1, j, n_tissue, n_tissue)] == 255 || \
                    grid[GET_IDX_2D(i - 1, j, n_tissue, n_tissue)] == 255 || \
                    grid[GET_IDX_2D(i, j + 1, n_tissue, n_tissue)] == 255 || \
                    grid[GET_IDX_2D(i, j - 1, n_tissue, n_tissue)] == 255)
                        p[GET_IDX_2D(i, j, n_tissue, n_tissue)] = 0.0;
                else if (grid[GET_IDX_2D(i + 1, j, n_tissue, n_tissue)] == 0 || \
                         grid[GET_IDX_2D(i - 1, j, n_tissue, n_tissue)] == 0 || \
                         grid[GET_IDX_2D(i, j + 1, n_tissue, n_tissue)] == 0 || \
                         grid[GET_IDX_2D(i, j - 1, n_tissue, n_tissue)] == 0)
                        p[GET_IDX_2D(i, j, n_tissue, n_tissue)] = 1.0;
            }
        }
    }
}

double *Solver::perm_het(double *p, int *grid)
{
    int i, j, m = 0;
    double max = 0., aux1, aux2, erro = 1.;

    /* ----- Calculo das pressões -----*/
    while (erro > 0.000001)
    {
        max = 0.;
        for (i = 1; i < n_tissue - 1; i++)
        {
            for (j = 1; j < n_tissue - 1; j++)
            {
                if ((grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 128 || grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 50) && !(grid[GET_IDX_2D(i + 1, j, n_tissue, n_tissue)] == 255 || grid[GET_IDX_2D(i - 1, j, n_tissue, n_tissue)] == 255 || grid[GET_IDX_2D(i, j - 1, n_tissue, n_tissue)] == 255 || grid[GET_IDX_2D(i, j + 1, n_tissue, n_tissue)] == 255) && !(grid[GET_IDX_2D(i + 1, j, n_tissue, n_tissue)] == 0 || grid[GET_IDX_2D(i - 1, j, n_tissue, n_tissue)] == 0 || grid[GET_IDX_2D(i, j - 1, n_tissue, n_tissue)] == 0 || grid[GET_IDX_2D(i, j + 1, n_tissue, n_tissue)] == 0))
                {
                    aux1 = p[GET_IDX_2D(i, j, n_tissue, n_tissue)];

                    p[GET_IDX_2D(i, j, n_tissue, n_tissue)] = (p[GET_IDX_2D(i, j + 1, n_tissue, n_tissue)] + p[GET_IDX_2D(i, j - 1, n_tissue, n_tissue)] + p[GET_IDX_2D(i + 1, j, n_tissue, n_tissue)] + p[GET_IDX_2D(i - 1, j, n_tissue, n_tissue)]) / 4.;

                    aux2 = fabs(p[GET_IDX_2D(i, j, n_tissue, n_tissue)] - aux1);

                    if (aux2 > max)
                        max = aux2;
                }
                else if (grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 255 || grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 0)
                {
                    p[GET_IDX_2D(i, j, n_tissue, n_tissue)] = -1.;
                }
            }
        }

        if (max < erro)
            erro = max;
        m++;
        //        printf("erro =  %f \n", erro);
    }
    //    printf("%d \n", m);
    return p;
}

double Solver::signal_intensity_healthy (int k, int * tecido, Cell * S, double * peso, double SI_health)
{
    int i, j;
//     double fator = 130.;
//     real SI = 0.;
//     real SI_intra = 0.;
//     real SI_extra = 0.;
//     real SI_intra_normal = 0.;
//     real SI_extra_normal = 0.;
//     real SI_fibrosis = 0.;
//     real SI_intra_extra = 0.;
//     real SI_intra_infarto = 0.;
//     real SI_extra_infarto = 0.;
//     real SI_extra_fibrose = 0.;
    SI_health = 0.;

    for(i = 0; i <= n_tissue-1; i++)
    {
        for(j = 0; j <= n_tissue-1; j++)
        {
            if(tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 128 || tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 50)
            {
                if(tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 50)
                {
//                     SI_ischemic += (1.0 - phi) * lambda_infarto * S[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration;
//                     SI_ischemic += (phi)*C[GET_IDX_3D(vswap,i,j,2,n,n)] + ((1.0 - phi)*lambda_infarto*S[GET_IDX_3D(vswap,i,j,2,n,n)]) + ((1.0 - phi)*lambda_infarto*lambda_fibrose*C_fibrose[GET_IDX_3D(vswap,i,j,2,n,n)]);
//                     SI_intra_infarto += phi*C[GET_IDX_3D(vswap,i,j,2,n,n)];
//                     SI_extra_infarto += (1.0 - phi)*S[GET_IDX_3D(vswap,i,j,2,n,n)];
// 
//                     SI += (phi)*C[GET_IDX_3D(vswap,i,j,2,n,n)] + ((1.0 - phi)*lambda_infarto*S[GET_IDX_3D(vswap,i,j,2,n,n)]) + ((1.0 - phi)*lambda_infarto*lambda_fibrose*C_fibrose[GET_IDX_3D(vswap,i,j,2,n,n)]);
//                     SI_intra += (phi*C[GET_IDX_3D(vswap,i,j,2,n,n)]); //global
//                     SI_extra += ((1.0 - phi)*lambda_infarto*S[GET_IDX_3D(vswap,i,j,2,n,n)]); // global
// 
//                     SI_extra_fibrose += (1.0 - phi)*lambda_infarto*S[GET_IDX_3D(vswap,i,j,2,n,n)] + (1.0 - phi)*lambda_infarto*lambda_fibrose*C_fibrose[GET_IDX_3D(vswap,i,j,2,n,n)]; //global
//                     SI_fibrosis += (1.0 - phi)*lambda_infarto*lambda_fibrose*(C_fibrose[GET_IDX_3D(vswap,i,j,2,n,n)]); // global
                }


                else if(peso[GET_IDX_2D(i,j,n_tissue,n_tissue)] >= 0.75 && j < n_tissue/3. && i < 10.*n_tissue/16. && i > 6.*n_tissue/16.)
                {
                    SI_health += /*(phi)*C[GET_IDX_3D(vswap,i,j,2,n,n)] + (1.0 - phi) * lambda **/ S[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration; /*+ ((1.0 - phi)*lambda*lambda_fibrose*C_fibrose[GET_IDX_3D(vswap,i,j,2,n,n)])*/
//                     SI_intra_normal += phi*C[GET_IDX_3D(vswap,i,j,2,n,n)];
//                     SI_extra_normal += (1-phi)*S[GET_IDX_3D(vswap,i,j,2,n,n)];
                }
            }
        }
    }

//     printf("%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f \n", k*dt, fator*SI, fator*SI_intra, fator*SI_extra, fator*SI_ischemic, fator*SI_health, fator*SI_fibrosis, fator*SI_extra_fibrose);
    
    
    
//     printf("%.4f %.4f \n", k*dt, fator*SI_health);

    return SI_health;
//     SI = 0.;
//     SI_intra = 0.;
//     SI_intra_infarto = 0.;
//     SI_intra_normal = 0.;
//     SI_extra = 0.;
//     SI_extra_infarto = 0.;
//     SI_extra_normal = 0.;
//     SI_health = 0.;
//     SI_ischemic = 0.;
//     SI_fibrosis = 0.;
//     SI_intra_extra = 0.;
//     SI_extra_fibrose = 0.;
}

double Solver::signal_intensity_healthy_v2 (int k, int *grid, double *c_new, double *peso, double SI_health) {
    int i, j;
    SI_health = 0.0;

    for (i = 0; i <= n_tissue - 1; i++) {
        for (j = 0; j <= n_tissue - 1; j++) {
            if (grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 128 || grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 50) {
                if (grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 50) {
                    // ???
                }
                else if (peso[GET_IDX_2D(i, j, n_tissue, n_tissue)] >= 0.75 && j < n_tissue / 3. && i < 10. * n_tissue / 16. && i > 7. * n_tissue / 16.) {
                    SI_health += c_new[GET_IDX_2D(i, j, n_tissue, n_tissue)];
                }
            }
        }
    }

    return SI_health;
}

double Solver::signal_intensity_ischemic(int k, int *grid, Cell *S, double *peso, double SI_ischemic)
{
    int i, j;
    //     double fator = 130.;
    //     real SI = 0.;
    //     real SI_intra = 0.;
    //     real SI_extra = 0.;
    //     real SI_intra_normal = 0.;
    //     real SI_extra_normal = 0.;
    //     real SI_fibrosis = 0.;
    //     real SI_intra_extra = 0.;
    //     real SI_intra_infarto = 0.;
    //     real SI_extra_infarto = 0.;
    //     real SI_extra_fibrose = 0.;
    SI_ischemic = 0.;

    for (i = 0; i <= n_tissue - 1; i++)
    {
        for (j = 0; j <= n_tissue - 1; j++)
        {
            if (grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 128 || grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 50)
            {
                //                 if(grid[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 50)//  && peso[GET_IDX_2D(i,j,n_tissue,n_tissue)] >= 0.95)

                //                 if(peso[GET_IDX_2D(i,j,n_tissue,n_tissue)] >= 0.85 && j > n_tissue/2. && i > 43.*n_tissue/100. && i < 57.*n_tissue/100.)

                if ((j > n_tissue / 2. && j < 74. * n_tissue / 100. && i > 43. * n_tissue / 100. && i < 57. * n_tissue / 100.) || (j == 74. * n_tissue / 100. && i > 44. * n_tissue / 100. && i < 56. * n_tissue / 100.) || (j == 75. * n_tissue / 100. && i > 45. * n_tissue / 100. && i < 55. * n_tissue / 100.) || (j == 76. * n_tissue / 100. && i > 46. * n_tissue / 100. && i < 54. * n_tissue / 100.) || (j == 77. * n_tissue / 100. && i > 48. * n_tissue / 100. && i < 54. * n_tissue / 100.))
                {
                    SI_ischemic += /*(1.0 - phi) * lambda_infarto * */ S[GET_IDX_2D(i, j, n_tissue, n_tissue)].concentration;
                    //                     SI_ischemic += (phi)*C[GET_IDX_3D(vswap,i,j,2,n,n)] + ((1.0 - phi)*lambda_infarto*S[GET_IDX_3D(vswap,i,j,2,n,n)]) + ((1.0 - phi)*lambda_infarto*lambda_fibrose*C_fibrose[GET_IDX_3D(vswap,i,j,2,n,n)]);
                    //                     SI_intra_infarto += phi*C[GET_IDX_3D(vswap,i,j,2,n,n)];
                    //                     SI_extra_infarto += (1.0 - phi)*S[GET_IDX_3D(vswap,i,j,2,n,n)];
                    //
                    //                     SI += (phi)*C[GET_IDX_3D(vswap,i,j,2,n,n)] + ((1.0 - phi)*lambda_infarto*S[GET_IDX_3D(vswap,i,j,2,n,n)]) + ((1.0 - phi)*lambda_infarto*lambda_fibrose*C_fibrose[GET_IDX_3D(vswap,i,j,2,n,n)]);
                    //                     SI_intra += (phi*C[GET_IDX_3D(vswap,i,j,2,n,n)]); //global
                    //                     SI_extra += ((1.0 - phi)*lambda_infarto*S[GET_IDX_3D(vswap,i,j,2,n,n)]); // global
                    //
                    //                     SI_extra_fibrose += (1.0 - phi)*lambda_infarto*S[GET_IDX_3D(vswap,i,j,2,n,n)] + (1.0 - phi)*lambda_infarto*lambda_fibrose*C_fibrose[GET_IDX_3D(vswap,i,j,2,n,n)]; //global
                    //                     SI_fibrosis += (1.0 - phi)*lambda_infarto*lambda_fibrose*(C_fibrose[GET_IDX_3D(vswap,i,j,2,n,n)]); // global
                }

                else if (peso[GET_IDX_2D(i, j, n_tissue, n_tissue)] >= 0.75 && j < n_tissue / 3. && i < 10. * n_tissue / 16. && i > 6. * n_tissue / 16.)
                {
                    //                     SI_health += /*(phi)*C[GET_IDX_3D(vswap,i,j,2,n,n)] + */((1.0 - phi) * lambda * S[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration); /*+ ((1.0 - phi)*lambda*lambda_fibrose*C_fibrose[GET_IDX_3D(vswap,i,j,2,n,n)])*/
                    //                     SI_intra_normal += phi*C[GET_IDX_3D(vswap,i,j,2,n,n)];
                    //                     SI_extra_normal += (1-phi)*S[GET_IDX_3D(vswap,i,j,2,n,n)];
                }
            }
        }
    }

    //     printf("%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f \n", k*dt, fator*SI, fator*SI_intra, fator*SI_extra, fator*SI_ischemic, fator*SI_health, fator*SI_fibrosis, fator*SI_extra_fibrose);

    //     printf("%.4f %.4f \n", k*dt, fator*SI_health);

    return SI_ischemic;
    //     SI = 0.;
    //     SI_intra = 0.;
    //     SI_intra_infarto = 0.;
    //     SI_intra_normal = 0.;
    //     SI_extra = 0.;
    //     SI_extra_infarto = 0.;
    //     SI_extra_normal = 0.;
    //     SI_health = 0.;
    //     SI_ischemic = 0.;
    //     SI_fibrosis = 0.;
    //     SI_intra_extra = 0.;
    //     SI_extra_fibrose = 0.;
}

double Solver::signal_intensity_ischemic_v2 (int k, int *grid, double *c_new, double *peso, double SI_ischemic) {
    int i, j;
    SI_ischemic = 0.;

    for (i = 0; i <= n_tissue - 1; i++) {
        for (j = 0; j <= n_tissue - 1; j++) {
            if (grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 128 || grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 50) {
                if ((j > n_tissue / 2. && j < 74. * n_tissue / 100. && i > 43. * n_tissue / 100. && i < 57. * n_tissue / 100.) || (j == 74. * n_tissue / 100. && i > 44. * n_tissue / 100. && i < 56. * n_tissue / 100.) || (j == 75. * n_tissue / 100. && i > 45. * n_tissue / 100. && i < 55. * n_tissue / 100.) || (j == 76. * n_tissue / 100. && i > 46. * n_tissue / 100. && i < 54. * n_tissue / 100.) || (j == 77. * n_tissue / 100. && i > 48. * n_tissue / 100. && i < 54. * n_tissue / 100.)) {
                    SI_ischemic += c_new[GET_IDX_2D(i, j, n_tissue, n_tissue)];
                }
                else if (peso[GET_IDX_2D(i, j, n_tissue, n_tissue)] >= 0.75 && j < n_tissue / 3. && i < 10. * n_tissue / 16. && i > 6. * n_tissue / 16.) {
                    
                }
            }
        }
    }

    return SI_ischemic;
}

double Solver::signal_intensity(int k, int *grid, Cell *S, Cell *C_fibrose, double *peso, double SI_fibrose, double SI_ischemic, double SI)
{
    int i, j;
    SI = 0.;
    //     real SI_intra = 0.;
    //     real SI_extra = 0.;
    //     real SI_intra_normal = 0.;
    //     real SI_extra_normal = 0.;
    //     real SI_fibrosis = 0.;
    //     real SI_intra_extra = 0.;
    //     real SI_intra_infarto = 0.;
    //     real SI_extra_infarto = 0.;
    //     real SI_extra_fibrose = 0.;
    SI_fibrose = 0.;
    SI_ischemic = 0.;

    for (i = 0; i <= n_tissue - 1; i++)
    {
        for (j = 0; j <= n_tissue - 1; j++)
        {
            if (grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 128 || grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 50)
            {
                //                 if(grid[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 50)//  && peso[GET_IDX_2D(i,j,n_tissue,n_tissue)] >= 0.95)
                //                 if(peso[GET_IDX_2D(i,j,n_tissue,n_tissue)] >= 0.85 && j > n_tissue/2. && i > 43.*n_tissue/100. && i < 57.*n_tissue/100.)
                if ((j > n_tissue / 2. && j < 74. * n_tissue / 100. && i > 43. * n_tissue / 100. && i < 57. * n_tissue / 100.) || (j == 74. * n_tissue / 100. && i > 44. * n_tissue / 100. && i < 56. * n_tissue / 100.) || (j == 75. * n_tissue / 100. && i > 45. * n_tissue / 100. && i < 55. * n_tissue / 100.) || (j == 76. * n_tissue / 100. && i > 46. * n_tissue / 100. && i < 54. * n_tissue / 100.) || (j == 77. * n_tissue / 100. && i > 48. * n_tissue / 100. && i < 54. * n_tissue / 100.))
                {
                    //                     SI_fibrose += /*(1.0 - phi) * lambda_infarto * */S[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration;
                    SI_ischemic += ((1.0 - phi) * lambda_infarto * S[GET_IDX_2D(i, j, n_tissue, n_tissue)].concentration);
                    //                     SI_intra_infarto += phi*C[GET_IDX_3D(vswap,i,j,2,n,n)];
                    //                     SI_extra_infarto += (1.0 - phi)*S[GET_IDX_3D(vswap,i,j,2,n,n)];
                    //
                    //                     SI += ((1.0 - phi)*lambda_infarto*S[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration) + ((1.0 - phi)*lambda_infarto*lambda_fibrose*C_fibrose[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration);

                    SI += S[GET_IDX_2D(i, j, n_tissue, n_tissue)].concentration + C_fibrose[GET_IDX_2D(i, j, n_tissue, n_tissue)].concentration;

                    //                     SI_intra += (phi*C[GET_IDX_3D(vswap,i,j,2,n,n)]); //global
                    //                     SI_extra += ((1.0 - phi)*lambda_infarto*S[GET_IDX_3D(vswap,i,j,2,n,n)]); // global
                    //
                    //                     SI_extra_fibrose += (1.0 - phi)*lambda_infarto*S[GET_IDX_3D(vswap,i,j,2,n,n)] + (1.0 - phi)*lambda_infarto*lambda_fibrose*C_fibrose[GET_IDX_3D(vswap,i,j,2,n,n)]; //global
                    SI_fibrose += (1.0 - phi) * lambda_infarto * lambda_fibrose * (C_fibrose[GET_IDX_2D(i, j, n_tissue, n_tissue)].concentration); // global
                }

                else if (peso[GET_IDX_2D(i, j, n_tissue, n_tissue)] >= 0.75 && j < n_tissue / 3. && i < 10. * n_tissue / 16. && i > 6. * n_tissue / 16.)
                {
                    //                     SI_health += /*(phi)*C[GET_IDX_3D(vswap,i,j,2,n,n)] + */((1.0 - phi) * lambda * S[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration); /*+ ((1.0 - phi)*lambda*lambda_fibrose*C_fibrose[GET_IDX_3D(vswap,i,j,2,n,n)])*/
                    //                     SI_intra_normal += phi*C[GET_IDX_3D(vswap,i,j,2,n,n)];
                    //                     SI_extra_normal += (1-phi)*S[GET_IDX_3D(vswap,i,j,2,n,n)];
                }
            }
        }
    }

    //     printf("%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f \n", k*dt, fator*SI, fator*SI_intra, fator*SI_extra, fator*SI_ischemic, fator*SI_health, fator*SI_fibrosis, fator*SI_extra_fibrose);

    //     printf("%.4f %.4f \n", k*dt, fator*SI_health);

    return SI;
    //     SI = 0.;
    //     SI_intra = 0.;
    //     SI_intra_infarto = 0.;
    //     SI_intra_normal = 0.;
    //     SI_extra = 0.;
    //     SI_extra_infarto = 0.;
    //     SI_extra_normal = 0.;
    //     SI_health = 0.;
    //     SI_ischemic = 0.;
    //     SI_fibrosis = 0.;
    //     SI_intra_extra = 0.;
    //     SI_extra_fibrose = 0.;
}

double Solver::signal_intensity_v2 (int k, int *grid, double *c_new, double *c_new_fibrosis, double *peso, double SI_fibrose, double SI_ischemic, double SI) {
    int i, j;
    SI = 0.;
    SI_fibrose = 0.;
    SI_ischemic = 0.;

    for (i = 0; i <= n_tissue - 1; i++) {
        for (j = 0; j <= n_tissue - 1; j++) {
            if (grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 128 || grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 50) {
                if ((j > n_tissue / 2. && j < 74. * n_tissue / 100. && i > 43. * n_tissue / 100. && i < 57. * n_tissue / 100.) || (j == 74. * n_tissue / 100. && i > 44. * n_tissue / 100. && i < 56. * n_tissue / 100.) || (j == 75. * n_tissue / 100. && i > 45. * n_tissue / 100. && i < 55. * n_tissue / 100.) || (j == 76. * n_tissue / 100. && i > 46. * n_tissue / 100. && i < 54. * n_tissue / 100.) || (j == 77. * n_tissue / 100. && i > 48. * n_tissue / 100. && i < 54. * n_tissue / 100.)) {
                    SI_ischemic += ((1.0 - phi) * lambda_infarto * c_new[GET_IDX_2D(i, j, n_tissue, n_tissue)]);
                    
                    SI += c_new[GET_IDX_2D(i, j, n_tissue, n_tissue)] + c_new_fibrosis[GET_IDX_2D(i, j, n_tissue, n_tissue)];

                    SI_fibrose += (1.0 - phi) * lambda_infarto * lambda_fibrose * (c_new_fibrosis[GET_IDX_2D(i, j, n_tissue, n_tissue)]); // global
                }
                else if (peso[GET_IDX_2D(i, j, n_tissue, n_tissue)] >= 0.75 && j < n_tissue / 3. && i < 10. * n_tissue / 16. && i > 6. * n_tissue / 16.) {

                }
            }
        }
    }

    return SI;
}

double Solver::signal_intensity_fibrose(int k, int *grid, Cell *S, Cell *C_fibrose, double *peso, double SI_fibrose, double SI_ischemic, double SI)
{
    int i, j;
    SI = 0.;
    //     real SI_intra = 0.;
    //     real SI_extra = 0.;
    //     real SI_intra_normal = 0.;
    //     real SI_extra_normal = 0.;
    //     real SI_fibrosis = 0.;
    //     real SI_intra_extra = 0.;
    //     real SI_intra_infarto = 0.;
    //     real SI_extra_infarto = 0.;
    //     real SI_extra_fibrose = 0.;
    SI_fibrose = 0.;
    SI_ischemic = 0.;

    for (i = 0; i <= n_tissue - 1; i++)
    {
        for (j = 0; j <= n_tissue - 1; j++)
        {
            if (grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 128 || grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 50)
            {
                //                 if(grid[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 50)//  && peso[GET_IDX_2D(i,j,n_tissue,n_tissue)] >= 0.95)
                //                 if(peso[GET_IDX_2D(i,j,n_tissue,n_tissue)] >= 0.85 && j > n_tissue/2. && i > 43.*n_tissue/100. && i < 57.*n_tissue/100.)
                if ((j > n_tissue / 2. && j < 74. * n_tissue / 100. && i > 43. * n_tissue / 100. && i < 57. * n_tissue / 100.) || (j == 74. * n_tissue / 100. && i > 44. * n_tissue / 100. && i < 56. * n_tissue / 100.) || (j == 75. * n_tissue / 100. && i > 45. * n_tissue / 100. && i < 55. * n_tissue / 100.) || (j == 76. * n_tissue / 100. && i > 46. * n_tissue / 100. && i < 54. * n_tissue / 100.) || (j == 77. * n_tissue / 100. && i > 48. * n_tissue / 100. && i < 54. * n_tissue / 100.))
                {
                    SI_fibrose += /*(1.0 - phi) * lambda_infarto * */ C_fibrose[GET_IDX_2D(i, j, n_tissue, n_tissue)].concentration;
                    SI_ischemic += ((1.0 - phi) * lambda_infarto * S[GET_IDX_2D(i, j, n_tissue, n_tissue)].concentration);
                    //                     SI_intra_infarto += phi*C[GET_IDX_3D(vswap,i,j,2,n,n)];
                    //                     SI_extra_infarto += (1.0 - phi)*S[GET_IDX_3D(vswap,i,j,2,n,n)];
                    //
                    //                     SI += ((1.0 - phi)*lambda_infarto*S[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration) + ((1.0 - phi)*lambda_infarto*lambda_fibrose*C_fibrose[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration);

                    //                     SI += S[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration + C_fibrose[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration;

                    //                     SI_intra += (phi*C[GET_IDX_3D(vswap,i,j,2,n,n)]); //global
                    //                     SI_extra += ((1.0 - phi)*lambda_infarto*S[GET_IDX_3D(vswap,i,j,2,n,n)]); // global
                    //
                    //                     SI_extra_fibrose += (1.0 - phi)*lambda_infarto*S[GET_IDX_3D(vswap,i,j,2,n,n)] + (1.0 - phi)*lambda_infarto*lambda_fibrose*C_fibrose[GET_IDX_3D(vswap,i,j,2,n,n)]; //global
                    //                     SI_fibrose += (1.0 - phi)*lambda_infarto*lambda_fibrose*(C_fibrose[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration); // global
                }

                else if (peso[GET_IDX_2D(i, j, n_tissue, n_tissue)] >= 0.75 && j < n_tissue / 3. && i < 10. * n_tissue / 16. && i > 6. * n_tissue / 16.)
                {
                    //                     SI_health += /*(phi)*C[GET_IDX_3D(vswap,i,j,2,n,n)] + */((1.0 - phi) * lambda * S[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration); /*+ ((1.0 - phi)*lambda*lambda_fibrose*C_fibrose[GET_IDX_3D(vswap,i,j,2,n,n)])*/
                    //                     SI_intra_normal += phi*C[GET_IDX_3D(vswap,i,j,2,n,n)];
                    //                     SI_extra_normal += (1-phi)*S[GET_IDX_3D(vswap,i,j,2,n,n)];
                }
            }
        }
    }

    //     printf("%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f \n", k*dt, fator*SI, fator*SI_intra, fator*SI_extra, fator*SI_ischemic, fator*SI_health, fator*SI_fibrosis, fator*SI_extra_fibrose);

    //     printf("%.4f %.4f \n", k*dt, fator*SI_health);

    return SI_fibrose;
    //     SI = 0.;
    //     SI_intra = 0.;
    //     SI_intra_infarto = 0.;
    //     SI_intra_normal = 0.;
    //     SI_extra = 0.;
    //     SI_extra_infarto = 0.;
    //     SI_extra_normal = 0.;
    //     SI_health = 0.;
    //     SI_ischemic = 0.;
    //     SI_fibrosis = 0.;
    //     SI_intra_extra = 0.;
    //     SI_extra_fibrose = 0.;
}

double Solver::signal_intensity_fibrose_v2 (int k, int *grid, double *c_new, double *c_new_fibrosis, double *peso, double SI_fibrose, double SI_ischemic, double SI) {
    int i, j;
    SI = 0.0;
    SI_fibrose = 0.0;
    SI_ischemic = 0.0;

    for (i = 0; i <= n_tissue - 1; i++) {
        for (j = 0; j <= n_tissue - 1; j++) {
            if (grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 128 || grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 50) {
                if ((j > n_tissue / 2. && j < 74. * n_tissue / 100. && i > 43. * n_tissue / 100. && i < 57. * n_tissue / 100.) || (j == 74. * n_tissue / 100. && i > 44. * n_tissue / 100. && i < 56. * n_tissue / 100.) || (j == 75. * n_tissue / 100. && i > 45. * n_tissue / 100. && i < 55. * n_tissue / 100.) || (j == 76. * n_tissue / 100. && i > 46. * n_tissue / 100. && i < 54. * n_tissue / 100.) || (j == 77. * n_tissue / 100. && i > 48. * n_tissue / 100. && i < 54. * n_tissue / 100.)) {
                    SI_fibrose += c_new_fibrosis[GET_IDX_2D(i, j, n_tissue, n_tissue)];
                    SI_ischemic += ((1.0 - phi) * lambda_infarto * c_new[GET_IDX_2D(i, j, n_tissue, n_tissue)]);
                }
                else if (peso[GET_IDX_2D(i, j, n_tissue, n_tissue)] >= 0.75 && j < n_tissue / 3. && i < 10. * n_tissue / 16. && i > 6. * n_tissue / 16.) {
                
                }
            }
        }
    }

    return SI_fibrose;
}

void Solver::read_and_allocate_tissue_geometry(std::string filename)
{
    int p = 0;
    ifstream in(filename.c_str());
    if (!in)
    {
        fprintf(stderr, "[-] ERROR! Cannot open tissue DAT file!\n");
        exit(EXIT_FAILURE);
    }

    int nx, ny;
    double h;
    in >> nx >> ny >> h;

    // Update variables
    this->n_tissue = nx;
    this->h_gadolinio = h;
    
    // Allocate memory for the tissue mesh
    allocate_memory_for_tissue();

    // Read the tissue mesh
    for (int j = 0; j < nx; j++)
    {
        for (int i = 0; i < ny; i++)
        {
            in >> grid[GET_IDX_2D(j, i, nx, ny)];
            if (grid[GET_IDX_2D(j, i, nx, ny)] == 128 ||
                grid[GET_IDX_2D(j, i, nx, ny)] == 50)
                p++;
        }
    }
    in.close();

    
}

double Solver::gaussiana(double t)
{
    return (1. / (sigma * sqrt(2 * M_PI))) * exp(-0.5 * pow((t - T_peak) / sigma, 2));
}

double Solver::eq_terminal (Node *ptr, int k, int i, int j, double *reaction, int graph_id) {
    int id_viz = ptr->edges->id;
    double flow = calcula_flow(ptr);
    double radius = calcula_radius(ptr);
    double C[1];

    calc_concentracao_vizinhos(ptr, graph_id, C);

    #ifdef PRESSURE
    double term1 = -((flow / FACTOR_FLOW) * dt / dx) / (M_PI * pow(radius, 2));
    double term4 = -((flow / FACTOR_FLOW) * dt / dx) / (M_PI * pow(radius, 2));
    #else
    double term1 = -1. * dt / dx;
    #endif

    double *c_old = NULL;
    if (graph_id == 0) {
        c_old = this->concentration_old_network_1;
    }
    else if (graph_id == 1) {
        c_old = this->concentration_old_network_2;
    }
    else {
        c_old = this->concentration_old_network_3;
    }

    double difusao = D * (C[0] - c_old[k]) / (dx * dx);
    double adveccao = ((flow / FACTOR_FLOW) / (M_PI * pow(radius, 2))) * (c_old[k] - C[0]) / dx;
    double troca = reaction[GET_IDX_2D(i, j, n_tissue, n_tissue)] / dx;
    return (difusao - adveccao + troca) * dt + c_old[k];
}

double Solver::eq_geral(Node *ptr, int i, int graph_id) {
    double soma_difusao = calc_vizinhos_difusao(ptr, graph_id);
    double soma_adveccao = calc_vizinhos_adveccao(ptr, graph_id);
    double C[2];

    double flow = calcula_flow(ptr);
    double radius = calcula_radius(ptr);

    calc_concentracao_vizinhos(ptr, graph_id, C);

    #ifdef PRESSURE
    double *c_old = NULL;
    if (graph_id == 0) {
        c_old = this->concentration_old_network_1;
    }
    else if (graph_id == 1) {
        c_old = this->concentration_old_network_2;
    }
    else {
        c_old = this->concentration_old_network_3;
    }

    double adveccao = ((flow / FACTOR_FLOW) / (M_PI * pow(radius, 2))) * (c_old[i] - C[0]) / dx;
    double difusao = D * (C[1] - 2 * c_old[i] + C[0]) / (dx * dx);
    return (difusao - adveccao) * dt + c_old[i];
    #endif

}

double Solver::eq_bifurcacao(Node *ptr, int i, int graph_id) {
    double soma_difusao = calc_vizinhos_difusao(ptr, graph_id);
    double soma_adveccao = calc_vizinhos_adveccao(ptr, graph_id);

    double flow[3];
    double radius[3];
    double C[3];

    double area_media = calcula_area_bifurcacao(ptr);

    calc_concentracao_vizinhos(ptr, graph_id, C);
    calcula_elementos_bifurcacao(ptr, flow, radius);

    #ifdef PRESSURE
    double *c_old = NULL;
    if (graph_id == 0) {
        c_old = this->concentration_old_network_1;
    }
    else if (graph_id == 1) {
        c_old = this->concentration_old_network_2;
    }
    else {
        c_old = this->concentration_old_network_3;
    }

    double difusao = (D * (C[0] - (1 + pow(radius[1], 2) / pow(radius[0], 2) + pow(radius[2], 2) / pow(radius[0], 2)) * c_old[i] + (pow(radius[1], 2) / pow(radius[0], 2)) * C[1] + (pow(radius[2], 2) / pow(radius[0], 2)) * C[2])) / (dx * dx);
    double adveccao = ((c_old[i] / pow(radius[0], 2)) * (((flow[1] / FACTOR_FLOW) / (M_PI * pow(radius[1], 2))) * pow(radius[1], 2) + ((flow[2] / FACTOR_FLOW) / (M_PI * pow(radius[2], 2))) * pow(radius[2], 2)) - ((flow[0] / FACTOR_FLOW) / (M_PI * pow(radius[0], 2))) * C[0]) / dx;
    return (difusao - adveccao) * dt + c_old[i];
    #else
    double *c_old = this->concentration_old_network_1;
    return ((D * ((soma_difusao - (ptr->num_edges * c_old[i])) / (dx * dx))) - ((1. / dx) * (((ptr->num_edges - 1) * c_old[i]) - soma_adveccao))) * dt + c_old[i];
    #endif

}

double Solver::calc_vizinhos_difusao (Node *noAtual, int graph_id) {
    Edge *ptrl = noAtual->edges;
    double soma = 0.0;
    int id;
    double *c_old = NULL;

    if (graph_id == 0)
        c_old = this->concentration_new_network_1;
    else if (graph_id == 1)
        c_old = this->concentration_new_network_2;
    else if (graph_id == 2)
        c_old = this->concentration_new_network_3;

    while (ptrl != NULL){
        id = ptrl->id;
        soma += c_old[id];
        ptrl = ptrl->next;
    }

    return soma;
}

double Solver::calc_vizinhos_adveccao (Node *noAtual, int graph_id) {
    Edge *ptrl = noAtual->edges;
    double soma = 0.;
    int id;
    double *c_old = NULL;

    if (graph_id == 0)
        c_old = this->concentration_old_network_1;
    else if (graph_id == 1)
        c_old = this->concentration_old_network_2;
    else if (graph_id == 2)
        c_old = this->concentration_old_network_3;

    while (ptrl != NULL) {
        id = ptrl->id;
        if (ptrl->entrada == 1) {
            soma += c_old[id];
        }
        ptrl = ptrl->next;
    }

    return soma;
}

void Solver::calc_concentracao_vizinhos (Node *noAtual, int graph_id, double p[]) {
    Edge *ptrl = noAtual->edges;
    double soma = 0.;
    int id;
    int m = 1;
    double *c_old = NULL;

    if (graph_id == 1) {
        c_old = this->concentration_old_network_1;
    }
    else if (graph_id == 2) {
        c_old = this->concentration_old_network_2;
    }
    else {
        c_old = this->concentration_old_network_3;
    }

    while (ptrl != NULL) {
        id = ptrl->id;
        if (ptrl->entrada == 1) {
            p[0] = c_old[id];
        }
        else {
            p[m] = c_old[id];
            m++;
        }
        ptrl = ptrl->next;
    }
}

double Solver::calcula_area_bifurcacao(Node *noAtual)
{
    Edge *ptrl = noAtual->edges;
    double area[3];
    int i = 0;

    while (ptrl != NULL)
    {
        area[i] = M_PI * pow(ptrl->radius, 2);
        //         printf("%f \n",area[i]);
        i++;
        ptrl = ptrl->next;
    }

    //     int j;
    //     for(j = 0; j < 3; j++)
    //     {
    //         printf("%d \n",area[j]);
    //
    //     }

    //     printf("Erro no nó %d ao se calcular o radius da bifurcação\n", noAtual->id);
    return (area[0] + area[1] + area[2]) / 3.;
}

double Solver::calcula_flow (Node *noAtual) {
    Edge *ptrl = noAtual->edges;
    while (ptrl != NULL) {
        if (ptrl->entrada == 1)
            return ptrl->flow;

        ptrl = ptrl->next;
    }
    printf("[-] ERROR em nó %d ao se calcular o flow!\n", noAtual->id);
    return 0.0;
}

void Solver::calcula_elementos_bifurcacao(Node *noAtual, double flow[], double radius[])
{
    Edge *ptrl = noAtual->edges;
    int p = 1;

    while (ptrl != NULL)
    {
        if (ptrl->entrada == 1)
        {
            flow[0] = ptrl->flow;
            radius[0] = ptrl->radius;
        }
        else
        {
            flow[p] = ptrl->flow;
            radius[p] = ptrl->radius;
            //             printf("%d\n",p);
            p++;
        }
        //         p++;
        ptrl = ptrl->next;
    }
    //     printf("Erro no nó %d ao se calcular o flow ou o radius da bifurcação \n", noAtual->id);
}

// void Solver::calcula_elementos_bifurcacao (Node * noAtual, double flow[], double area[])
// {
//     Edge * ptrl = noAtual->edges;
//     int p = 1;
//
//     while(ptrl != NULL)
//     {
//         if (ptrl->entrada == 1)
//         {
//             flow[0] = ptrl->flow;
//             area[0] = M_PI*pow(ptrl->radius,2);
//         }
//         else
//         {
//             flow[p] = ptrl->flow;
//             area[p] = M_PI*pow(ptrl->radius,2);
// //             printf("%d\n",p);
//             p++;
//         }
// //         p++;
//         ptrl = ptrl->next;
//     }
// //     printf("Erro no nó %d ao se calcular o flow ou o radius da bifurcação \n", noAtual->id);
// }

double Solver::calcula_radius (Node *noAtual) {
    Edge *ptrl = noAtual->edges;
    while (ptrl != NULL) {
        if (ptrl->entrada == 1)
            return ptrl->radius;

        ptrl = ptrl->next;
    }
    printf("[-] ERROR em nó %d ao se calcular o radius\n", noAtual->id);
    return 0.0;
}

// double Solver::calcula_area_bifurcacao (Node * noAtual)
// {
//     Edge * ptrl = noAtual->edges;
//     double area[3];
//     int i = 0;
//
//     while(ptrl != NULL)
//     {
//         area[i] = M_PI * pow(ptrl->edges,2);
//         i++;
//         ptrl = ptrl->next;
//     }
//     printf("Erro no nó %d ao se calcular o radius da bifurcação\n", noAtual->id);
//     return (area[0]+area[1]+area[2])/3;
// }

// TODO: Mudar a copia para swap de ponteiro
void Solver::next_timestep() {
    int i, j;
    int np;

    np = g1->getTotalNodes();
    for (i = 0; i < np; i++)
        this->concentration_old_network_1[i] = this->concentration_new_network_1[i];

    np = g2->getTotalNodes();
    for (i = 0; i < np; i++)
        this->concentration_old_network_2[i] = this->concentration_new_network_2[i];

    np = g3->getTotalNodes();
    for (i = 0; i < np; i++)
        this->concentration_old_network_3[i] = this->concentration_new_network_3[i];

    for (i = 0; i < n_tissue; i++) {
        for (j = 0; j < n_tissue; j++) {
            this->concentration_old_tissue[GET_IDX_2D(i, j, n_tissue, n_tissue)] = this->concentration_new_tissue[GET_IDX_2D(i, j, n_tissue, n_tissue)];
            this->concentration_old_fibrosis[GET_IDX_2D(i, j, n_tissue, n_tissue)] = this->concentration_new_fibrosis[GET_IDX_2D(i, j, n_tissue, n_tissue)];
        }
    }
}

// double calcNorm (double x1, double y1, double z1, double x2, double y2, double z2)"
// {
//     return sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2));
// }

void Solver::print()
{
    //     int np = g1->getTotalNodes();
    //     for (int i = 0; i < np; i++)
    //         printf("id = %d -- C = %.10lf\n",i,cOld[i].concentration);
    Edge *ptrl;
    Node *ptr = g1->getListNodes();
    printf("======================= PRINTING GRAPH ================================\n");
    while (ptr != NULL)
    {
        //         #ifdef PRESSURE
        // 	    printf("|| %d (%d) (%.4lf %.4lf %.4lf) [%.4lf] %d ||",ptr->id,ptr->type,ptr->x,ptr->y,ptr->z,ptr->pressure,ptr->num_edges);
        // 		#else
        //         //printf("|| %d (%d) (%.4lf %.4lf %.4lf) %d ||",ptr->id,ptr->type,ptr->x,ptr->y,ptr->z,ptr->num_edges);
        //         printf("|| %d ||", ptr->id);
        // 		#endif

        ptrl = ptr->edges;
        while (ptrl != NULL)
        {
            // 			printf(" --> || %d %.4lf (%.4lf %.4lf %.4lf) ||",ptrl->id,ptrl->w,ptrl->dest->x,ptrl->dest->y,ptrl->dest->z);
            //             printf(" --> || %d (%d) ||",ptrl->id,ptrl->entrada);
            printf(" -> | Dist = %.4lf | Vel = %.4lf |\n", calc_norm(ptrl->dest->x, ptrl->dest->y, ptrl->dest->z, ptr->x, ptr->y, ptr->z), ptr->velocity);
            ptrl = ptrl->next;
        }
        printf("\n");
        ptr = ptr->next;
    }
    printf("=======================================================================\n");
}

// double calcNorm (double x1, double y1, double z1, double x2, double y2, double z2)
// {
//     return sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2));
// }

void Solver::allocate_memory_for_network_cells(int net_id, int num_cells)
{
    if (net_id == 1)
    {
        //this->cell_types_network_1 = (int *)malloc(sizeof(int) * num_cells);
        //this->roi_network_1 = (int *)malloc(sizeof(int) * num_cells);
        //this->concentration_old_network_1 = (double *)malloc(sizeof(double) * num_cells);
        //this->concentration_new_network_1 = (double *)malloc(sizeof(double) * num_cells);
        
        this->cell_types_network_1 = (int*)calloc(num_cells, sizeof(int));
        this->roi_network_1 = (int*)calloc(num_cells, sizeof(int));
        this->concentration_old_network_1 = (double*)calloc(num_cells, sizeof(double));
        this->concentration_new_network_1 = (double*)calloc(num_cells, sizeof(double));
    }
    else if (net_id == 2)
    {
        //this->cell_types_network_2 = (int *)malloc(sizeof(int) * num_cells);
        //this->roi_network_2 = (int *)malloc(sizeof(int) * num_cells);
        //this->concentration_old_network_2 = (double *)malloc(sizeof(double) * num_cells);
        //this->concentration_new_network_2 = (double *)malloc(sizeof(double) * num_cells);

        this->cell_types_network_2 = (int*)calloc(num_cells, sizeof(int));
        this->roi_network_2 = (int*)calloc(num_cells, sizeof(int));
        this->concentration_old_network_2 = (double*)calloc(num_cells, sizeof(double));
        this->concentration_new_network_2 = (double*)calloc(num_cells, sizeof(double));
    }
    else
    {
        //this->cell_types_network_3 = (int *)malloc(sizeof(int) * num_cells);
        //this->roi_network_3 = (int *)malloc(sizeof(int) * num_cells);
        //this->concentration_old_network_3 = (double *)malloc(sizeof(double) * num_cells);
        //this->concentration_new_network_3 = (double *)malloc(sizeof(double) * num_cells);

        this->cell_types_network_3 = (int*)calloc(num_cells, sizeof(int));
        this->roi_network_3 = (int*)calloc(num_cells, sizeof(int));
        this->concentration_old_network_3 = (double*)calloc(num_cells, sizeof(double));
        this->concentration_new_network_3 = (double*)calloc(num_cells, sizeof(double));
    }
}

void Solver::allocate_memory_for_tissue()
{
    this->grid = (int *)malloc(sizeof(int) * this->n_tissue * this->n_tissue);
    this->cell_types_tissue = (int *)malloc(sizeof(int) * this->n_tissue * this->n_tissue);
    this->concentration_old_tissue = (double *)malloc(sizeof(double) * this->n_tissue * this->n_tissue);
    this->concentration_new_tissue = (double *)malloc(sizeof(double) * this->n_tissue * this->n_tissue);
    this->concentration_old_fibrosis = (double *)malloc(sizeof(double) * this->n_tissue * this->n_tissue);
    this->concentration_new_fibrosis = (double *)malloc(sizeof(double) * this->n_tissue * this->n_tissue);
}

void Solver::set_type_cells_using_network(Graph *g, int *cell_type, double *c_old, double *c_new)
{
    int i = 0;
    Node *ptr = g->getListNodes();

    while (ptr != NULL)
    {
        // Normal
        if (ptr->num_edges == 2)
        {
            cell_type[i] = 2;
        }
        // Bifurcação
        else if (ptr->num_edges > 2)
        {
            cell_type[i] = 3;
        }
        // Raiz (Pode ser que tenha que alterar depois)
        else if (ptr->num_edges == 1 && ptr->id == 0)
        {
            cell_type[i] = 0;
        }
        // Terminal
        else if (ptr->num_edges == 1)
        {
            cell_type[i] = 1;
        }

        // Inicialmente as concentracoes sao nulas
        c_old[i] = 0.0;
        c_new[i] = 0.0;

        // Avancar o ponteiro do grafo e o indice do vetor
        ptr = ptr->next;
        i++;
    }
}

void Solver::set_type_cells_tissue_using_network(int net_id)
{
    int np, p[2];
    int *cell_types_network = NULL;
    int *cell_types_tissue = this->cell_types_tissue;

    if (net_id == 1)
    {
        np = this->g1->getTotalNodes();
        cell_types_network = this->cell_types_network_1;
        grid_to_net_1.assign(np, std::vector<int>());
    }
    else if (net_id == 2)
    {
        np = this->g2->getTotalNodes();
        cell_types_network = this->cell_types_network_2;
        grid_to_net_2.assign(np, std::vector<int>());
    }
    else
    {
        np = this->g3->getTotalNodes();
        cell_types_network = this->cell_types_network_3;
        grid_to_net_3.assign(np, std::vector<int>());
    }

    // TODO: Change this in the refined version ... Load the connections from a file instead!
    // Tag the closest tissue cell to each terminal from the networks
    for (int k = 0; k < np; k++)
    {
        if (cell_types_network[k] == 1)
        {
            busca_proximal(k, p, net_id - 1);
            cell_types_tissue[GET_IDX_2D(p[0], p[1], n_tissue, n_tissue)] = 1; // Terminal
            if (net_id == 1) {
                grid_to_net_1[k].push_back(p[0]);
                grid_to_net_1[k].push_back(p[1]);
            }
            else if (net_id == 2) {
                grid_to_net_2[k].push_back(p[0]);
                grid_to_net_2[k].push_back(p[1]);
            }
            else if (net_id == 3) {
                grid_to_net_3[k].push_back(p[0]);
                grid_to_net_3[k].push_back(p[1]);
            }
        }
    }
}

void Solver::write_vtk_file_tissue(std::string filename, int nx, int ny)
{

    std::ofstream arqvtk;
    int *cell_types_tissue = this->cell_types_tissue;
    double *c_old_tissue = this->concentration_old_tissue;
    double *c_old_fibrosis = this->concentration_old_fibrosis;

    arqvtk.open(filename.c_str());
    arqvtk << "# vtk DataFile Version 3.0\n";
    arqvtk << "vtk output\n";
    arqvtk << "ASCII\n";
    arqvtk << "DATASET RECTILINEAR_GRID\n";
    arqvtk << "DIMENSIONS " << nx << " " << ny << " 1\n";

    arqvtk << "X_COORDINATES " << nx << " double \n";
    for (int i = 0; i < nx; i++)
        arqvtk << i * h_gadolinio << " ";
    arqvtk << "\n";

    arqvtk << "Y_COORDINATES " << ny << " double \n";
    for (int j = 0; j < ny; j++)
        arqvtk << j * h_gadolinio << " ";
    arqvtk << "\n";

    arqvtk << "Z_COORDINATES 1 double \n";
    arqvtk << "0\n";

    arqvtk << "POINT_DATA " << nx * ny << "\n";
    arqvtk << "FIELD FieldData 3\n";
    arqvtk << "Concentration 1 " << nx * ny << " double\n";
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            arqvtk << c_old_tissue[GET_IDX_2D(j, i, nx, ny)] + c_old_fibrosis[GET_IDX_2D(j, i, nx, ny)] << " ";
            //arqvtk << c_old_tissue[GET_IDX_2D(j, i, nx, ny)] << " ";
        }
    }
    arqvtk << "\nMETADATA\n";
    arqvtk << "INFORMATION 0\n\n";

    arqvtk << "CellType 1 " << nx * ny << " double\n";
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            arqvtk << cell_types_tissue[GET_IDX_2D(j, i, nx, ny)] << " ";
        }
    }
    arqvtk << "\nMETADATA\n";
    arqvtk << "INFORMATION 0\n\n";

    arqvtk << "Grid 1 " << nx * ny << " double\n";
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            arqvtk << grid[GET_IDX_2D(j, i, nx, ny)] << " ";
        }
    }
    arqvtk << "\nMETADATA\n";
    arqvtk << "INFORMATION 0\n\n";

    arqvtk << "\n";
    arqvtk.close();
}

int Solver::count_number_active_tissue_cells()
{
    int i, j;
    int N = 0;

    for (i = 0; i < this->n_tissue; i++)
    {
        for (j = 0; j < this->n_tissue; j++)
        {
            if (grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 128 ||
                grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 50)
                N++;
        }
    }
    return N;
}

void Solver::set_initial_conditions_tissue(double *reaction, double *peso)
{

    double *c_old_tissue = this->concentration_old_tissue;
    double *c_new_tissue = this->concentration_new_tissue;
    double *c_old_fibrosis = this->concentration_old_fibrosis;
    double *c_new_fibrosis = this->concentration_new_fibrosis;

    for (int i = 0; i < this->n_tissue; i++)
    {
        for (int j = 0; j < this->n_tissue; j++)
        {
            if (grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 255 ||
                grid[GET_IDX_2D(i, j, n_tissue, n_tissue)] == 0)
            {

                c_old_tissue[GET_IDX_2D(i, j, n_tissue, n_tissue)] = -0.0;
                c_new_tissue[GET_IDX_2D(i, j, n_tissue, n_tissue)] = -0.0;
                c_old_fibrosis[GET_IDX_2D(i, j, n_tissue, n_tissue)] = 0.0;
                c_new_fibrosis[GET_IDX_2D(i, j, n_tissue, n_tissue)] = 0.0;
            }
            else
            {
                c_old_tissue[GET_IDX_2D(i, j, n_tissue, n_tissue)] = 0.0;
                c_new_tissue[GET_IDX_2D(i, j, n_tissue, n_tissue)] = 0.0;
                c_old_fibrosis[GET_IDX_2D(i, j, n_tissue, n_tissue)] = 0.0;
                c_new_fibrosis[GET_IDX_2D(i, j, n_tissue, n_tissue)] = 0.0;
            }
            reaction[GET_IDX_2D(i, j, n_tissue, n_tissue)] = 0.0;
            peso[GET_IDX_2D(i, j, n_tissue, n_tissue)] = 0.0;
        }
    }
}

double Solver::integral_intra_v2 () {
    int j, np;

    int *cell_types;
    double *c_new;

    double V = 0.0;
    int n_terminais = 0;

    // Network 1
    np = g1->getTotalNodes();
    cell_types = this->cell_types_network_1;
    c_new = this->concentration_new_network_1;
    for (int j = 0; j < np; j++) {
        // Terminal
        if (cell_types[j] == 1) {
            V += c_new[j];
            n_terminais++;
        }
    }

    // Network 2
    np = g2->getTotalNodes();
    cell_types = this->cell_types_network_2;
    c_new = this->concentration_new_network_2;
    for (int j = 0; j < np; j++) {
        // Terminal
        if (cell_types[j] == 1) {
            V += c_new[j];
            n_terminais++;
        }
    }

    // Network 3
    np = g3->getTotalNodes();
    cell_types = this->cell_types_network_3;
    c_new = this->concentration_new_network_3;
    for (int j = 0; j < np; j++) {
        // Terminal
        if (cell_types[j] == 1) {
            V += c_new[j];
            n_terminais++;
        }
    }

    // Mean
    return V / n_terminais;
}

void Solver::solve1D_source_v2 (double *C_1D, double V, double adv_1D, double difusao_1D) {
    int i; 
    double dif_1D_global = 0.0;

    //double *new_C_1D = new double[n_1D];
    
    //C_1D[0] = V;
    this->C_1D_old[0] = V;
    for (i = 1; i < n_1D; i++) {
        if (i != n_1D - 1) {
            //adv_1D = (v_1D * C_1D[i] - v_1D * C_1D[i - 1]) / h_1D;
            //difusao_1D = (dif_1D_global * (C_1D[i + 1] - C_1D[i]) - dif_1D_global * (C_1D[i] - C_1D[i - 1])) / (h_1D * h_1D);
            adv_1D = (v_1D * this->C_1D_old[i] - v_1D * this->C_1D_old[i - 1]) / h_1D;
            difusao_1D = (dif_1D_global * (this->C_1D_old[i + 1] - this->C_1D_old[i]) - dif_1D_global * (this->C_1D_old[i] - this->C_1D_old[i - 1])) / (h_1D * h_1D);
        }
        else {
            //adv_1D = (v_1D * C_1D[i] - v_1D * C_1D[i - 1]) / h_1D;
            //difusao_1D = 0.;
            adv_1D = (v_1D * this->C_1D_old[i] - v_1D * this->C_1D_old[i - 1]) / h_1D;
            difusao_1D = 0.;
        }

        if (i == n_1D - 1)
            //new_C_1D[i] = (-adv_1D - rins * C_1D[i] - (dif_1D_global * (C_1D[i] - C_1D[i - 1])) / (h_1D * h_1D)) * dt / (phi_1D) + C_1D[i];
            this->C_1D_new[i] = (-adv_1D - rins * this->C_1D_old[i] - (dif_1D_global * (this->C_1D_old[i] - this->C_1D_old[i - 1])) / (h_1D * h_1D)) * dt / (phi_1D) + this->C_1D_old[i];
        else
            //new_C_1D[i] = (difusao_1D - adv_1D - rins * C_1D[i]) * dt / (phi_1D) + C_1D[i];
            this->C_1D_new[i] = (difusao_1D - adv_1D - rins * this->C_1D_old[i]) * dt / (phi_1D) + this->C_1D_old[i];
    }

    for (i = 1; i < n_1D; i++) {
        //C_1D[i] = new_C_1D[i];
        this->C_1D_old[i] = this->C_1D_new[i];
    }

    //delete [] new_C_1D;
}

void Solver::write_vtk_file_networks (Graph *g, int *cell_type, int *roi, double *c_new, int net_id, int iter) {
    FILE *file;
    char filename[500];

    Node *ptr = g->getListNodes();
    int np = g->getTotalNodes();
    int ne = g->getTotalEdges();

    // Escrever o cabecalho do VTK
    sprintf(filename, "%s/graph%d/sol%d.vtk", this->output_dir.c_str(), net_id+1, iter);
    file = fopen(filename, "w+");
    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Equacao Transporte\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET POLYDATA\n");
    fprintf(file, "POINTS %d float\n", np);

    // Escrever pontos
    while (ptr != NULL) {
        fprintf(file, "%.10lf %.10lf %.10lf\n", ptr->x, ptr->y, ptr->z);
        ptr = ptr->next;
    }

    // Escrever linhas
    ptr = g->getListNodes();

    fprintf(file, "LINES %d %d\n", ne, ne * 3);
    while (ptr != NULL) {
        Edge *ptrl = ptr->edges;
        while (ptrl != NULL) {
            fprintf(file, "2 %d %d\n", ptr->id, ptrl->dest->id);
            ptrl = ptrl->next;
        }
        ptr = ptr->next;
    }

    fprintf(file, "POINT_DATA %d\n", np);
    fprintf(file, "SCALARS concentracao float 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    // Escrever concentracoes
    ptr = g->getListNodes();
    while (ptr != NULL) {
        fprintf(file, "%.10lf\n", c_new[ptr->id]);
        ptr = ptr->next;
    }

    fprintf(file, "SCALARS celltype float 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    // Escrever tipos de celulas
    ptr = g->getListNodes();
    while (ptr != NULL) {
        fprintf(file, "%d\n", cell_type[ptr->id]);
        ptr = ptr->next;
    }

    fprintf(file, "SCALARS roi float 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    // Escrever tipos de celulas
    ptr = g->getListNodes();
    while (ptr != NULL) {
        fprintf(file, "%d\n", roi[ptr->id]);
        ptr = ptr->next;
    }

    // Escrever os raios
    ptr = g->getListNodes();
    fprintf(file, "CELL_DATA %d\n", ne);
    fprintf(file, "scalars radius float\nLOOKUP_TABLE default\n");

    while (ptr != NULL) {
        Edge *ptrl = ptr->edges;
        while (ptrl != NULL) {
            fprintf(file, "%.10lf\n", ptrl->radius);
            ptrl = ptrl->next;
        }
        ptr = ptr->next;
    }

    fclose(file);
}

void Solver::solve_concentration_for_network(Graph *g, int net_id, double *c_old, double *c_new, \
                                                       int *cell_types, int *roi, \
                                                       double *reaction, double fonte_new, \
                                                       double &SI_health_intra) {
    int u[2];
    int np = g->getTotalNodes();
    Node *ptr = g->getListNodes();
    for (int j = 0; j < np; j++) {
        // Raiz
        if (cell_types[j] == 0) {
            c_new[j] = fonte_new;
        }
        // Terminal
        else if (cell_types[j] == 1) {
            //busca_proximal(j, u, net_id);
            search_proximal_in_map(j, u, net_id);
            calcula_reaction(reaction, j, u[0], u[1], net_id);
            c_new[j] = eq_terminal(ptr, j, u[0], u[1], reaction, net_id);
        }
        // Normal
        else if (cell_types[j] == 2) {
            c_new[j] = eq_geral(ptr, j, net_id);
        }
        // Bifurcacao
        else if (cell_types[j] == 3) 
            c_new[j] = eq_bifurcacao(ptr, j, net_id);

        // Healthy
        if (roi[j] == 0)
            SI_health_intra += c_new[j];

        ptr = ptr->next;
    }

// DEBUG
/*
    if (net_id == 2) {
        char filename[500];
        sprintf(filename, "outputs/sol_Lucas_graph3_iter_%d.txt", counter);
        //sprintf(filename, "outputs/reaction_Lucas_tiss_iter_%d.txt", counter);
        FILE *file = fopen(filename, "w+");
        for (int j = 0; j < np; j++) {
            fprintf(file, "%g\n", c_new[j]);
        }
        //for (int j = 0; j < n_tissue; j++) {
        //    for (int i = 0; i < n_tissue; i++) {
        //        fprintf(file, "%g\n", reaction[GET_IDX_2D(j, i, n_tissue, n_tissue)]);
        //    }
        //}
        fclose(file);
        
        // DEBUG
        if (counter == 10) {
            exit(1);
        }
        else {
            counter++;
        }
    }
*/
}

void Solver::solve_concentration_for_tissue (double *reaction) {
    
    int p, q;
    double dif_tissue;
    double *c_old_tissue = this->concentration_old_tissue;
    double *c_old_tissue_fibrosis = this->concentration_old_fibrosis;
    double *c_new_tissue = this->concentration_new_tissue;
    double *c_new_tissue_fibrosis = this->concentration_new_fibrosis;
    int *cell_types_tissue = this->cell_types_tissue;
    for (p = 0; p < this->n_tissue; p++) {
        for (q = 0; q < this->n_tissue; q++) {
            // Check for active tissue cells
            if (grid[GET_IDX_2D(p, q, n_tissue, n_tissue)] == 128 || \
                grid[GET_IDX_2D(p, q, n_tissue, n_tissue)] == 50) {
                    // Case 1
                    if ((grid[GET_IDX_2D(p + 1, q, n_tissue, n_tissue)] == 255 && grid[GET_IDX_2D(p, q + 1, n_tissue, n_tissue)] == 255) || \
                        (grid[GET_IDX_2D(p + 1, q, n_tissue, n_tissue)] == 0 && grid[GET_IDX_2D(p, q + 1, n_tissue, n_tissue)] == 0))
                            dif_tissue = (D / (h_gadolinio * h_gadolinio)) * \
                                         (-(c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q - 1, n_tissue, n_tissue)]) - (c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p - 1, q, n_tissue, n_tissue)]));
                    // Case 2
                    else if ((grid[GET_IDX_2D(p + 1, q, n_tissue, n_tissue)] == 255 && grid[GET_IDX_2D(p, q - 1, n_tissue, n_tissue)] == 255) || \
                            (grid[GET_IDX_2D(p + 1, q, n_tissue, n_tissue)] == 0 && grid[GET_IDX_2D(p, q - 1, n_tissue, n_tissue)] == 0))
                            dif_tissue = (D / (h_gadolinio * h_gadolinio)) * \
                                    (+(c_old_tissue[GET_IDX_2D(p, q + 1, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)]) - (c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p - 1, q, n_tissue, n_tissue)]));
                    // Case 3
                    else if ((grid[GET_IDX_2D(p - 1, q, n_tissue, n_tissue)] == 255 && grid[GET_IDX_2D(p, q + 1, n_tissue, n_tissue)] == 255) ||
                            (grid[GET_IDX_2D(p - 1, q, n_tissue, n_tissue)] == 0 && grid[GET_IDX_2D(p, q + 1, n_tissue, n_tissue)] == 0))
                            dif_tissue = (D / (h_gadolinio * h_gadolinio)) * \
                                    (-(c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q - 1, n_tissue, n_tissue)]) + (c_old_tissue[GET_IDX_2D(p + 1, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)]));
                    // Case 4
                    else if ((grid[GET_IDX_2D(p - 1, q, n_tissue, n_tissue)] == 255 && grid[GET_IDX_2D(p, q - 1, n_tissue, n_tissue)] == 255) ||
                            (grid[GET_IDX_2D(p - 1, q, n_tissue, n_tissue)] == 0 && grid[GET_IDX_2D(p, q - 1, n_tissue, n_tissue)] == 0))
                            dif_tissue = (D / (h_gadolinio * h_gadolinio)) * \
                                    (+(c_old_tissue[GET_IDX_2D(p, q + 1, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)]) + (c_old_tissue[GET_IDX_2D(p + 1, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)]));
                    // Case 5
                    else if ((grid[GET_IDX_2D(p, q + 1, n_tissue, n_tissue)] == 255) || (grid[GET_IDX_2D(p, q + 1, n_tissue, n_tissue)] == 0))
                            dif_tissue = (D / (h_gadolinio * h_gadolinio)) * \
                                    (-(c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q - 1, n_tissue, n_tissue)]) + (c_old_tissue[GET_IDX_2D(p + 1, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)]) - (c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p - 1, q, n_tissue, n_tissue)]));
                    // Case 6
                    else if ((grid[GET_IDX_2D(p, q - 1, n_tissue, n_tissue)] == 255) || (grid[GET_IDX_2D(p, q - 1, n_tissue, n_tissue)] == 0))
                            dif_tissue = (D / (h_gadolinio * h_gadolinio)) * \
                                    (+(c_old_tissue[GET_IDX_2D(p, q + 1, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)]) + (c_old_tissue[GET_IDX_2D(p + 1, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)]) - (c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p - 1, q, n_tissue, n_tissue)]));
                    // Case 7
                    else if ((grid[GET_IDX_2D(p + 1, q, n_tissue, n_tissue)] == 255) || (grid[GET_IDX_2D(p + 1, q, n_tissue, n_tissue)] == 0))
                            dif_tissue = (D / (h_gadolinio * h_gadolinio)) * \
                                    (+(c_old_tissue[GET_IDX_2D(p, q + 1, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)]) - (c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q - 1, n_tissue, n_tissue)]) - (c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p - 1, q, n_tissue, n_tissue)]));
                    // Case 8
                    else if ((grid[GET_IDX_2D(p - 1, q, n_tissue, n_tissue)] == 255) || (grid[GET_IDX_2D(p - 1, q, n_tissue, n_tissue)] == 0))
                            dif_tissue = (D / (h_gadolinio * h_gadolinio)) * \
                                    (+(c_old_tissue[GET_IDX_2D(p, q + 1, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)]) - (c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q - 1, n_tissue, n_tissue)]) + (c_old_tissue[GET_IDX_2D(p + 1, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)]));
                    // Case 9
                    else
                            dif_tissue = (D / (h_gadolinio * h_gadolinio)) *
                                    (+(c_old_tissue[GET_IDX_2D(p, q + 1, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)]) - (c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q - 1, n_tissue, n_tissue)]) + (c_old_tissue[GET_IDX_2D(p + 1, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)]) - (c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] - c_old_tissue[GET_IDX_2D(p - 1, q, n_tissue, n_tissue)]));


                // Not a Terminal
                if (cell_types_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] != 1) 
                    reaction[GET_IDX_2D(p, q, n_tissue, n_tissue)] = 0.0;

                if ((q > n_tissue / 2. && q < 74. * n_tissue / 100. && p > 43. * n_tissue / 100. && p < 57. * n_tissue / 100.) || \
                    (q == 74. * n_tissue / 100. && p > 44. * n_tissue / 100. && p < 56. * n_tissue / 100.) || \
                    (q == 75. * n_tissue / 100. && p > 45. * n_tissue / 100. && p < 55. * n_tissue / 100.) || \
                    (q == 76. * n_tissue / 100. && p > 46. * n_tissue / 100. && p < 54. * n_tissue / 100.) || \
                    (q == 77. * n_tissue / 100. && p > 48. * n_tissue / 100. && p < 54. * n_tissue / 100.)) {
                        c_new_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] =
                            (-lambda_fibrose * rf * c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] + dif_tissue - reaction[GET_IDX_2D(p, q, n_tissue, n_tissue)] / ((1.0 - phi) * lambda_infarto) - decaimento_infarto * c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] / lambda_infarto) * dt \
                            + c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)];

                        if (lambda_fibrose == 0.0)
                            c_new_tissue_fibrosis[GET_IDX_2D(p, q, n_tissue, n_tissue)] = 0.0;
                        else
                            c_new_tissue_fibrosis[GET_IDX_2D(p, q, n_tissue, n_tissue)] =
                                (rf * c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] - decaimento_fibrose * c_old_tissue_fibrosis[GET_IDX_2D(p, q, n_tissue, n_tissue)] / lambda_infarto) * dt + c_old_tissue_fibrosis[GET_IDX_2D(p, q, n_tissue, n_tissue)];
                }
                else {
                    c_new_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] =
                        (dif_tissue - reaction[GET_IDX_2D(p, q, n_tissue, n_tissue)] / ((1.0 - phi) * lambda) - decaimento * c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)] / lambda) * dt \
                        + c_old_tissue[GET_IDX_2D(p, q, n_tissue, n_tissue)];
                }
            }
        }
    }
}

void Solver::calculate_signal_intensity (int i, double *peso, double SI_ischemic, double &SI_health_tissue, double &SI_fibrose, double &SI, \
                                         double *SI_ischemic_vec, double *SI_fibrose_vec, double *SI_health_vec) {
    if (i % SAVEEACH == 0) {
        SI_health_tissue = signal_intensity_healthy_v2(i, grid, this->concentration_new_tissue, peso, SI_health_tissue); 
        SI_fibrose = signal_intensity_fibrose_v2(i, grid, this->concentration_new_tissue, this->concentration_new_fibrosis, peso, SI_fibrose, SI_ischemic, SI);
        SI = signal_intensity_v2(i, grid, this->concentration_new_tissue, this->concentration_new_fibrosis, peso, SI_fibrose, SI_ischemic, SI);

        SI_ischemic_vec[i / SAVEEACH] = SI;
        SI_fibrose_vec[i / SAVEEACH] = SI_fibrose;
        SI_health_vec[i / SAVEEACH] = SI_health_tissue;
    }
}

void Solver::post_processing (double *SI_health_vec, double *SI_ischemic_vec, double *SI_fibrose_vec) {
    FILE *gnuplot;
    char filename[500];
    double aux_1 = 0.0;
    double aux_2 = 0.0;

    for (int i = 0; i < M; i++) {
        if (i % SAVEEACH == 0 || i == 0) {
            aux_1 = SI_health_vec[i / SAVEEACH];
            
            if (aux_1 > aux_2)
                aux_2 = aux_1;
        }
    }

    double fator = 83.63226015538544 / aux_2;
    //printf("\n \n %f %f \n \n", aux_2, fator);

    for (int i = 0; i < M; i++) {
        if (i % SAVEEACH == 0 || i == 0) {
            sprintf(filename, "%s/gauss.txt", this->output_dir.c_str());
            gnuplot = fopen(filename, "a");
            //             fprintf(gnuplot,"%f %f %f %f\n", i*dt, fator_SI*(phi*SI_health_intra + SI_health_tissue), fator_SI*phi*SI_health_intra, fator_SI*SI_health_tissue);
            //             fprintf(gnuplot,"%f %f %f \n", i*dt, fator_SI*SI_health_tissue, fator_SI*SI_ischemic_tissue);
            fprintf(gnuplot, "%.4f %.4f %.4f %.4f\n", i * dt, fator * SI_health_vec[i / SAVEEACH], fator * SI_ischemic_vec[i / SAVEEACH], fator * SI_fibrose_vec[i / SAVEEACH]);
        }
    }
}

void Solver::search_proximal_in_map (int j, int u[], int net_id) {
    if (net_id == 0 && grid_to_net_1[j].size() == 2) {
        u[0] = grid_to_net_1[j][0];
        u[1] = grid_to_net_1[j][1];
    }
    else if (net_id == 1 && grid_to_net_2[j].size() == 2) {
        u[0] = grid_to_net_2[j][0];
        u[1] = grid_to_net_2[j][1];
    }
    else if (net_id == 2 && grid_to_net_3[j].size() == 2) {
        u[0] = grid_to_net_3[j][0];
        u[1] = grid_to_net_3[j][1];
    }
}

void Solver::create_folders () {
    char folder_path[500];
    printf("[!] Simulation files will be saved at :> '%s'\n", this->output_dir.c_str());
    if (mkdir(this->output_dir.c_str(), 0777) == -1) {
        printf("\t[-] ERROR! Creating folder '%s'! Probably it already exists!\n", this->output_dir.c_str());
    }
    else {
        printf("\t[+] Folder '%s' created!\n", this->output_dir.c_str());
    }

    sprintf(folder_path, "%s/tissue", this->output_dir.c_str());
    if (mkdir(folder_path, 0777) == -1) {
        printf("\t[-] ERROR! Creating folder '%s'! Probably it already exists!\n", folder_path);
    }
    else {
        printf("\t[+] Folder '%s' created!\n", folder_path);
    }

    for (int i = 0; i < 3; i++) {
        sprintf(folder_path, "%s/graph%d", this->output_dir.c_str(), i+1);
        if (mkdir(folder_path, 0777) == -1) {
            printf("\t[-] ERROR! Creating folder '%s'! Probably it already exists!\n", folder_path);
        }
        else {
            printf("\t[+] Folder '%s' created!\n", folder_path);
        }
    }
}