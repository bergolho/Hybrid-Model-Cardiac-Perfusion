#include "../include/solver_joao.h"

// Number of threads to solve the system of ODEs
//static constexpr int nthreads = 4;

// TODO: Try to put this inside the Solver_Joao class
int joao_counter = 0;
std::vector<std::vector<int>> joao_grid_to_net_1;
std::vector<std::vector<int>> joao_grid_to_net_2;
std::vector<std::vector<int>> joao_grid_to_net_3;

Solver_Joao::Solver_Joao (int argc, char *argv[], int load_map)
{
    // Parse dos argumentos de entrada
    dt = atof(argv[1]);
    tmax = atof(argv[2]);
    network_filename_1 = argv[3];
    network_filename_2 = argv[4];
    network_filename_3 = argv[5];
    tissue_filename = argv[6];
    output_dir = argv[7];
    
    createFolders();

    // Numero de passos de tempo
    M = nearbyint(tmax/dt);
    // Cria o grafo a partir do arquivo de malha
    g = new Graph(network_filename_1, dx);
    g2 = new Graph(network_filename_2, dx);
    g3 = new Graph(network_filename_3, dx);

    leGeometria(tissue_filename);

    // Alocando memoria para as celulas
    setControlVolumes();
    
    // Inicilizar os tipos de celulas
    setTypeCell();
    // OK
     
    setTypeTissue();
    //setTypeTissue_v2();        
    
    // TO DO
    //setPlot();
}

void Solver_Joao::setTypeCell ()
{
    Node *ptr;
    int i;
    
    // Grafo 1
    ptr = g->getListNodes();
    i = 0;
    while (ptr != NULL)
    {
        // Tipo Normal
        if (ptr->num_edges == 2)
        {
            cOld[i].type = 2; 
            cNew[i].type = 2;
        }
        // Tipo bifurcação
        else if (ptr->num_edges > 2)
        {
            cOld[i].type = 3; 
            cNew[i].type = 3;
        }
        // Tipo raiz (Pode ser que tenha que alterar depois)
        else if (ptr->num_edges == 1 && ptr->id == 0)
        {
            cOld[i].type = 0;
            cNew[i].type = 0;
//             printf("cu\n");
        }
        // Tipo terminal
        else if (ptr->num_edges == 1)
        {
            cOld[i].type = 1;
            cNew[i].type = 1;
        }
        
        // Inicialmente as concentracoes sao nulas
        cOld[i].concentration = 0;
        cNew[i].concentration = 0;
        
        // Avancar o ponteiro do grafo e o indice do vetor
        ptr = ptr->next;
        i++;
    }
    
    // Grafo 2
    ptr = g2->getListNodes();
    i = 0;
    while (ptr != NULL)
    {
        // Tipo Normal
        if (ptr->num_edges == 2)
        {
            cOld2[i].type = 2; 
            cNew2[i].type = 2;
        }
        // Tipo bifurcação
        else if (ptr->num_edges > 2)
        {
            cOld2[i].type = 3; 
            cNew2[i].type = 3;
        }
        // Tipo raiz (Pode ser que tenha que alterar depois)
        else if (ptr->num_edges == 1 && ptr->id == 0)
        {
            cOld2[i].type = 0;
            cNew2[i].type = 0;
//             printf("cu2\n");
        }
        // Tipo terminal
        else if (ptr->num_edges == 1)
        {
            cOld2[i].type = 1;
            cNew2[i].type = 1;
        }
        
        // Inicialmente as concentracoes sao nulas
        cOld2[i].concentration = 0;
        cNew2[i].concentration = 0;
        
        // Avancar o ponteiro do grafo e o indice do vetor
        ptr = ptr->next;
        i++;
    }
    
    // Grafo 3
    ptr = g3->getListNodes();
    i = 0;
    while (ptr != NULL)
    {
        // Tipo Normal
        if (ptr->num_edges == 2)
        {
            cOld3[i].type = 2; 
            cNew3[i].type = 2;
        }
        // Tipo bifurcação
        else if (ptr->num_edges > 2)
        {
            cOld3[i].type = 3; 
            cNew3[i].type = 3;
        }
        // Tipo raiz (Pode ser que tenha que alterar depois)
        else if (ptr->num_edges == 1 && ptr->id == 0)
        {
            cOld3[i].type = 0;
            cNew3[i].type = 0;
//             printf("cu3\n");
        }
        // Tipo terminal
        else if (ptr->num_edges == 1)
        {
            cOld3[i].type = 1;
            cNew3[i].type = 1;
        }
        
        // Inicialmente as concentracoes sao nulas
        cOld3[i].concentration = 0;
        cNew3[i].concentration = 0;
        
        // Avancar o ponteiro do grafo e o indice do vetor
        ptr = ptr->next;
        i++;
    }
}

void Solver_Joao::setControlVolumes ()
{
    // Capturar o numero de volumes do primeiro grafo
    int np = g->getTotalNodes();
    
    cOld = (Cell_Joao*)malloc(sizeof(Cell_Joao)*np);
    cNew = (Cell_Joao*)malloc(sizeof(Cell_Joao)*np);
    
    np = g2->getTotalNodes();
    
    cOld2 = (Cell_Joao*)malloc(sizeof(Cell_Joao)*np);
    cNew2 = (Cell_Joao*)malloc(sizeof(Cell_Joao)*np);
    
    np = g3->getTotalNodes();
    
    cOld3 = (Cell_Joao*)malloc(sizeof(Cell_Joao)*np);
    cNew3 = (Cell_Joao*)malloc(sizeof(Cell_Joao)*np);
    
    cOld_tissue = (Cell_Joao*)malloc(sizeof(Cell_Joao)*n_tissue*n_tissue);
    cNew_tissue = (Cell_Joao*)malloc(sizeof(Cell_Joao)*n_tissue*n_tissue);
}

void Solver_Joao::writeVTKFile_tissue(char * nome, int nx, int ny)
{
    std::ofstream arqvtk;
    arqvtk.open(nome);
    arqvtk << "# vtk DataFile Version 3.0\n";
    arqvtk << "vtk output\n";
    arqvtk << "ASCII\n";
    arqvtk << "DATASET RECTILINEAR_GRID\n";
    arqvtk << "DIMENSIONS " << nx << " " << ny << " 1\n";

    arqvtk << "X_COORDINATES " << nx << " double \n";
    for(int i=0; i<nx; i++)
    arqvtk << i*h_gadolinio << " ";
    arqvtk << "\n";

    arqvtk << "Y_COORDINATES " << ny << " double \n";
    for(int j=0; j<ny; j++)
    arqvtk << j*h_gadolinio << " ";
    arqvtk << "\n";

    arqvtk << "Z_COORDINATES 1 double \n";
    arqvtk << "0\n";

    arqvtk << "POINT_DATA " << nx*ny << "\n";
    arqvtk << "FIELD FieldData 3\n";
//   arqvtk << "Pressure(kPa) 1 " << nx*ny << " double \n";
    arqvtk << "C 1 " << nx*ny << " double \n";

    for(int j=0; j<ny; j++)
        for(int i=0; i<nx; i++)
            {
                arqvtk << cOld_tissue[GET_IDX_2D(j,i,nx,ny)].concentration << " ";
            }

    arqvtk << "\nMETADATA\n";
    arqvtk << "INFORMATION 0\n\n";

    arqvtk << "CellType 1 " << nx * ny << " double\n";
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            arqvtk << cOld_tissue[GET_IDX_2D(j, i, nx, ny)].type << " ";
        }
    }
    arqvtk << "\nMETADATA\n";
    arqvtk << "INFORMATION 0\n\n";

    arqvtk << "Grid 1 " << nx * ny << " double\n";
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            arqvtk << tecido[GET_IDX_2D(j, i, nx, ny)] << " ";
        }
    }
    arqvtk << "\nMETADATA\n";
    arqvtk << "INFORMATION 0\n\n";

//    arqvtk << "VECTORS vectors double \n";
//    for(int i = 0; i < nx; i++)
//        for(int j=0; j<ny; j++)
//            {
//                arqvtk << vxGET_IDX_2D(i,j,nx,ny) << " "<< vyGET_IDX_2D(i,j,nx,ny) << " "<< 0.0 << "    ";
//            }

    arqvtk << "\n";
    arqvtk.close();
}

void Solver_Joao::writeVTKFile (int iter, int graph_id)
{
    FILE *file;
    int np, ne;
    char filename[50];
    Node *ptr;
    
    if (graph_id == 0)
    {
        ptr = g->getListNodes();
        np = g->getTotalNodes();
        ne = g->getTotalEdges();
    }
    else if (graph_id == 1)
    {
        ptr = g2->getListNodes();
        np = g2->getTotalNodes();
        ne = g2->getTotalEdges();
    }
    else if (graph_id == 2)
    {
        ptr = g3->getListNodes();
        np = g3->getTotalNodes();
        ne = g3->getTotalEdges();
    }

    // Escrever o cabecalho do VTK
//     sprintf(filename,"VTK/sol%d.vtk",iter);
    sprintf(filename,"%s/graph%d/sol%d.vtk",this->output_dir.c_str(),graph_id+1,iter);
    file = fopen(filename,"w+");
    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Equacao Transporte\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %d float\n",np);
    
    // Escrever pontos
    while (ptr != NULL)
    {
        fprintf(file,"%.10lf %.10lf %.10lf\n",ptr->x, ptr->y, ptr->z);
        ptr = ptr->next;
    }
    
    // Escrever linhas
    if (graph_id == 0)
    {
        ptr = g->getListNodes();
    }
    else if (graph_id == 1)
    {
        ptr = g2->getListNodes();
    }
    else if (graph_id == 2)
    {
        ptr = g3->getListNodes();
    }
    
    
    fprintf(file,"LINES %d %d\n",ne,ne*3);
    while (ptr != NULL)
    {
        Edge *ptrl = ptr->edges;
        while (ptrl != NULL)
        {
            fprintf(file,"2 %d %d\n",ptr->id,ptrl->dest->id);
            ptrl = ptrl->next;
        }
        ptr = ptr->next;
    }
    
    fprintf(file,"POINT_DATA %d\n",np);
    fprintf(file,"SCALARS concentracao float 1\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    
    // Escrever concentracoes
    if (graph_id == 0)
    {
        ptr = g->getListNodes();
//         int np = g->getTotalNodes();
    }
    else if (graph_id == 1)
    {
        ptr = g2->getListNodes();
//         int np = g2->getTotalNodes();
    }
    else if (graph_id == 2)
    {
        ptr = g3->getListNodes();
//         int np = g3->getTotalNodes();
    }
    
//     for (int j = 0; j < np; j++)
//     {
//         if (graph_id == 0)
//         {
//             fprintf(file,"%.10lf\n",cNew[j].concentration);
// //             printf("%f\n",cOld[j].concentration);
//         }
//         else if (graph_id == 1)
//         {
// //             printf("%f\n",cOld2[j].concentration);
//             fprintf(file,"%.10lf\n",cNew2[j].concentration);
//         }
//         else if (graph_id == 2)
//             fprintf(file,"%.10lf\n",cNew3[j].concentration);
//     }
    
    while (ptr != NULL)
    {
        if (graph_id == 0)
        {
            fprintf(file,"%.10lf\n",cNew[ptr->id].concentration);
//             printf("%f\n",cOld[ptr->id].concentration);
        }
        else if (graph_id == 1)
        {
//             printf("%f\n",cOld2[ptr->id].concentration);
            fprintf(file,"%.10lf\n",cNew2[ptr->id].concentration);
        }
        else if (graph_id == 2)
            fprintf(file,"%.10lf\n",cNew3[ptr->id].concentration);
        ptr = ptr->next;
    }
    
    // Escrever os raios
    if (graph_id == 0)
    {
        ptr = g->getListNodes();
//         int np = g->getTotalNodes();
    }
    else if (graph_id == 1)
    {
        ptr = g2->getListNodes();
//         int np = g2->getTotalNodes();
    }
    else if (graph_id == 2)
    {
        ptr = g3->getListNodes();
//         int np = g3->getTotalNodes();
    }
    
    fprintf(file,"CELL_DATA %d\n",(int)ne);
    fprintf(file,"scalars radius float\nLOOKUP_TABLE default\n");
    
    while (ptr != NULL)
    {
        Edge *ptrl = ptr->edges;
        while (ptrl != NULL)
        {
            fprintf(file,"%.10lf\n",ptrl->radius);
            ptrl = ptrl->next;
        }
        ptr = ptr->next;
    }
    
//     for (int i = 0; i < ; i++)
//         fprintf(outFile,"%.10lf\n",);
    
    fclose(file);
}

void Solver_Joao::writeVTKFile_static ()
{
    FILE *file;
    int np, ne;
    char filename[50];
    Node *ptr = g->getListNodes();
    np = g->getTotalNodes();
    ne = g->getTotalEdges();

    // Escrever o cabecalho do VTK
//     sprintf(filename,"VTK/sol%d.vtk",iter);
    sprintf(filename,"flow.vtk");
    file = fopen(filename,"w+");
    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Equacao Transporte\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %d float\n",np);
    while (ptr != NULL)
    {
        fprintf(file,"%.10lf %.10lf %.10lf\n",ptr->x, ptr->y, ptr->z);
        ptr = ptr->next;
    }
    
    ptr = g->getListNodes();
    fprintf(file,"LINES %d %d\n",ne,ne*3);
    while (ptr != NULL)
    {
        Edge *ptrl = ptr->edges;
        while (ptrl != NULL)
        {
            fprintf(file,"2 %d %d\n",ptr->id,ptrl->dest->id);
            ptrl = ptrl->next;
        }
        ptr = ptr->next;
    }
    
    fprintf(file,"CELL_DATA %d \n",ne);
    fprintf(file,"scalars flow float\nLOOKUP_TABLE default\n");
//     for (int i = 0; i < (int)np; i++)
//         fprintf(outFile,"%.10lf\n",mesh->elements[i].flow);
        
    ptr = g->getListNodes();
    while (ptr != NULL)
    {
        Edge *ptrl = ptr->edges;
        while (ptrl != NULL)
        {
            fprintf(file,"%f\n",ptrl->flow);
            ptrl = ptrl->next;
        }
        ptr = ptr->next;
    }
//     Edge * ptrl = noAtual->edges;
    fclose(file);
}


// void Solver_Joao::calc_velocity ()
// {
//     Node *ptr = g->getListNodes();
//     
//     while (ptr != NULL)
//     {
//         Edge *ptrl = ptr->edges;
//         while (ptrl != NULL)
//         {
//             double aux = calcNorm (ptrl->dest->x,ptrl->dest->y,ptrl->dest->z,ptr->x,ptr->y,ptr->z);
// //             fprintf(file,"2 %d %d\n",ptr->id,ptrl->dest->id);
//             ptrl->velocity = (0.133322387415*(ptrl->dest->pressure - ptr->pressure))/(aux);
//             ptrl = ptrl->next;
//         }
//         ptr = ptr->next;
//     }
// }

void Solver_Joao::setTypeTissue ()
{
    int i, j, k;
    
    for(i = 0; i < n_tissue; i++)
    {
        for(j = 0; j < n_tissue; j++)
        {
            cOld_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].type = 0;
            cNew_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].type = 0;
        }
    }
    
    int np;
    int p[2];
    
    // Verificar as conexoes do Grafo 1
    np = g->getTotalNodes();
    joao_grid_to_net_1.assign(np, std::vector<int>());
    
    for(k = 0; k < np; k++)
    {
        if(cOld[k].type == 1)
        {
            busca_proximal(k,p,0);
            cNew_tissue[GET_IDX_2D(p[0],p[1],n_tissue,n_tissue)].type = 1; // terminal
            cOld_tissue[GET_IDX_2D(p[0],p[1],n_tissue,n_tissue)].type = 1; // terminal
            joao_grid_to_net_1[k].push_back(p[0]);
            joao_grid_to_net_1[k].push_back(p[1]);
        }
    }
    
    // Verificar as conexoes do Grafo 2
    np = g2->getTotalNodes();
    joao_grid_to_net_2.assign(np, std::vector<int>());
    
    for(k = 0; k < np; k++)
    {
        if(cOld2[k].type == 1)
        {
            busca_proximal(k,p,1);
            cNew_tissue[GET_IDX_2D(p[0],p[1],n_tissue,n_tissue)].type = 1; // terminal
            cOld_tissue[GET_IDX_2D(p[0],p[1],n_tissue,n_tissue)].type = 1; // terminal
            joao_grid_to_net_2[k].push_back(p[0]);
            joao_grid_to_net_2[k].push_back(p[1]);
        }
    }
    
    // Verificar as conexoes do Grafo 3
    np = g3->getTotalNodes();
    joao_grid_to_net_3.assign(np, std::vector<int>());
    
    for(k = 0; k < np; k++)
    {
        if(cOld3[k].type == 1)
        {
            busca_proximal(k,p,2);
            cNew_tissue[GET_IDX_2D(p[0],p[1],n_tissue,n_tissue)].type = 1; // terminal
            cOld_tissue[GET_IDX_2D(p[0],p[1],n_tissue,n_tissue)].type = 1; // terminal
            joao_grid_to_net_3[k].push_back(p[0]);
            joao_grid_to_net_3[k].push_back(p[1]);
        }
    }
    
}

// Only executed for refined grids ...
void Solver_Joao::setTypeTissue_v2 ()
{
    int i, j, k;
    
    for(i = 0; i < n_tissue; i++)
    {
        for(j = 0; j < n_tissue; j++)
        {
            cOld_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].type = 0;
            cNew_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].type = 0;
        }
    }
    
    char filename[500];
    FILE *file;
    int np, nx, ny, value;
    int p[2];
    double h;
    
    // Verificar as conexoes do Grafo 1
    np = g->getTotalNodes();
    joao_grid_to_net_1.assign(np, std::vector<int>());
    
    //sprintf(filename, "coupling_map/graph_1_map_vent50x50.dat");
    sprintf(filename, "coupling_map/graph_1_map_vent25x25.dat");
    file = fopen(filename, "r");
    fscanf(file, "%d %d %lf", &nx, &ny, &h);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            fscanf(file, "%d", &value);
            if (value != 0) {
                joao_grid_to_net_1[value].push_back(i);
                joao_grid_to_net_1[value].push_back(j);
                cNew_tissue[GET_IDX_2D(i, j, n_tissue,n_tissue)].type = 1; // terminal
                cOld_tissue[GET_IDX_2D(i, j, n_tissue,n_tissue)].type = 1; // terminal
            }
        }
    }
    fclose(file);
    
    // Verificar as conexoes do Grafo 2
    np = g2->getTotalNodes();
    joao_grid_to_net_2.assign(np, std::vector<int>());
    
    //sprintf(filename, "coupling_map/graph_2_map_vent50x50.dat");
    sprintf(filename, "coupling_map/graph_2_map_vent25x25.dat");
    file = fopen(filename, "r");
    fscanf(file, "%d %d %lf", &nx, &ny, &h);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            fscanf(file, "%d", &value);
            if (value != 0) {
                joao_grid_to_net_2[value].push_back(i);
                joao_grid_to_net_2[value].push_back(j);
                cNew_tissue[GET_IDX_2D(i, j, n_tissue,n_tissue)].type = 1; // terminal
                cOld_tissue[GET_IDX_2D(i, j, n_tissue,n_tissue)].type = 1; // terminal
            }
        }
    }
    fclose(file);
    
    // Verificar as conexoes do Grafo 3
    np = g3->getTotalNodes();
    joao_grid_to_net_3.assign(np, std::vector<int>());
    
    //sprintf(filename, "coupling_map/graph_3_map_vent50x50.dat");
    sprintf(filename, "coupling_map/graph_3_map_vent25x25.dat");
    file = fopen(filename, "r");
    fscanf(file, "%d %d %lf", &nx, &ny, &h);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            fscanf(file, "%d", &value);
            if (value != 0) {
                joao_grid_to_net_3[value].push_back(i);
                joao_grid_to_net_3[value].push_back(j);
                cNew_tissue[GET_IDX_2D(i, j, n_tissue,n_tissue)].type = 1; // terminal
                cOld_tissue[GET_IDX_2D(i, j, n_tissue,n_tissue)].type = 1; // terminal
            }
        }
    }
    fclose(file);   
}

void Solver_Joao::busca_proximal(int id, int p[], int graph_id)
{
    int i, j;
    double aux, dist = 100000.;
    
    Node * ptr;
    
    if(graph_id == 0)
        ptr = g->searchNode(id);
    else if(graph_id == 1)
        ptr = g2->searchNode(id);
    else if(graph_id == 2)
        ptr = g3->searchNode(id);
    
    for(i = (ptr->y - h_gadolinio)/h_gadolinio; i <= (ptr->y + h_gadolinio)/h_gadolinio; i++)
    {
        for(j = (ptr->x - h_gadolinio)/h_gadolinio; j <= (ptr->x + h_gadolinio)/h_gadolinio; j++)
        {
            aux = calcNorm(ptr->x,ptr->y,ptr->z,j*h_gadolinio,i*h_gadolinio,0.);
            if(aux < dist)
            {
                dist = aux;
                p[0] = i;
                p[1] = j;
            }
        }
    }
}

void Solver_Joao::search_region(int id, int p[], int graph_id, double * peso)
{
    int i, j;
    double aux, dist = 100000.;
    
    Node * ptr;
    
    if(graph_id == 0)
        ptr = g->searchNode(id);
    else if(graph_id == 1)
        ptr = g2->searchNode(id);
    else if(graph_id == 2)
        ptr = g3->searchNode(id);
    
    for(i = (ptr->y - h_gadolinio)/h_gadolinio; i <= (ptr->y + h_gadolinio)/h_gadolinio; i++)
    {
        for(j = (ptr->x - h_gadolinio)/h_gadolinio; j <= (ptr->x + h_gadolinio)/h_gadolinio; j++)
        {
            aux = calcNorm(ptr->x,ptr->y,ptr->z,j*h_gadolinio,i*h_gadolinio,0.);
            if(aux < dist)
            {
                dist = aux;
                p[0] = i;
                p[1] = j;
            }
        }
    }
    
    if(graph_id == 0)
    {
        if((peso[GET_IDX_2D(p[0],p[1],n_tissue,n_tissue)] >= 0.75 && p[1] < n_tissue/3. && p[0] < 10.*n_tissue/16. && p[0] > 6.*n_tissue/16.) && (cOld[id].type == 1)/*terminal*/)
            cOld[id].roi = 0;
        else
            cOld[id].roi = 1;
    }
    
    else if(graph_id == 1)
    {
        if((peso[GET_IDX_2D(p[0],p[1],n_tissue,n_tissue)] >= 0.75 && p[1] < n_tissue/3. && p[0] < 10.*n_tissue/16. && p[0] > 6.*n_tissue/16.) && (cOld2[id].type == 1)/*terminal*/)
            cOld2[id].roi = 0;
        else
            cOld2[id].roi = 1;
    }
    
    else if(graph_id == 2)
    {
        if((peso[GET_IDX_2D(p[0],p[1],n_tissue,n_tissue)] >= 0.75 && p[1] < n_tissue/3. && p[0] < 10.*n_tissue/16. && p[0] > 6.*n_tissue/16.) && (cOld3[id].type == 1)/*terminal*/)
            cOld3[id].roi = 0;
        else
            cOld3[id].roi = 1;
    }
//     else if(tecido[GET_IDX_2D(p[0],p[1],n_tissue,n_tissue)] == 50)
//         cOld[id].roi = 1;
//     else
//         cOld[id].roi = 2;
}

int Solver_Joao::number_tissue()
{
    int i, j;
    int N = 0;
    
    for(i = 0; i < n_tissue; i++)
    {
        for(j = 0; j < n_tissue; j++)
        {
            if(tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 128 || tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 50)
                N++;
        }
    }
    return N;
}

//double Solver_Joao::integral_extra (int N)
double Solver_Joao::integral_extra ()
{
    int i, j;

    double V = 0.;
    int n_active = 0;

    for (i = 1; i < n_tissue-1; i++)
    {
        for (j = 1; j < n_tissue-1; j++)
        {
            if(tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 128 || \
            tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 50) // tecido sadio ou tecido infartado
            {
                V += cNew_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration;
                n_active++;
            }
        }
    }

    // V = V/10.; // soma das concentrações dividida pelo número de terminais
    //V = V/N;
    V = V / (double)n_active;
    return V;
}

double Solver_Joao::integral_intra ()
{
    int j;
    double V = 0.;
    
    int np;
    int n_terminais = 0;;
    
    np = g->getTotalNodes();    
        
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

    return V/n_terminais;
}

void Solver_Joao::solve1D_source(double * C_1D, double V, double adv_1D, double difusao_1D)
{
    int i;
    double * new_C_1D = new double[n_1D];

    C_1D[0] = V; //+= ou +?
//    fou1D(C_1D, P_1D, vswap, C_face_1D, v_1D);
    double dif_1D_global = 0.0;

    for (i = 1; i < n_1D; i++)
    {
        if(i != n_1D-1)
        {
            adv_1D = (v_1D * C_1D[i] - v_1D * C_1D[i-1])/h_1D;

            difusao_1D = (dif_1D_global * (C_1D[i+1] - C_1D[i] )
                        - dif_1D_global * (C_1D[i] - C_1D[i-1] )) / (h_1D*h_1D);
        }

        else
        {
            adv_1D = (v_1D * C_1D[i] - v_1D * C_1D[i-1])/h_1D;
            difusao_1D = 0.;
        }

        if(i == n_1D-1)
            new_C_1D[i] = (- adv_1D - rins * C_1D[i] - (dif_1D_global * (C_1D[i] - C_1D[i-1])) / (h_1D*h_1D)) * dt / (phi_1D) + C_1D[i];

        else
            new_C_1D[i] = (difusao_1D - adv_1D - rins * C_1D[i]) * dt / (phi_1D) + C_1D[i];
    }

    for (i = 1; i < n_1D; i++)
    {
        C_1D[i] = new_C_1D[i];
    }
}

double Solver_Joao::calcula_fonte(double * C_1D, double t)
{
    double fonte_new = (1./(sigma*sqrt(2*M_PI))) * exp(-0.5*pow((t - T_peak)/sigma,2));
    double sumidouro_new;
    
    if (C_1D[n_1D-1] - fonte_new < 0.)
        fonte_new += 0.;
    else
        fonte_new += - r3 * (C_1D[n_1D-1] - fonte_new);

    return fonte_new;
}

double Solver_Joao::calcula_fonte_arteria(double * C_1D, int i, double dt, double C)
{
    double fonte_new = (1./(sigma*sqrt(2*M_PI))) * exp(-0.5*pow((i*dt - T_peak)/sigma,2));
    
    if (C_1D[n_1D-1] - fonte_new < 0.)
        C = 0.;
    else
        C = r3 * (C_1D[n_1D-1] - fonte_new);
    
    fonte_new += C;

    return fonte_new;
}

void Solver_Joao::solve ()
{
    char filename_Ce[500];
    sprintf(filename_Ce, "%s/Ce.txt", this->output_dir.c_str());
    FILE *file_Ce = fopen(filename_Ce, "w+");
    Node *ptr= g->getListNodes();

    writeVTKFile_static ();

    double dif_tissue = 0.;
    int i, j, p, q;

    int N = number_tissue();
//     printf("%d\n",N);
    #ifdef OUTPUT
    printf("[!] Solving transient problem ... \n");
    printf("[!] Progress\n");
    fflush(stdout);
    #endif

//     calc_velocity();
//     print();

    int np;
    double SI_health_tissue = 0.;
    double SI_ischemic_tissue = 0.;
    double SI_health_intra = 0.;
    
    // Inicializa os termos de comunicacao do tecido
    double * reaction = (double*)malloc(sizeof(double)*n_tissue*n_tissue);
    double * peso =  new double[n_tissue*n_tissue];

    int u[2];
     
    for(i = 0; i < n_tissue; i++)
    {
        for(j = 0; j < n_tissue; j++)
        {
            if(tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 255 || tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 0)
            {
                cNew_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration = -0.;
                cOld_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration = -0.;
            }

            else 
            {
                cNew_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration = 0.;
                cOld_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration = 0.;
            }
            reaction[GET_IDX_2D(i,j,n_tissue,n_tissue)] = 0.;
            peso[GET_IDX_2D(i,j,n_tissue,n_tissue)] = 0.;
        }
    }
    
    // Parte 1D da recirculação
//     int m = 0;
    double difusao_1D = 0., adv_1D = 0., V = 0.;
    double * C_1D =  new double[n_1D*2];
    for(i = 0; i < 2*n_1D; i++)
    {
        C_1D[i] = 0.;
    }
//     double C = 0.;

    
    char tissue[100];
    
    FILE *gnuplot;
    char filename[50];
    
    contorno_perm_het(peso, tecido);
    peso = perm_het(peso, tecido);
    
    for(i = 0; i < 3; i++)
    {
        if(i == 0)
        {
//             ptr = g->getListNodes();
            np = g->getTotalNodes();
        }
        else if(i == 1)
        {
//             ptr = g2->getListNodes();
            np = g2->getTotalNodes();
        }
        else if(i == 2)
        {
//             ptr = g3->getListNodes();
            np = g3->getTotalNodes();
        }
        
        for(j = 0; j < np; j++)
        {
            search_region(j, u, i, peso);
        }
    }
    
    double aux = 0.;
    double * SI_health_vec = new double[M/SAVEEACH];
    double * SI_ischemic_vec = new double[M/SAVEEACH];
//         
    // Time loop
    for (int i = 0; i < M; i++)
    {
//         SI_health_tissue = 0.;
//         SI_ischemic_tissue = 0.;
        SI_health_intra = 0.;
        double t = i*dt;
        //recirculação
        // Calcula a média da concentração no domínio
        V = integral_intra();
//         V = C/N;
        
        // Resolve a parte 1D (Recirculação), com V sendo a condição de contorno prescrita no lado esquerdo
        solve1D_source(C_1D, V, adv_1D, difusao_1D);
        
        // Calcula o termo de fonte que será atribudo como fluxo no contorno
        double fonte_new = calcula_fonte(C_1D, t);
//         double fonte_new = calcula_fonte_arteria(C_1D, i, dt, C);
        
        // Print the progress of the solution
        #ifdef OUTPUT
        printProgress(i,M);
        #endif

        // Write the solution data to a file
        //writePlotData(t);

        // Write the solution to .vtk file
        #ifdef VTK
        for (int k = 0; k < 3; k++)
        {
            if (i % SAVEEACH == 0)
            {
                writeVTKFile(i,k);
            }
        }
        
        if (i % SAVEEACH == 0)
        {
            sprintf(tissue, "%s/tissue/tissue%d.vtk", this->output_dir.c_str(), i);
            writeVTKFile_tissue (tissue, n_tissue, n_tissue);
        }
        #endif                
//         C = 0.;
        
        // Resolver as concentracoes
        
        // Grafo 1
        ptr = g->getListNodes();
        np = g->getTotalNodes();
        
        for (int j = 0; j < np; j++)
        {
//             if(j == 217)
//                 printf("%d\n",cOld[j].type);
            if (cOld[j].type == 0) // raíz
            {
                cNew[j].concentration = fonte_new;
//                 printf("%.5f - ", cNew[j].concentration);
            }
//                 cNew[j].concentration = gaussiana(t);
            else if (cOld[j].type == 1) // terminal
            {
                // Mapa com as coordenadas (i,j) dos vizinhos no tecido
                // do terminal 'j'
                std::vector<std::pair<int,int>> u_arr;  

                // BEFORE REVISION - Malha grossa
                //busca_proximal(j,u,0);
                search_proximal_in_map(j, u, 0);    // Busca O(1)
                calcula_reaction(reaction, j, u[0], u[1],0);
                cNew[j].concentration = eq_terminal(ptr, j, u[0], u[1], reaction, 0);

                // REVISION CODE - Somente para malha refinada
                //search_proximal_in_map_v2(j, u_arr, 0);
                //calcula_reaction_v2(reaction, j, u_arr, 0);
                //cNew[j].concentration = eq_terminal_v2(ptr, j, u_arr, reaction, 0);
                                
            }
            else if (cOld[j].type == 2) // normal
                cNew[j].concentration = eq_geral(ptr,j,0);
            
            else if (cOld[j].type == 3) // bifurcação
                cNew[j].concentration = eq_bifurcacao(ptr,j,0);

            // Passar para a proxima celula    
            // CAUSA DO BUG ! FALTOU COLOCAR ISSO !
            
            if(cOld[j].roi == 0)
                    SI_health_intra += cNew[j].concentration;  
            
            
            ptr = ptr->next;
        }
        
        // Grafo 2
        ptr = g2->getListNodes();
        np = g2->getTotalNodes();
        
        for (int j = 0; j < np; j++)
        {
//             if(j == 217)
//                 printf("%d\n",cOld[j].type);
            if (cOld2[j].type == 0) // raíz
            {
                cNew2[j].concentration = fonte_new;
//                 printf("%.5f - ", cNew2[j].concentration);
            }
//                 cNew[j].concentration = gaussiana(t);
            else if (cOld2[j].type == 1) // terminal
            {
                std::vector<std::pair<int,int>> u_arr;

                // BEFORE REVISION
                //busca_proximal(j,u,1);
                search_proximal_in_map(j, u, 1);
                calcula_reaction(reaction, j, u[0], u[1],1);
                cNew2[j].concentration = eq_terminal(ptr, j, u[0], u[1], reaction, 1);

                //search_proximal_in_map_v2(j, u_arr, 1);
                //calcula_reaction_v2(reaction, j, u_arr, 1);
                //cNew2[j].concentration = eq_terminal_v2(ptr, j, u_arr, reaction, 1);
            }
            else if (cOld2[j].type == 2) // normal
                cNew2[j].concentration = eq_geral(ptr,j,1);
            
            else if (cOld2[j].type == 3) // bifurcação
                cNew2[j].concentration = eq_bifurcacao(ptr,j,1);


            if(cOld2[j].roi == 0)
                    SI_health_intra += cNew2[j].concentration;
            
            // Passar para a proxima celula    
            // CAUSA DO BUG ! FALTOU COLOCAR ISSO !
            ptr = ptr->next;
        }
        
        // Grafo 3
        ptr = g3->getListNodes();
        np = g3->getTotalNodes();
        
        for (int j = 0; j < np; j++)
        {
//             if(j == 217)
//                 printf("%d\n",cOld[j].type);
            if (cOld3[j].type == 0) // raíz
            {
                cNew3[j].concentration = fonte_new;
//                 printf("%.5f\n", cNew3[j].concentration);
            }
//                 cNew[j].concentration = gaussiana(t);
            else if (cOld3[j].type == 1) // terminal
            {
                std::vector<std::pair<int,int>> u_arr;

                // BEFORE REVISION
                //busca_proximal(j,u,2);
                search_proximal_in_map(j, u, 2);
                calcula_reaction(reaction, j, u[0], u[1],2);
                cNew3[j].concentration = eq_terminal(ptr, j, u[0], u[1], reaction, 2);

                //search_proximal_in_map_v2(j, u_arr, 2);
                //calcula_reaction_v2(reaction, j, u_arr, 2);
                //cNew3[j].concentration = eq_terminal_v2(ptr, j, u_arr, reaction, 2);
            }
            else if (cOld3[j].type == 2) // normal
                cNew3[j].concentration = eq_geral(ptr,j,2);
            
            else if (cOld3[j].type == 3) // bifurcação
                cNew3[j].concentration = eq_bifurcacao(ptr,j,2);
            
            if(cOld3[j].roi == 0)
                SI_health_intra += cNew3[j].concentration;

            
            // Passar para a proxima celula    
            // CAUSA DO BUG ! FALTOU COLOCAR ISSO !
            ptr = ptr->next;
        }
/*
// DEBUG
        char filename[500];
        //sprintf(filename, "outputs/sol_Joao_graph3_iter_%d.txt", joao_counter);
        sprintf(filename, "outputs/reaction_Joao_tiss_iter_%d.txt", joao_counter);
        FILE *file = fopen(filename, "w+");
        //for (int j = 0; j < np; j++) {
        //    fprintf(file, "%g\n", c_new[j]);
        //}
        for (int j = 0; j < n_tissue; j++) {
            for (int i = 0; i < n_tissue; i++) {
                fprintf(file, "%g\n", reaction[GET_IDX_2D(j, i, n_tissue, n_tissue)]);
            }
        }
        fclose(file);

        // DEBUG
        if (joao_counter == 10) {
            exit(1);
        }
        else {
            joao_counter++;
        }
*/

        for(p = 0; p < n_tissue; p++)
        {
            for(q = 0; q < n_tissue; q++)
            {
                if(tecido[GET_IDX_2D(p,q,n_tissue,n_tissue)] == 128 || tecido[GET_IDX_2D(p,q,n_tissue,n_tissue)] == 50) 
                {
                    if((tecido[GET_IDX_2D(p+1,q,n_tissue,n_tissue)] == 255 && tecido[GET_IDX_2D(p,q+1,n_tissue,n_tissue)] == 255) ||
                       (tecido[GET_IDX_2D(p+1,q,n_tissue,n_tissue)] == 0 && tecido[GET_IDX_2D(p,q+1,n_tissue,n_tissue)] == 0))
                        dif_tissue = (D/(h_gadolinio*h_gadolinio)) * 
                        (- (cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q-1,n_tissue,n_tissue)].concentration) 
                         - (cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p-1,q,n_tissue,n_tissue)].concentration));
                    
                    else if((tecido[GET_IDX_2D(p+1,q,n_tissue,n_tissue)] == 255 && tecido[GET_IDX_2D(p,q-1,n_tissue,n_tissue)] == 255) ||
                            (tecido[GET_IDX_2D(p+1,q,n_tissue,n_tissue)] == 0 && tecido[GET_IDX_2D(p,q-1,n_tissue,n_tissue)] == 0))
                        dif_tissue = (D/(h_gadolinio*h_gadolinio)) * 
                        (+ (cOld_tissue[GET_IDX_2D(p,q+1,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration) 
                         - (cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p-1,q,n_tissue,n_tissue)].concentration));
                        
                    else if((tecido[GET_IDX_2D(p-1,q,n_tissue,n_tissue)] == 255 && tecido[GET_IDX_2D(p,q+1,n_tissue,n_tissue)] == 255) ||
                            (tecido[GET_IDX_2D(p-1,q,n_tissue,n_tissue)] == 0 && tecido[GET_IDX_2D(p,q+1,n_tissue,n_tissue)] == 0))
                        dif_tissue = (D/(h_gadolinio*h_gadolinio)) * 
                        (- (cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q-1,n_tissue,n_tissue)].concentration) 
                         + (cOld_tissue[GET_IDX_2D(p+1,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration));
                    
                    else if((tecido[GET_IDX_2D(p-1,q,n_tissue,n_tissue)] == 255 && tecido[GET_IDX_2D(p,q-1,n_tissue,n_tissue)] == 255) ||
                            (tecido[GET_IDX_2D(p-1,q,n_tissue,n_tissue)] == 0 && tecido[GET_IDX_2D(p,q-1,n_tissue,n_tissue)] == 0))
                        dif_tissue = (D/(h_gadolinio*h_gadolinio)) * 
                        (+ (cOld_tissue[GET_IDX_2D(p,q+1,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration) 
                         + (cOld_tissue[GET_IDX_2D(p+1,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration));
                    
                    else if((tecido[GET_IDX_2D(p,q+1,n_tissue,n_tissue)] == 255) || (tecido[GET_IDX_2D(p,q+1,n_tissue,n_tissue)] == 0))
                        dif_tissue = (D/(h_gadolinio*h_gadolinio)) * 
                        (- (cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q-1,n_tissue,n_tissue)].concentration) 
                         + (cOld_tissue[GET_IDX_2D(p+1,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration)
                         - (cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p-1,q,n_tissue,n_tissue)].concentration));
                    
                    else if((tecido[GET_IDX_2D(p,q-1,n_tissue,n_tissue)] == 255) || (tecido[GET_IDX_2D(p,q-1,n_tissue,n_tissue)] == 0))
                        dif_tissue = (D/(h_gadolinio*h_gadolinio)) * 
                        (+ (cOld_tissue[GET_IDX_2D(p,q+1,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration) 
                         + (cOld_tissue[GET_IDX_2D(p+1,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration)
                         - (cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p-1,q,n_tissue,n_tissue)].concentration));
                        
                    else if((tecido[GET_IDX_2D(p+1,q,n_tissue,n_tissue)] == 255) || (tecido[GET_IDX_2D(p+1,q,n_tissue,n_tissue)] == 0))
                        dif_tissue = (D/(h_gadolinio*h_gadolinio)) * 
                        (+ (cOld_tissue[GET_IDX_2D(p,q+1,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration) 
                         - (cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q-1,n_tissue,n_tissue)].concentration)
                         - (cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p-1,q,n_tissue,n_tissue)].concentration));
                    
                    else if((tecido[GET_IDX_2D(p-1,q,n_tissue,n_tissue)] == 255) || (tecido[GET_IDX_2D(p-1,q,n_tissue,n_tissue)] == 0))
                        dif_tissue = (D/(h_gadolinio*h_gadolinio)) * 
                        (+ (cOld_tissue[GET_IDX_2D(p,q+1,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration) 
                         - (cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q-1,n_tissue,n_tissue)].concentration)
                         + (cOld_tissue[GET_IDX_2D(p+1,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration));
                    
                    else
                        dif_tissue = (D/(h_gadolinio*h_gadolinio)) * 
                       (+ (cOld_tissue[GET_IDX_2D(p,q+1,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration)
                        - (cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q-1,n_tissue,n_tissue)].concentration)
                        + (cOld_tissue[GET_IDX_2D(p+1,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration)
                        - (cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration - cOld_tissue[GET_IDX_2D(p-1,q,n_tissue,n_tissue)].concentration));
                    
                    // Este IF exclui todos os pontos do contorno.
                    // Dentro dele, prescreve-se fluxo pras faces dos nós na "segunda camada" do contorno
//                     if(!(((tecido[GET_IDX_2D(p+1,q,n_tissue,n_tissue)] == 0 && tecido[GET_IDX_2D(p,q+1,n_tissue,n_tissue)] == 0) || (tecido[GET_IDX_2D(p+1,q,n_tissue,n_tissue)] == 255 && tecido[GET_IDX_2D(p,q+1,n_tissue,n_tissue)] == 255)) &&
//                          ((tecido[GET_IDX_2D(p+1,q,n_tissue,n_tissue)] == 0 && tecido[GET_IDX_2D(p,q-1,n_tissue,n_tissue)] == 0) || (tecido[GET_IDX_2D(p+1,q,n_tissue,n_tissue)] == 255 &&    tecido[GET_IDX_2D(p,q-1,n_tissue,n_tissue)] == 255)) &&
//                          ((tecido[GET_IDX_2D(p-1,q,n_tissue,n_tissue)] == 0 && tecido[GET_IDX_2D(p,q+1,n_tissue,n_tissue)] == 0) || (tecido[GET_IDX_2D(p-1,q,n_tissue,n_tissue)] == 255 && tecido[GET_IDX_2D(p,q+1,n_tissue,n_tissue)] == 255)) &&
//                          ((tecido[GET_IDX_2D(p-1,q,n_tissue,n_tissue)] == 0 && tecido[GET_IDX_2D(p,q-1,n_tissue,n_tissue)] == 0) || (tecido[GET_IDX_2D(p-1,q,n_tissue,n_tissue)] == 255 && tecido[GET_IDX_2D(p,q-1,n_tissue,n_tissue)] == 255)) &&
//                           (tecido[GET_IDX_2D(p,q+1,n_tissue,n_tissue)] == 255) && (tecido[GET_IDX_2D(p,q+1,n_tissue,n_tissue)] == 0) && (tecido[GET_IDX_2D(p,q-1,n_tissue,n_tissue)] == 255) && (tecido[GET_IDX_2D(p,q-1,n_tissue,n_tissue)] == 0) &&
//                           (tecido[GET_IDX_2D(p+1,q,n_tissue,n_tissue)] == 255) && (tecido[GET_IDX_2D(p+1,q,n_tissue,n_tissue)] == 0) && (tecido[GET_IDX_2D(p-1,q,n_tissue,n_tissue)] == 255) && (tecido[GET_IDX_2D(p-1,q,n_tissue,n_tissue)] == 0) ))
//                     {
                        
                        if(cNew_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].type != 1) // não terminal
                            reaction[GET_IDX_2D(p,q,n_tissue,n_tissue)] = 0.;
                                        
                        cNew_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration = (dif_tissue - reaction[GET_IDX_2D(p,q,n_tissue,n_tissue)]/((1.0-phi)*lambda)
                                        - decaimento*cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration/lambda) * dt /*/ ((1.0 - phi)*lambda)*/
                                        + cOld_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration;
                                        
//                         cNew_tissue[GET_IDX_2D(p,q,n_tissue,n_tissue)].concentration *= (1.0 - phi)*lambda;
                                                
//                     }
                }                
            }
        }
        
        // Cálculo da Intensidade do Sinal do Tecido (extravascular)
        

        if(i % SAVEEACH == 0)
        {
            SI_health_tissue = signal_intensity_healthy(i, tecido, cNew_tissue, peso, SI_health_tissue);//, SI_ischemic_tissue);
            SI_ischemic_tissue = signal_intensity_ischemic(i, tecido, cNew_tissue, peso, SI_ischemic_tissue);
            
//             aux = SI_health_tissue;
            
            SI_ischemic_vec[i/SAVEEACH] = SI_ischemic_tissue;
            SI_health_vec[i/SAVEEACH] = SI_health_tissue;
        }
        
//         if(SI_health_tissue > aux)
//            aux = SI_health_tissue;
        
//         if(i % SAVEEACH == 0)
//         {
//             sprintf(filename,"sol/gauss%.3f.txt",decaimento);
//             gnuplot = fopen(filename,"a");
// //             fprintf(gnuplot,"%f %f %f %f\n", i*dt, fator_SI*(phi*SI_health_intra + SI_health_tissue), fator_SI*phi*SI_health_intra, fator_SI*SI_health_tissue);
//             fprintf(gnuplot,"%f %f %f \n", i*dt, fator_SI*SI_health_tissue, fator_SI*SI_ischemic_tissue);
//         }
        
        // Computar o valor da concentracao extravascular
        double Ce = integral_extra();
        fprintf(file_Ce, "%g %g\n", t, Ce);

        // Pular para a proxima iteracao
        nextTimestep();
            
    }

        
    double aux_1 = 0.;
    double aux_2 = 0.;
    
    for(i = 0; i < M; i++)
    {
        if(i % SAVEEACH == 0 || i == 0)
        {
            aux_1 = SI_health_vec[i/SAVEEACH];
            
            if(aux_1 > aux_2)
                aux_2 = aux_1;
        }
    }
    
    double fator = 83.63226015538544/aux_2;
//     double fator = 28.386583;
    printf("\n \n %f %f \n \n", aux_2, fator); 
    
    for(i = 0; i < M; i++)
    {
        if(i % SAVEEACH == 0 || i == 0)
        {
            sprintf(filename,"%s/gauss.txt", this->output_dir.c_str());
            gnuplot = fopen(filename,"a");
//             fprintf(gnuplot,"%f %f %f %f\n", i*dt, fator_SI*(phi*SI_health_intra + SI_health_tissue), fator_SI*phi*SI_health_intra, fator_SI*SI_health_tissue);
//             fprintf(gnuplot,"%f %f %f \n", i*dt, fator_SI*SI_health_tissue, fator_SI*SI_ischemic_tissue);
            
            fprintf(gnuplot, "%.4f %.4f %.4f \n", i*dt, fator*SI_health_vec[i/SAVEEACH], fator*SI_ischemic_vec[i/SAVEEACH]);
        }
    }
    
    free(SI_health_vec);
    free(SI_ischemic_vec);
    fclose(gnuplot);
    fclose(file_Ce);
}

void Solver_Joao::calcula_reaction(double * reaction, int k, int i, int j, int graph_id)
{
    if(graph_id == 0)
    {
        if(cOld[k].concentration < cOld_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration)
            reaction[GET_IDX_2D(i,j,n_tissue,n_tissue)] = 0.;

        else
        {
            if(!(tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 50))
                reaction[GET_IDX_2D(i,j,n_tissue,n_tissue)] = - r * (cOld[k].concentration - cOld_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration);

            else
                reaction[GET_IDX_2D(i,j,n_tissue,n_tissue)] = - r_problema * (cOld[k].concentration - cOld_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration);
        }
    }
    
    else if(graph_id == 1)
    {
        if(cOld2[k].concentration < cOld_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration)
            reaction[GET_IDX_2D(i,j,n_tissue,n_tissue)] = 0.;

        else
        {
            if(!(tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 50))
                reaction[GET_IDX_2D(i,j,n_tissue,n_tissue)] = - r * (cOld2[k].concentration - cOld_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration);

            else
                reaction[GET_IDX_2D(i,j,n_tissue,n_tissue)] = - r_problema * (cOld2[k].concentration - cOld_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration);
        }
    }
    
    else if(graph_id == 2)
    {
        if(cOld3[k].concentration < cOld_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration)
            reaction[GET_IDX_2D(i,j,n_tissue,n_tissue)] = 0.;

        else
        {
            if(!(tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 50))
                reaction[GET_IDX_2D(i,j,n_tissue,n_tissue)] = - r * (cOld3[k].concentration - cOld_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration);

            else
                reaction[GET_IDX_2D(i,j,n_tissue,n_tissue)] = - r_problema * (cOld3[k].concentration - cOld_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration);
        }
    }
}

// void Solver_Joao::leGeometria()
// {
//   int p = 0;
//   std::fstream fileIn;
//   char nome[100];
// 
//   sprintf(nome, "vent%dx%dInfarto.dat",n_tissue, n_tissue);
// 
//   fileIn.open(nome, std::ios_base::in);
//     for(int j = 0; j < n_tissue; j++)
//     {
//         for(int i = 0; i < n_tissue; i++)
//         {
//             fileIn>>tecido[GET_IDX_2D(j,i,n_tissue, n_tissue)];
// //             printf("%d ",tecido[GET_IDX_2D(j,i,n_tissue, n_tissue)]);
//             if (tecido[GET_IDX_2D(j,i,n_tissue,n_tissue)] == 128 || tecido[GET_IDX_2D(j,i,n_tissue,n_tissue)] == 50)
//                 p++;
//         }
//     }
//     printf("\nnúmero de pontos = %d\n",p);
//     fileIn.close();
// }

double * Solver_Joao::contorno_perm_het(double * p, int * tecido)
{
    int i, j;

    for(i = 1; i < n_tissue-1; i++)
    {
        for(j = 1; j < n_tissue-1; j++)
        {
            if(tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 128 || tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 50)
            {
                if(tecido[GET_IDX_2D(i+1,j,n_tissue,n_tissue)] == 255 || tecido[GET_IDX_2D(i-1,j,n_tissue,n_tissue)] == 255 || tecido[GET_IDX_2D(i,j+1,n_tissue,n_tissue)] == 255 || tecido[GET_IDX_2D(i,j-1,n_tissue,n_tissue)] == 255)
                    p[GET_IDX_2D(i,j,n_tissue,n_tissue)] = 0.;

                else if (tecido[GET_IDX_2D(i+1,j,n_tissue,n_tissue)] == 0 || tecido[GET_IDX_2D(i-1,j,n_tissue,n_tissue)] == 0 || tecido[GET_IDX_2D(i,j+1,n_tissue,n_tissue)] == 0 || tecido[GET_IDX_2D(i,j-1,n_tissue,n_tissue)] == 0)
                    p[GET_IDX_2D(i,j,n_tissue,n_tissue)] = 1.;
            }
        }
    }
    return p;
}

double * Solver_Joao::perm_het (double * p, int * tecido)
{
    int i, j, m = 0;
    double max = 0., aux1, aux2, erro = 1.;

    /* ----- Calculo das pressões -----*/
    while(erro > 0.000001)
    {
        max = 0.;
        for(i = 1; i < n_tissue-1; i++)
        {
            for(j = 1; j < n_tissue-1; j++)
            {
              if((tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 128 || tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 50) && !(tecido[GET_IDX_2D(i+1,j,n_tissue,n_tissue)] == 255 || tecido[GET_IDX_2D(i-1,j,n_tissue,n_tissue)] == 255 || tecido[GET_IDX_2D(i,j-1,n_tissue,n_tissue)] == 255 || tecido[GET_IDX_2D(i,j+1,n_tissue,n_tissue)] == 255)
                && !(tecido[GET_IDX_2D(i+1,j,n_tissue,n_tissue)] == 0 || tecido[GET_IDX_2D(i-1,j,n_tissue,n_tissue)] == 0 || tecido[GET_IDX_2D(i,j-1,n_tissue,n_tissue)] == 0 || tecido[GET_IDX_2D(i,j+1,n_tissue,n_tissue)] == 0))
              {
                  aux1 = p[GET_IDX_2D(i,j,n_tissue,n_tissue)];

                  p[GET_IDX_2D(i,j,n_tissue,n_tissue)] = (p[GET_IDX_2D(i,j+1,n_tissue,n_tissue)] + p[GET_IDX_2D(i,j-1,n_tissue,n_tissue)] + p[GET_IDX_2D(i+1,j,n_tissue,n_tissue)] + p[GET_IDX_2D(i-1,j,n_tissue,n_tissue)])/4.;

                  aux2 = fabs(p[GET_IDX_2D(i,j,n_tissue,n_tissue)] - aux1);

                  if(aux2 > max)
                      max = aux2;
              }
              else if(tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 255 || tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 0)
              {
                  p[GET_IDX_2D(i,j,n_tissue,n_tissue)] = -1.;
              }
            }
        }

        if(max < erro)
            erro = max;
        m++;
//        printf("erro =  %f \n", erro);
    }
//    printf("%d \n", m);
    return p;
}


double Solver_Joao::signal_intensity_healthy(int k, int * tecido, Cell_Joao * S, double * peso, double SI_health)
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

double Solver_Joao::signal_intensity_ischemic(int k, int * tecido, Cell_Joao * S, double * peso, double SI_ischemic)
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

    for(i = 0; i <= n_tissue-1; i++)
    {
        for(j = 0; j <= n_tissue-1; j++)
        {
            if(tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 128 || tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 50)
            {
                if(tecido[GET_IDX_2D(i,j,n_tissue,n_tissue)] == 50  && peso[GET_IDX_2D(i,j,n_tissue,n_tissue)] >= 0.95)
                {
                    SI_ischemic += /*(1.0 - phi) * lambda_infarto * */S[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration;
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

void Solver_Joao::leGeometria (string filename)
{
    int p = 0;
    int nx, ny;
    double h;
    ifstream in(filename.c_str());
//     if (!in) error("Cannot open MSH file");

    in >> nx >> ny >> h;
    n_tissue = nx;
    h_gadolinio = h;

    tecido = (int*)malloc(sizeof(int)*n_tissue*n_tissue);
    for(int j = 0; j < n_tissue; j++)
    {
        for(int i = 0; i < n_tissue; i++)
        {
            in >> tecido[GET_IDX_2D(j,i,n_tissue,n_tissue)];
//                 printf("%d ",tecido[GET_IDX_2D(j,i,n_tissue, n_tissue)]);
            if (tecido[GET_IDX_2D(j,i,n_tissue,n_tissue)] == 128 || tecido[GET_IDX_2D(j,i,n_tissue,n_tissue)] == 50)
                p++;
        }
    }
//     printf("\nnúmero de pontos = %d\n",p);

    in.close();
}

double Solver_Joao::gaussiana (double t)
{
    return (1./(sigma*sqrt(2*M_PI))) * exp(-0.5*pow((t - T_peak)/sigma,2));
}

double Solver_Joao::eq_terminal (Node *ptr, int k, int i, int j, double * reaction, int graph_id)
{
    int id_viz = ptr->edges->id;
    double flow = calcula_flow(ptr);
    double radius = calcula_radius(ptr);

    #ifdef PRESSURE
//     double term1 = - ptr->velocity * dt / dx;
    double term1 = - ((flow/FACTOR_FLOW) * dt / dx)/(M_PI*pow(radius,2));
    #else 
    double term1 = - 1. * dt / dx;
    #endif
    
    if(graph_id == 0)
    {
        double term2 = cOld[k].concentration - cOld[id_viz].concentration;
        double term3 = reaction[GET_IDX_2D(i,j,n_tissue,n_tissue)] * dt / (dx/**phi*/);

        return (term1 * term2 + term3 + cOld[k].concentration);
    }
    else if(graph_id == 1)
    {
        double term2 = cOld2[k].concentration - cOld2[id_viz].concentration;
        double term3 = reaction[GET_IDX_2D(i,j,n_tissue,n_tissue)] * dt / (dx/**phi*/);

        return (term1 * term2 + term3 + cOld2[k].concentration);
    }
    else if(graph_id == 2)
    {
        double term2 = cOld3[k].concentration - cOld3[id_viz].concentration;
        double term3 = reaction[GET_IDX_2D(i,j,n_tissue,n_tissue)] * dt / (dx/**phi*/);

        return (term1 * term2 + term3 + cOld3[k].concentration);
    }
}

double Solver_Joao::eq_geral (Node *ptr, int i, int graph_id)
{
    double soma_difusao = calcVizinhosDifusao(ptr, graph_id);
    double soma_adveccao = calcVizinhosAdveccao(ptr, graph_id);
    
    double flow = calcula_flow(ptr);
    double radius = calcula_radius(ptr);

    //printf("Soma difusao = %.10lf\n", soma_difusao);                            
    //printf("Soma adveccao = %.10lf\n", soma_adveccao);
    //printf("-------------------------------------------------------\n");
//     if(i == 218)
//         printf("%f\n", M_PI*pow(radius,2));
    
//     return (( D * ( (soma_difusao - (ptr->num_edges * cOld[i].concentration)) / (dx*dx) ) ) 
//              - ( ( ptr->velocity / dx) * ( ( (ptr->num_edges - 1) * cOld[i].concentration ) - soma_adveccao ))) * dt 
//             + cOld[i].concentration;
    
    #ifdef PRESSURE
        if(graph_id == 0)
            return (( D * ( (soma_difusao - (ptr->num_edges * cOld[i].concentration)) / (dx*dx) ) ) \
                    - ( ( (flow/FACTOR_FLOW) / dx) * ( ( (ptr->num_edges - 1) * cOld[i].concentration ) - soma_adveccao )) / (M_PI*pow(radius,2))) * dt 
                    + cOld[i].concentration;
        else if(graph_id == 1)
            return (( D * ( (soma_difusao - (ptr->num_edges * cOld2[i].concentration)) / (dx*dx) ) ) \
                    - ( ( (flow/FACTOR_FLOW) / dx) * ( ( (ptr->num_edges - 1) * cOld2[i].concentration ) - soma_adveccao )) / (M_PI*pow(radius,2))) * dt 
                    + cOld2[i].concentration;
        else if(graph_id == 2)
            return (( D * ( (soma_difusao - (ptr->num_edges * cOld3[i].concentration)) / (dx*dx) ) ) \
                    - ( ( (flow/FACTOR_FLOW) / dx) * ( ( (ptr->num_edges - 1) * cOld3[i].concentration ) - soma_adveccao )) / (M_PI*pow(radius,2))) * dt 
                    + cOld3[i].concentration;
    #else 
    return (( D * ( (soma_difusao - (ptr->num_edges * cOld[i].concentration)) / (dx*dx) ) ) \
             - ( ( 1. / dx) * ( ( (ptr->num_edges - 1) * cOld[i].concentration ) - soma_adveccao ))) * dt 
            + cOld[i].concentration;
    #endif
}

double Solver_Joao::eq_bifurcacao (Node *ptr, int i, int graph_id)
{
    double soma_difusao = calcVizinhosDifusao(ptr, graph_id);
    double soma_adveccao = calcVizinhosAdveccao(ptr, graph_id);
    
    double flow[3];
    double area[3];
    
    double area_media = calcula_area_bifurcacao(ptr);
    

    calcula_elementos_bifurcacao(ptr, flow, area);
    
//     if(i == 61)
//         printf("%f %f %f \n",flow[0], flow[1], flow[2]);    
//     printf("%f \n", area_media);
    //printf("Soma difusao = %.10lf\n", soma_difusao);                            
    //printf("Soma adveccao = %.10lf\n", soma_adveccao);
    //printf("-------------------------------------------------------\n");
    
    #ifdef PRESSURE
    if(graph_id == 0)
        return (( D * ( (soma_difusao - (ptr->num_edges * cOld[i].concentration)) / (dx*dx) ) ) \
                - ( ( ((flow[1]/FACTOR_FLOW)/*/area[1]*/) + ((flow[2]/FACTOR_FLOW)/*/area[2]*/) ) * cOld[i].concentration  - soma_adveccao*((flow[0]/FACTOR_FLOW)/*/area[0]*/) ) / (area_media*dx)) * dt 
                + cOld[i].concentration;
    else if(graph_id == 1)
        return (( D * ( (soma_difusao - (ptr->num_edges * cOld2[i].concentration)) / (dx*dx) ) ) \
                - ( ( (flow[1]/FACTOR_FLOW/*/area[1]*/) + (flow[2]/FACTOR_FLOW/*/area[2]*/) ) * cOld2[i].concentration  - soma_adveccao*(flow[0]/FACTOR_FLOW/*/area[0]*/) ) / (area_media*dx)) * dt 
                + cOld2[i].concentration;
    else if(graph_id == 2)
        return (( D * ( (soma_difusao - (ptr->num_edges * cOld3[i].concentration)) / (dx*dx) ) ) \
                - ( ( (flow[1]/FACTOR_FLOW/*/area[1]*/) + (flow[2]/FACTOR_FLOW/*/area[2]*/) ) * cOld3[i].concentration  - soma_adveccao*(flow[0]/FACTOR_FLOW/*/area[0]*/) ) / (area_media*dx)) * dt 
                + cOld3[i].concentration;
//     return (( D * ( (soma_difusao - (ptr->num_edges * cOld[i].concentration)) / (dx*dx) ) ) 
//              - ( ( flow / dx) * ( ( (ptr->num_edges - 1) * cOld[i].concentration ) - soma_adveccao)) /*/ area_media*/) * dt 
//             + cOld[i].concentration;
    #else 
    return (( D * ( (soma_difusao - (ptr->num_edges * cOld[i].concentration)) / (dx*dx) ) ) \
             - ( ( 1. / dx) * ( ( (ptr->num_edges - 1) * cOld[i].concentration ) - soma_adveccao ))) * dt 
            + cOld[i].concentration;
    #endif
}

double Solver_Joao::calcVizinhosDifusao(Node * noAtual, int graph_id)
{    
    Edge * ptrl = noAtual->edges;
    double soma = 0.;
    int id;
    
    while(ptrl != NULL)
    {
        id = ptrl->id;
        if(graph_id == 0)
            soma += cOld[id].concentration;
        else if(graph_id == 1)
            soma += cOld2[id].concentration;
        else if(graph_id == 2)
            soma += cOld3[id].concentration;
        ptrl = ptrl->next;
    }
    
    return soma;
}

double Solver_Joao::calcVizinhosAdveccao (Node * noAtual, int graph_id)
{
    Edge * ptrl = noAtual->edges;
    double soma = 0.;
    int id;
    
    while(ptrl != NULL)
    {
        id = ptrl->id;
        if (ptrl->entrada == 1)
        {
            if(graph_id == 0)
                soma += cOld[id].concentration;
            else if(graph_id == 1)
                soma += cOld2[id].concentration;
            else if(graph_id == 2)
                soma += cOld3[id].concentration;
        }
        ptrl = ptrl->next;
    }
    
    return soma;
}

double Solver_Joao::calcula_area_bifurcacao (Node * noAtual)
{
    Edge * ptrl = noAtual->edges;
    double area[3];
    int i = 0;
    
    while(ptrl != NULL)
    {
        area[i] = M_PI * pow(ptrl->radius,2); 
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
    return (area[0]+area[1]+area[2])/3;
}

double Solver_Joao::calcula_flow (Node * noAtual)
{
    Edge * ptrl = noAtual->edges;
    
    while(ptrl != NULL)
    {
        if (ptrl->entrada == 1)
            return ptrl->flow;
        
        ptrl = ptrl->next;
    }
    printf("Erro no nó %d ao se calcular o flow\n", noAtual->id);
    return 0.;
}

void Solver_Joao::calcula_elementos_bifurcacao (Node * noAtual, double flow[], double area[])
{
    Edge * ptrl = noAtual->edges;
    int p = 1;
    
    while(ptrl != NULL)
    {
        if (ptrl->entrada == 1)
        {
            flow[0] = ptrl->flow;
            area[0] = M_PI*pow(ptrl->radius,2);
        }
        else
        {
            flow[p] = ptrl->flow;
            area[p] = M_PI*pow(ptrl->radius,2);
//             printf("%d\n",p);
            p++;            
        }
//         p++;
        ptrl = ptrl->next;
    }
//     printf("Erro no nó %d ao se calcular o flow ou o radius da bifurcação \n", noAtual->id);
}

double Solver_Joao::calcula_radius (Node * noAtual)
{
    Edge * ptrl = noAtual->edges;
    
    while(ptrl != NULL)
    {
        if (ptrl->entrada == 1)
            return ptrl->radius;
        
        ptrl = ptrl->next;
    }
    printf("Erro no nó %d ao se calcular o radius\n", noAtual->id);
    return 0.;
}

// double Solver_Joao::calcula_area_bifurcacao (Node * noAtual)
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

void Solver_Joao::nextTimestep ()
{
    int i, j;
    
    // TO DO: Mudar copia para swap de ponteiro
    //swap(&vol[i].yOld,&vol[i].yNew);
    
    int np;
    
    np = g->getTotalNodes();
    for (i = 0; i < np; i++)
        cOld[i].concentration = cNew[i].concentration;
    
    np = g2->getTotalNodes();
    for (i = 0; i < np; i++)
        cOld2[i].concentration = cNew2[i].concentration;
    
    np = g3->getTotalNodes();
    for (i = 0; i < np; i++)
        cOld3[i].concentration = cNew3[i].concentration;

    
    for(i = 0; i < n_tissue; i++)
    {
        for(j = 0; j < n_tissue; j++)
            cOld_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration = cNew_tissue[GET_IDX_2D(i,j,n_tissue,n_tissue)].concentration;
    }
}

void printProgress (int iter, int max_iter)
{
    double progress = iter / (double)max_iter;
    int barWidth = 100;

    cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) 
    {
        if (i < pos) cout << "=";
        else if (i == pos) cout << ">";
        else cout << " ";
    }
    cout << "] " << int(progress * 100.0) << " %\r";
    cout.flush();
}

// double calcNorm (double x1, double y1, double z1, double x2, double y2, double z2)"
// {
//     return sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2));
// }


void Solver_Joao::print ()
{
//     int np = g->getTotalNodes();
//     for (int i = 0; i < np; i++)
//         printf("id = %d -- C = %.10lf\n",i,cOld[i].concentration);
    Edge *ptrl;
    Node *ptr = g->getListNodes();
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
            printf(" -> | Dist = %.4lf | Vel = %.4lf |\n", calcNorm (ptrl->dest->x,ptrl->dest->y,ptrl->dest->z,ptr->x,ptr->y,ptr->z), ptr->velocity);
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

void Solver_Joao::search_proximal_in_map (int j, int u[], int net_id) {
    if (net_id == 0 && joao_grid_to_net_1[j].size() % 2 == 0) {
        u[0] = joao_grid_to_net_1[j][0];
        u[1] = joao_grid_to_net_1[j][1];
    }
    else if (net_id == 1 && joao_grid_to_net_2[j].size() % 2 == 0) {
        u[0] = joao_grid_to_net_2[j][0];
        u[1] = joao_grid_to_net_2[j][1];
    }
    else if (net_id == 2 && joao_grid_to_net_3[j].size() % 2 == 0) {
        u[0] = joao_grid_to_net_3[j][0];
        u[1] = joao_grid_to_net_3[j][1];
    }
}

void Solver_Joao::search_proximal_in_map_v2 (int j, std::vector<std::pair<int,int>> &u_arr, int net_id) {
    if (net_id == 0 && joao_grid_to_net_1[j].size() % 2 == 0) {
        int k = 0;
        int n_link = joao_grid_to_net_1[j].size() / 2;
        for (int i = 0; i < n_link; i++) {
            int ii = joao_grid_to_net_1[j][k];
            int jj = joao_grid_to_net_1[j][k+1];
            k += 2;
            u_arr.push_back(std::make_pair(ii, jj));
        }
    }
    else if (net_id == 1 && joao_grid_to_net_2[j].size() % 2 == 0) {
        int k = 0;
        int n_link = joao_grid_to_net_2[j].size() / 2;
        for (int i = 0; i < n_link; i++) {
            int ii = joao_grid_to_net_2[j][k];
            int jj = joao_grid_to_net_2[j][k+1];
            k += 2;
            u_arr.push_back(std::make_pair(ii, jj));
        }
    }
    else if (net_id == 2 && joao_grid_to_net_3[j].size() % 2 == 0) {
        int k = 0;
        int n_link = joao_grid_to_net_3[j].size() / 2;
        for (int i = 0; i < n_link; i++) {
            int ii = joao_grid_to_net_3[j][k];
            int jj = joao_grid_to_net_3[j][k+1];
            k += 2;
            u_arr.push_back(std::make_pair(ii, jj));
        }
    }
}

void Solver_Joao::createFolders () {
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

void Solver_Joao::calcula_reaction_v2(double * reaction, int k, std::vector<std::pair<int,int>> u_arr, int graph_id) {
    if(graph_id == 0) {
        // Pass through the tissue neighbour cells from terminal 'k'
        for (int i = 0; i < u_arr.size(); i++) {
            int ii = u_arr[i].first;
            int jj = u_arr[i].second;

            if(cOld[k].concentration < cOld_tissue[GET_IDX_2D(ii,jj,n_tissue,n_tissue)].concentration) {
                reaction[GET_IDX_2D(ii,jj,n_tissue,n_tissue)] = 0.;
            }
            else {
                if(!(tecido[GET_IDX_2D(ii,jj,n_tissue,n_tissue)] == 50))
                    reaction[GET_IDX_2D(ii,jj,n_tissue,n_tissue)] = - r * (cOld[k].concentration - cOld_tissue[GET_IDX_2D(ii,jj,n_tissue,n_tissue)].concentration);
                else
                    reaction[GET_IDX_2D(ii,jj,n_tissue,n_tissue)] = - r_problema * (cOld[k].concentration - cOld_tissue[GET_IDX_2D(ii,jj,n_tissue,n_tissue)].concentration);
            }
        }        
    }
    
    else if(graph_id == 1) {
        // Pass through the tissue neighbour cells from terminal 'k'
        for (int i = 0; i < u_arr.size(); i++) {
            int ii = u_arr[i].first;
            int jj = u_arr[i].second;

            if(cOld2[k].concentration < cOld_tissue[GET_IDX_2D(ii,jj,n_tissue,n_tissue)].concentration) {
                reaction[GET_IDX_2D(ii,jj,n_tissue,n_tissue)] = 0.;
            }
            else {
                if(!(tecido[GET_IDX_2D(ii,jj,n_tissue,n_tissue)] == 50))
                    reaction[GET_IDX_2D(ii,jj,n_tissue,n_tissue)] = - r * (cOld2[k].concentration - cOld_tissue[GET_IDX_2D(ii,jj,n_tissue,n_tissue)].concentration);
                else
                    reaction[GET_IDX_2D(ii,jj,n_tissue,n_tissue)] = - r_problema * (cOld2[k].concentration - cOld_tissue[GET_IDX_2D(ii,jj,n_tissue,n_tissue)].concentration);
            }
        }
    }
    else if(graph_id == 2) {
        // Pass through the tissue neighbour cells from terminal 'k'
        for (int i = 0; i < u_arr.size(); i++) {
            int ii = u_arr[i].first;
            int jj = u_arr[i].second;

            if(cOld3[k].concentration < cOld_tissue[GET_IDX_2D(ii,jj,n_tissue,n_tissue)].concentration) {
                reaction[GET_IDX_2D(ii,jj,n_tissue,n_tissue)] = 0.;
            }
            else {
                if(!(tecido[GET_IDX_2D(ii,jj,n_tissue,n_tissue)] == 50))
                    reaction[GET_IDX_2D(ii,jj,n_tissue,n_tissue)] = - r * (cOld3[k].concentration - cOld_tissue[GET_IDX_2D(ii,jj,n_tissue,n_tissue)].concentration);
                else
                    reaction[GET_IDX_2D(ii,jj,n_tissue,n_tissue)] = - r_problema * (cOld3[k].concentration - cOld_tissue[GET_IDX_2D(ii,jj,n_tissue,n_tissue)].concentration);
            }
        }
    }
}

double Solver_Joao::eq_terminal_v2 (Node *ptr, int k, std::vector<std::pair<int,int>> u_arr, double * reaction, int graph_id) {
    int id_viz = ptr->edges->id;
    double flow = calcula_flow(ptr);
    double radius = calcula_radius(ptr);

    #ifdef PRESSURE
//     double term1 = - ptr->velocity * dt / dx;
    double term1 = - ((flow/FACTOR_FLOW) * dt / dx)/(M_PI*pow(radius,2));
    #else 
    double term1 = - 1. * dt / dx;
    #endif
    
    if(graph_id == 0)
    {
        double term2 = cOld[k].concentration - cOld[id_viz].concentration;
        // Compute the mean between the neighbour tissue cells 
        double mean_reaction = 0.0;
        for (int i = 0; i < u_arr.size(); i++) {
            int ii = u_arr[i].first;
            int jj = u_arr[i].second;
            mean_reaction += reaction[GET_IDX_2D(ii,jj,n_tissue,n_tissue)];
        }
        mean_reaction /= (double)u_arr.size();
        double term3 = mean_reaction * dt / (dx/**phi*/);

        return (term1 * term2 + term3 + cOld[k].concentration);
    }
    else if(graph_id == 1)
    {
        double term2 = cOld2[k].concentration - cOld2[id_viz].concentration;
        // Compute the mean between the neighbour tissue cells 
        double mean_reaction = 0.0;
        for (int i = 0; i < u_arr.size(); i++) {
            int ii = u_arr[i].first;
            int jj = u_arr[i].second;
            mean_reaction += reaction[GET_IDX_2D(ii,jj,n_tissue,n_tissue)];
        }
        mean_reaction /= (double)u_arr.size();
        double term3 = mean_reaction * dt / (dx/**phi*/);

        return (term1 * term2 + term3 + cOld2[k].concentration);
    }
    else if(graph_id == 2)
    {
        double term2 = cOld3[k].concentration - cOld3[id_viz].concentration;
        // Compute the mean between the neighbour tissue cells 
        double mean_reaction = 0.0;
        for (int i = 0; i < u_arr.size(); i++) {
            int ii = u_arr[i].first;
            int jj = u_arr[i].second;
            mean_reaction += reaction[GET_IDX_2D(ii,jj,n_tissue,n_tissue)];
        }
        mean_reaction /= (double)u_arr.size();
        double term3 = mean_reaction * dt / (dx/**phi*/);

        return (term1 * term2 + term3 + cOld3[k].concentration);
    }
}