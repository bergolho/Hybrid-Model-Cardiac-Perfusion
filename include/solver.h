#ifndef SOLVER_H_
#define SOLVER_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>        // sqrt(), pow(), M_PI, ...
#include <string>
#include <cstring>
#include <map>
#include <omp.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../include/graph.h"
#include "../include/utils.h"

#define GET_IDX_2D(x, y, nx, ny) (ny)*(x) + y
#define SCALE_FACTOR 1.0

using namespace std;

// Function pointer
typedef double (*Func) (int type, int point, double t, double y[]);

struct ponto_tecido
{
  int type;                                    // 1 = terminal || 0 = não-terminal
  double x,y,z;
}typedef ponto_tecido;

// Estrutura de uma celula
struct Cell
{
  int type;                                   // 0 = raiz || 1 = terminal || 2 = normal || 3 = bifurcacao
  double concentration;                       // Concentration
  int roi;                                    // 0 = ROI healthy || 1 = ROI ischemic
}typedef Cell;

class Solver {
    
// Constantes da equacao de transporte
     // static constexpr int n_tissue = 100;
     static constexpr double T_peak = 24.0;
     static constexpr double fator_SI = 130.0;//0.05;
     static constexpr double sigma = 6.0;
     static constexpr double lambda = 0.25;

//   ultima árvore: alpha = 10 para infarto e 8 para isquemia

//   static constexpr int SAVEEACH = 1000;
//   static constexpr double lambda_infarto = 1.0; // referencia para isquemia: 0.25
//   static constexpr double lambda_fibrose = 0.5; // referencia para isquemia: 0.0
//   static constexpr double rf = 0.5; // referencia para isquemia: 0.0
//   static constexpr double decaimento = 0.007;
//   static constexpr double decaimento_infarto = 0.007;//.001;
//   static constexpr double decaimento_fibrose = .0009; // referencia para isquemia: 0.0; referencia para infarto: 0.0009

     static constexpr int SAVEEACH = 1000;
     static constexpr double lambda_infarto = 1.0; // referencia para isquemia: 0.25
     static constexpr double lambda_fibrose = 0.5; // referencia para isquemia: 0.0
     static constexpr double rf = 0.01; // referencia para isquemia: 0.0
     static constexpr double decaimento = 0.007;
     static constexpr double decaimento_infarto = 0.002;
     static constexpr double decaimento_fibrose = 0.001; // referencia para isquemia: 0.0; referencia para infarto: 0.0009
 //  static constexpr double rf = 0.5; // referencia para isquemia: 0.0

    static constexpr double r = 0.01;
    static constexpr double r_problema = 0.01;

    static constexpr double phi = 0.10; // Porosidade
    static constexpr double D = 0.05; // Difusao do sistema recirculatorio
    //static constexpr double h_gadolinio = 0.5*SCALE_FACTOR; // discretização espacial
    static constexpr double v_1D = 0.060; // Velocidade do sistema recirculatorio
    static constexpr int n_1D = 50;
    static constexpr double L_1D = 1.0;
    static constexpr double phi_1D = 1.0;
    static constexpr double h_1D = L_1D/n_1D;
    static constexpr double r3 = - 1.0;
    static constexpr double rins = 0.010;//1.01;
    static constexpr double FACTOR_FLOW = 1.0/16.6667;
public:
    Solver (int argc, char *argv[], int load_map);
    ~Solver ();
    double calc_vizinhos_difusao (Node * noAtual, int graph_id);
    double calc_vizinhos_adveccao (Node * noAtual, int graph_id);
    double gaussiana (double t);
    double eq_terminal (Node *ptr, int k, int i, int j, double * reaction, int graph_id);
    double eq_geral (Node *ptr, int i, int graph_id);
    double eq_bifurcacao (Node *ptr, int i, int graph_id);
    int number_tissue();
    void solve1D_source(double * C_1D, double V, double adv_1D, double difusao_1D);
    void solve1D_source_v2 (double *C_1D, double V, double adv_1D, double difusao_1D);
    double calcula_fonte(double * C_1D, double t);
    double calcula_fonte_arteria(double * C_1D, int i, double dt, double C);
    double calcula_flow (Node * noAtual);
    void calcula_elementos_bifurcacao (Node * noAtual, double flow[], double area[]);
    double calcula_radius (Node * noAtual);
    double calcula_area_bifurcacao (Node * noAtual);


//     ponto_tecido busca_proximal(Graph *g, int id, ponto_tecido P);
//     void calc_velocity();
    void solve (int load_map);
    void print ();
    void print_parameters ();
    void print_cells_network (int net_id);
    void print_cells_tissue ();
    void print_graph ();
    void error (const char msg[]);
private:
    Graph *g1;                            // Graph representing the CCO network 1
    int *cell_types_network_1;            // 0 = raiz, 1 = terminal, 2 = normal, 3 = bifurcacao
    int *roi_network_1;                   // 0 = healthy, 1 = ischemia
    double *concentration_old_network_1;  // Vetor de volumes de controle no tempo (n)
    double *concentration_new_network_1;  // Vetor de volumes de controle no tempo (n+1)
    
    Graph *g2;                            // Graph representing the CCO network 2
    int *cell_types_network_2;            // 0 = raiz, 1 = terminal, 2 = normal, 3 = bifurcacao
    int *roi_network_2;                   // 0 = healthy, 1 = ischemia
    double *concentration_old_network_2;  // Vetor de volumes de controle no tempo (n)
    double *concentration_new_network_2;  // Vetor de volumes de controle no tempo (n+1)

    Graph *g3;                            // Graph representing the CCO network 3
    int *cell_types_network_3;            // 0 = raiz, 1 = terminal, 2 = normal, 3 = bifurcacao
    int *roi_network_3;                   // 0 = healthy, 1 = ischemia
    double *concentration_old_network_3;  // Vetor de volumes de controle no tempo (n)
    double *concentration_new_network_3;  // Vetor de volumes de controle no tempo (n+1)

    int *cell_types_tissue;               // 0 = normal, 1 = ligado a arvore
    double *concentration_old_tissue;     // Vetor de volumes de controle no tempo (n)
    double *concentration_new_tissue;     // Vetor de volumes de controle no tempo (n+1)
    double *concentration_old_fibrosis;   // Vetor de volumes de controle no tempo (n)
    double *concentration_new_fibrosis;   // Vetor de volumes de controle no tempo (n+1)

    double *C_1D_old;                       // Vetor da solucao da fonte 1D (n)
    double *C_1D_new;                       // Vetor da solucao da fonte 1D (n+1)

    bool load_map_from_file;
    string map_filename;

  // OLD CODE
    Cell *cOld;                         // Vetor de volumes de controle no tempo (n)
    Cell *cNew;                         // Vetor de volumes de controle no tempo (n+1)
    Cell *cOld2;
    Cell *cNew2;
    Cell *cOld3;
    Cell *cNew3;

    Cell *cOld_tissue;                  // Vetor de volumes de controle no tecido no tempo (n)
    Cell *cNew_tissue;                  // Vetor de volumes de controle no tecido no tempo (n+1)

    Cell *cOld_fibrose;                  // Vetor de volumes de controle na fibrose no tempo (n)
    Cell *cNew_fibrose;                  // Vetor de volumes de controle na fibrose no tempo (n+1)
  // OLD CODE

    //int *tecido;                        // Grid of the tissue elements
    int *grid;                          // Grid of the tissue elements
    ponto_tecido *PT;
    int M;                              // Number of timesteps in time
    double dt;                          // Size timestep in time
    double tmax;                        // Maximum time of the simulation
    double dx;                          // Size timestep in space
    string network_filename_1;          // Mesh filename 1
    string network_filename_2;          // Mesh filename 2
    string network_filename_3;          // Mesh filename 3
    string tissue_filename;             // Tissue (0,128 or 255) filename
    string output_dir;                  // Output directory to store the simulation
    int n_tissue;                       // number of points at the tissue file
    int n_threads;                      // Number of OpenMP threads
    double h_gadolinio;                  // Space discretization

//  Node* searchNode (int id);
    void create_folders ();
    void allocate_memory_for_networks ();
    void set_type_cells_networks ();
    void set_type_cells_using_network (Graph *g, int *cell_type, double *c_old, double *c_new);
    void set_type_cells_tissue ();
    void set_type_cells_tissue_using_network (int net_id);
    void allocate_memory_for_tissue ();
    void allocate_memory_for_network_cells(int net_id, int num_cells);
    void busca_proximal(int id, int p[], int graph_id);
    void read_and_allocate_tissue_geometry (std::string filename);
    void set_initial_conditions_tissue (double *reaction, double *peso);
    int count_number_active_tissue_cells();
    double integral_extra (int N);
    double integral_intra ();
    double integral_intra_v2 ();
    void calcula_reaction(double * reaction, int k, int i, int j, int graph_id);
    void next_timestep ();
    void post_processing (double *SI_health_vec, double *SI_ischemic_vec, double *SI_fibrose_vec);
    void write_vtk_file_networks (Graph *g, int *cell_type, int *roi, double *c_new, int net_id, int iter);
    void write_vtk_file_tissue (std::string filename, int nx, int ny);
    void contorno_perm_het(double * p, int * tecido);
    double * perm_het (double * p, int * tecido);
    void calculate_signal_intensity (int i, double *peso, double SI_ischemic, double &SI_health_tissue, double &SI_fibrose, double &SI, \
                                         double *SI_ischemic_vec, double *SI_fibrose_vec, double *SI_health_vec);
    double signal_intensity_healthy(int k, int *grid, Cell * S, double *peso, double SI_health);
    double signal_intensity_healthy_v2(int k, int *grid, double *c_new, double *peso, double SI_health);
    double signal_intensity_ischemic (int k, int * tecido, Cell * S, double * peso, double SI_ischemic_tissue);
    double signal_intensity_ischemic_v2 (int k, int * tecido, double *c_new, double * peso, double SI_ischemic_tissue);
    double signal_intensity(int k, int * tecido, Cell * S, Cell * C_fibrose, double * peso, double SI_fibrose, double SI_ischemic, double SI);
    double signal_intensity_v2(int k, int * tecido, double *c_new, double *c_new_fibrosis, double * peso, double SI_fibrose, double SI_ischemic, double SI);
    double signal_intensity_fibrose(int k, int * tecido, Cell * S, Cell * C_fibrose, double * peso, double SI_fibrose, double SI_ischemic, double SI);
    double signal_intensity_fibrose_v2(int k, int * tecido, double *c_new, double *c_new_fibrose, double * peso, double SI_fibrose, double SI_ischemic, double SI);
    void search_region(int id, int p[], int graph_id, double * peso);
    void calc_concentracao_vizinhos (Node * noAtual, int graph_id, double p[]);
    void solve_concentration_for_network(Graph *g, int net_id, double *c_old, double *c_new, \
                                                       int *cell_types, int *roi, \
                                                       double *reaction, double fonte_new, \
                                                       double &SI_health_intra);
    void solve_concentration_for_tissue (double *reaction);
    void search_proximal_in_map (int j, int u[], int net_id);

    void parse_input_arguments (int argc, char *argv[]);
    void parse_input_arguments_v2 (int argc, char *argv[], int load_map);
};

#endif
