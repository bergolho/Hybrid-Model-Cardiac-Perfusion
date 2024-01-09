#ifndef SOLVER_JOAO_H_
#define SOLVER_JOAO_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>        // sqrt(), pow(), M_PI, ...
#include <string>
#include <omp.h>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>

#include "../include/graph.h"

#define GET_IDX_2D(x, y, nx, ny) (ny)*(x) + y
#define SCALE_FACTOR 1.0

using namespace std;
//using namespace Eigen;

// Macros for the Eigen library structures
//typedef SparseMatrix<double> SpMat;SI_health_tissue
//typedef Triplet<double> T;

// Function pointer
typedef double (*Func) (int type, int point, double t, double y[]);

// Estrutura de uma celula
struct Cell_Joao
{
  int type;                                   // 0 = raiz || 1 = terminal || 2 = normal || 3 = bifurcacao
  double concentration;                       // Concentration
  int roi;                                    // 0 = ROI healthy || 1 = ROI ischemic
}typedef Cell_Joao;

struct Derivative
{
  double t;                                   // Time of the maximum derivative
  double value;                               // Derivative value
}typedef Derivative;

struct Velocity
{
  FILE *velocityFile;                         // Reference to the file where the velocity will be stored
  int np;                                     // Number of volumes that the velocity will be calculated
  int id_source;                              // Identifier of the source volume
  int *ids;                                   // Identifier of the sink points (ids[0] -> source)
  double t1;                                  // Time when the despolarization occur on the source volume
  double *t2;                                 // Time when the despolarization occur on the sink volumes
}typedef Velocity;

// Structure for the plot volumes
struct Plot
{
  FILE **plotFile;                            // Reference to the file of the plot volume
  int np;                                     // Number of plot volumes
  int *ids;                                   // Identifier of the plot volume
}typedef Plot;

class Solver_Joao
{
    // Constantes da equacao de transporte
    //static constexpr double BETA = 0.14;
//     static constexpr int n_tissue = 100;
    static constexpr double T_peak = 24.0;
    static constexpr double fator_SI = 130.0;//0.05;
    static constexpr double sigma = 6.0;
    static constexpr double lambda = 0.25;

    //ultima árvore: alpha = 10 para infarto e 8 para isquemia


//    static constexpr int SAVEEACH = 1000;
//    static constexpr double lambda_infarto = 1.0; // referencia para isquemia: 0.25
//    static constexpr double lambda_fibrose = 0.5; // referencia para isquemia: 0.0
//    static constexpr double rf = 0.5; // referencia para isquemia: 0.0
//    static constexpr double decaimento = 0.007;
//    static constexpr double decaimento_infarto = 0.007;//.001;
//    static constexpr double decaimento_fibrose = .0009; // referencia para isquemia: 0.0; referencia para infarto: 0.0009


     static constexpr int SAVEEACH = 1000;
     static constexpr double lambda_infarto = 1.0; // referencia para isquemia: 0.25
     static constexpr double lambda_fibrose = 0.5; // referencia para isquemia: 0.0
     static constexpr double rf = 0.01; // referencia para isquemia: 0.0
     static constexpr double decaimento = 0.007;
     static constexpr double decaimento_infarto = 0.002;
     static constexpr double decaimento_fibrose = 0.001; // referencia para isquemia: 0.0; referencia para infarto: 0.0009
 //     static constexpr double rf = 0.5; // referencia para isquemia: 0.0

    static constexpr double r = 0.01;
    static constexpr double r_problema = 0.01;
    static constexpr double phi = 0.10;
    static constexpr double D = 0.05;             // Difusao
    //static constexpr double h_gadolinio = 0.5*SCALE_FACTOR;      // discretização espacial
    static constexpr double v_1D = 0.060;
    static constexpr int n_1D = 50;
    static constexpr double L_1D = 1.0;
    static constexpr double phi_1D = 1.0;
    static constexpr double h_1D = L_1D/n_1D;
    static constexpr double r3 = - 1.0;
    static constexpr double rins = 0.010;//1.01;
    static constexpr double FACTOR_FLOW = 1.0/16.6667;
public:
    Solver_Joao (int argc, char *argv[], int load_map);
    double calcVizinhosDifusao (Node * noAtual, int graph_id);
    double calcVizinhosAdveccao (Node * noAtual, int graph_id);
    double gaussiana (double t);
    double eq_terminal (Node *ptr, int k, int i, int j, double * reaction, int graph_id);
    double eq_terminal_v2 (Node *ptr, int k, std::vector<std::pair<int,int>> u_arr, double * reaction, int graph_id);
    double eq_geral (Node *ptr, int i, int graph_id);
    double eq_bifurcacao (Node *ptr, int i, int graph_id);
    int number_tissue();
    void solve1D_source(double * C_1D, double V, double adv_1D, double difusao_1D);
    double calcula_fonte(double * C_1D, double t);
    double calcula_fonte_arteria(double * C_1D, int i, double dt, double C);
    double calcula_flow (Node * noAtual);
    void calcula_elementos_bifurcacao (Node * noAtual, double flow[], double area[]);
    double calcula_radius (Node * noAtual);
    double calcula_area_bifurcacao (Node * noAtual);


//     ponto_tecido busca_proximal(Graph *g, int id, ponto_tecido P);
//     void calc_velocity();
    void solve ();
    void print ();
    void error (const char msg[]);
private:
    Graph *g;                           // Graph representing the Purkinje network 1
    Graph *g2;                           // Graph representing the Purkinje network 2
    Graph *g3;                           // Graph representing the Purkinje network 3

    Cell_Joao *cOld;                         // Vetor de volumes de controle no tempo (n)
    Cell_Joao *cNew;                         // Vetor de volumes de controle no tempo (n+1)
    Cell_Joao *cOld2;
    Cell_Joao *cNew2;
    Cell_Joao *cOld3;
    Cell_Joao *cNew3;

    Cell_Joao *cOld_tissue;                  // Vetor de volumes de controle no tecido no tempo (n)
    Cell_Joao *cNew_tissue;                  // Vetor de volumes de controle no tecido no tempo (n+1)

    Cell_Joao *cOld_fibrose;                  // Vetor de volumes de controle na fibrose no tempo (n)
    Cell_Joao *cNew_fibrose;                  // Vetor de volumes de controle na fibrose no tempo (n+1)

    int * tecido;
    Derivative *dvdt;                   // Vector with the maximum derivative for each volume
    Plot *plot;                         // Vector with the plot ids
    Velocity *vel;                      // Vector with the propagation velocity of the plot ids
    Func *func;                         // Vector of function of the celullar model
    int M;                              // Number of timesteps in time
    double dt;                          // Size timestep in time
    double tmax;                        // Maximum time of the simulation
    double dx;                          // Size timestep in space
    double h_gadolinio;
    string network_filename_1;          // Mesh filename 1
    string network_filename_2;          // Mesh filename 2
    string network_filename_3;          // Mesh filename 3
    string steady_filename;             // Input Steady-State filename
    string plot_filename;               // Plot id filename
    string tissue_filename;             // Tissue (0,128 or 255) filename
    string output_dir;                  // Output directory to store the simulation
    int n_tissue;    // number of points at the tissue file
    int n_threads;

//     Node* searchNode (int id);
    double alfa;                        // Parameter: R_pmj * Vol_pmj
    double d1;                          // Parameter: d1
    double BETA;                        // Surface / Volume ratio
    double SIGMA;                       // Conductivity Gap + Citoplasm
    void busca_proximal(int id, int p[], int graph_id);
    void setSensibilityParam (int argc, char *argv[]);
    void setTypeCell ();
    void setTypeTissue ();
    void setTypeTissue_v2 ();
    void setControlVolumes ();
    void setFunctions ();
    void setInitCondFromFile ();
    void setVelocityPoints ();
    void setDerivative ();
    void setPlot ();
    void setTerm ();
    void leGeometria(string filename);
    double integral_extra ();
    double integral_intra ();
    void calcula_reaction(double * reaction, int k, int i, int j, int graph_id);
    void calcula_reaction_v2(double * reaction, int k, std::vector<std::pair<int,int>> u_arr, int graph_id);
    void writeVTKFile_static ();
    //void setMatrix (SpMat &a);
    //void setMatrix2 (SpMat &a);
    //void assembleLoadVector (VectorXd &b);
    //void moveVstar (const VectorXd vm);
    bool isConnToPMJ (Edge *e);
    void solveODE (double t);
    void nextTimestep ();
    void calcMaxDerivative (double t);
    void calcVelocity ();
    void writeVTKFile_tissue (char * tissue, int nx, int ny);
    void writeVTKFile (int iter, int graph_id);
    void writePlotData (double t);
    double * contorno_perm_het(double * p, int * tecido);
    double * perm_het (double * p, int * tecido);
    double signal_intensity_healthy (int k, int * tecido, Cell_Joao * S, double * peso, double SI_health_tissue);
    double signal_intensity_ischemic (int k, int * tecido, Cell_Joao * S, double * peso, double SI_ischemic_tissue);
    double signal_intensity(int k, int * tecido, Cell_Joao * S, Cell_Joao * C_fibrose, double * peso, double SI_fibrose, double SI_ischemic, double SI);
    double signal_intensity_fibrose(int k, int * tecido, Cell_Joao * S, Cell_Joao * C_fibrose, double * peso, double SI_fibrose, double SI_ischemic, double SI);
    void search_region(int id, int p[], int graph_id, double * peso);
    void calc_Concentracao_Vizinhos (Node * noAtual, int graph_id, double p[]);
    void search_proximal_in_map (int j, int u[], int net_id);
    void search_proximal_in_map_v2 (int j, std::vector<std::pair<int,int>> &u_arr, int net_id);
    void createFolders ();
};



void swap (double **a, double **b);
void printProgress (int iter, int max_iter);

#endif
