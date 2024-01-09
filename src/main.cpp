// ===================================================================================================================
// Authors: Lucas Berg and Joao Rafael
// Paper: A Hybrid Model for Cardiac Perfusion: Coupling a Discrete Coronary Arterial Tree Model with a Continuous Porous-Media Flow Model of the Myocardium
// 		Link to paper :> https://doi.org/10.3390/e25081229
// Last change: 04/07/2023
// -------------------------------------------------------------------------------------------------------------------
// Program to solve the advection-reaction equation associated to the gadolino propagation coming from 3
// different vascular trees, generated via the Constrained Constructive Optimization (CCO) algorithm, 
// in a 2D tissue with fibrosis.
// --------------------------------------------------------------------------------------------------------
// This code has two versions one before the revision (Solver_Joao) and another for the revision (Solver).
//===================================================================================================================

#include <cstdio>
#include "../include/timer.h"
#include "../include/model.h"
#include "../include/utils.h"

using namespace std;

// TODO: Move this to a different location (input parameters ?)
const int SOLVER_TO_USE = 1;        // 1 = Lucas, 2 = Joao Rafael
const int LOAD_MAP_FROM_FILE = 1;   // 1 = No load, 2 = Load map file

int main (int argc, char *argv[]) {
  
  if (argc-1 != 7) {
    usage(argv[0]);
    return EXIT_FAILURE;
  }
  else {
    double start, finish, elapsed;

    GET_TIME(start);
    
    // Option 1 (explicit allocation and deallocation)
    //Model *model = new Model(argc, argv);
    //model->solve();
    //delete model;

    // Option 2 (allocation and deallocation is hidden)
    solve_model(argc,argv, SOLVER_TO_USE, LOAD_MAP_FROM_FILE);
    
    GET_TIME(finish);
    elapsed = finish - start;
  
    printf("%s\n", LINE_1.c_str());
    printf("[!] Tempo gasto = %.10lf segundos\n", elapsed);
    printf("%s\n", LINE_1.c_str());
    
    return EXIT_SUCCESS;
  }
}
