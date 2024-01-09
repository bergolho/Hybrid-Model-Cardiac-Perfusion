#ifndef MODEL_H_
#define MODEL_H_

#include <cstdio>
#include <cstdlib>
#include <string>
#include "../include/solver.h"
#include "../include/solver_joao.h"

using namespace std;

class Model
{
public:
    Model (int argc, char *argv[], int solver_id, int load_map);
    ~Model ();
    void solve (int solver_id, int load_map);
    void error (const char msg[]);
private:
    Solver *sol;
    Solver_Joao *sol_joao;
};

// Function to call the solver from main
void solve_model (int argc, char *argv[], int solver_id, int load_map);

#endif
