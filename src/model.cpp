#include "../include/model.h"

void solve_model (int argc, char *argv[], int solver_id, int load_map) {
    Model *model = new Model(argc, argv, solver_id, load_map);
    model->solve(solver_id, load_map);
    delete model;
}

Model::Model (int argc, char *argv[], int solver_id, int load_map) {
    if (solver_id == 1) {
        printf("[!] Using Lucas`s solver!\n");
        sol = new Solver(argc,argv,load_map);
    }
    else {
        printf("[!] Using Joao`s solver!\n");
        sol_joao = new Solver_Joao(argc,argv,load_map);
    }
}

Model::~Model () {
    if (sol)
        delete sol;
    if (sol_joao)
        delete sol_joao;
}

void Model::solve (int solver_id, int load_map) {
    if (solver_id == 1) {
        sol->solve(load_map);
    }
    else {
        sol_joao->solve();
    }
}

void Model::error (const char msg[]) {
    printf("[-] ERROR on Model ! %s !\n",msg);
    exit(EXIT_FAILURE);
}
