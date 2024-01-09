#include "../include/utils.h"

void usage (const char pName[]) {
  printf("%s\n", LINE_1.c_str());
  printf("Usage:> %s <input_mesh_file> <input_h> <output_mesh_file>\n",pName);
  printf("%s\n", LINE_2.c_str());
  printf("<input_mesh_file> = Input mesh with the tissue grid\n");
  printf("<input_h> = Input space discretization\n");
  printf("<output_mesh_file> = Output mesh with the refined tissue grid\n");
  printf("%s\n", LINE_2.c_str());
  printf("Solver Example: %s inputs/vent100x100.dat 0.5 outputs/vent50x50.dat\n",pName);
  printf("%s\n", LINE_1.c_str());
}

void print_progress (int iter, int max_iter) {
    double progress = iter / (double)max_iter;
    int barWidth = 100;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) 
    {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}