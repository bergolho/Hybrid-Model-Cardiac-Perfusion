#include "../include/utils.h"

void usage (const char pName[]) {
  printf("%s\n", LINE_1.c_str());
  printf("Usage:> %s <dt> <tmax> <net_mesh_file_1> <net_mesh_file_2> <net_mesh_file_3> <tissue_file> <output_dir>\n",pName);
  //printf("Usage:> %s <dt> <tmax> <net_mesh_file_1> <net_mesh_file_2> <net_mesh_file_3> <tissue_file> <output_dir> [input_graph_to_tissue_map]\n",pName);
  printf("%s\n", LINE_2.c_str());
  printf("<dt> = Size of the time discretization\n");
  printf("<tmax> = Maximum simulation time\n");
  printf("<net_mesh_file_1> = First network file with the mesh points and conections\n");
  printf("<net_mesh_file_2> = Second network file with the mesh points and conections\n");
  printf("<net_mesh_file_3> = Third network file with the mesh points and conections\n");
  printf("<tissue_file> = file that has the information of the tissue map\n");
  printf("<output_dir> = Output directory to store the simulation files\n");
  //printf("[input_graph_to_tissue_map] = File with a mapping of the tissue grid\n\tspecifying the index of the closest tree terminal [OPTIONAL]\n");
  printf("%s\n", LINE_2.c_str());
  printf("Solver Example: %s 0.001 100 meshes/arvore_1.msh meshes/arvore_2.msh meshes/arvore_3.msh meshes/new_mesh_format/vent100x100_with_nx_ny_h.dat vent100x100\n",pName);
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