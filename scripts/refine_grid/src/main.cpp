// ===================================================================================================================
// Authors: Lucas Berg
// Paper: ????
// Last change: 04/07/2023
// -------------------------------------------------------------------------------------------------------------------
// Program to refine a 2D tissue grid by of half the input space discretization. 
//===================================================================================================================

#include <iostream>
#include <fstream>
#include "../include/utils.h"

using namespace std;

void write_grid_to_vtk (std::string filename, int *grid, int nx, int ny, double h) {
    std::ofstream arqvtk;

    arqvtk.open(filename.c_str());
    arqvtk << "# vtk DataFile Version 3.0\n";
    arqvtk << "vtk output\n";
    arqvtk << "ASCII\n";
    arqvtk << "DATASET RECTILINEAR_GRID\n";
    arqvtk << "DIMENSIONS " << nx << " " << ny << " 1\n";

    arqvtk << "X_COORDINATES " << nx << " double \n";
    for (int i = 0; i < nx; i++)
        arqvtk << i * h << " ";
    arqvtk << "\n";

    arqvtk << "Y_COORDINATES " << ny << " double \n";
    for (int j = 0; j < ny; j++)
        arqvtk << j * h << " ";
    arqvtk << "\n";

    arqvtk << "Z_COORDINATES 1 double \n";
    arqvtk << "0\n";

    arqvtk << "POINT_DATA " << nx * ny << "\n";
    arqvtk << "FIELD FieldData 1\n";
    arqvtk << "Grid 1 " << nx * ny << " double\n";
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            arqvtk << grid[i*nx+j] << " ";
        }
    }
    arqvtk << "\nMETADATA\n";
    arqvtk << "INFORMATION 0\n\n";

    arqvtk << "\n";
    arqvtk.close();
}

void write_grid_to_dat (std::string filename, int *grid, int nx, int ny, double h) {
    std::ofstream arqdat;
    arqdat.open(filename.c_str());
    arqdat << nx << " " << ny << " " << h << std::endl;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            arqdat << grid[i*nx+j] << " ";
        }
        arqdat << "\n";
    }
    arqdat.close();
}   

int main (int argc, char *argv[]) {
    if (argc-1 != 3) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    
    std::string input_filename = argv[1];
    double input_h = atof(argv[2]);
    std::string output_filename = argv[3];

    // Read input mesh
    FILE *input_file = fopen(input_filename.c_str(), "r");
    int nx, ny;
    int *grid = NULL;
    fscanf(input_file, "%d %d", &nx, &ny);
    grid = (int*)malloc(sizeof(int)*nx*ny);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            fscanf(input_file, "%d", &grid[i*nx+j]);
        }
    }
    fclose(input_file);

    // Refine the grid by half
    int *refined_grid = (int*)malloc(sizeof(int)*(2*nx-1)*(2*ny-1));
    for (int i = 0; i < (2*nx-1); i++) {
        for (int j = 0; j < (2*ny-1); j++) {
            int ii = i/2;
            int jj = j/2;
            refined_grid[i*(2*nx-1)+j] = grid[ii*nx+jj];
        }
    }

    // Original grid
    //write_grid_to_vtk("outputs/vent100x100.vtk", grid, nx, ny, input_h);
    //write_grid_to_dat("outputs/vent100x100.dat", grid, nx, ny, input_h);

    // Refined grid
    //write_grid_to_vtk("outputs/vent50x50_with_nx_ny_h.vtk", refined_grid, (2*nx-1), (2*ny-1), input_h/2.0);
    //write_grid_to_dat("outputs/vent50x50_with_nx_ny_h.dat", refined_grid, (2*nx-1), (2*ny-1), input_h/2.0);

    write_grid_to_vtk("outputs/vent25x25_with_nx_ny_h.vtk", refined_grid, (2*nx-1), (2*ny-1), input_h/2.0);
    write_grid_to_dat(output_filename.c_str(), refined_grid, (2*nx-1), (2*ny-1), input_h/2.0);

    // Free memory
    free(grid);
    free(refined_grid);

    return 0;
}
