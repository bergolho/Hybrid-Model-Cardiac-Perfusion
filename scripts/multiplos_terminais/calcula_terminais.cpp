#include <cmath>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <sys/time.h>

using namespace std;

#define real double
#define SAVEEACH 1000

int n = 100;
real L = 50.; // milimetros
real h = L/n; //

int n_1D = 50;
real L_1D = 50.;//4300.; // milimetros
real h_1D = L_1D/n_1D; //
real phi = 0.10; // porosidade do meio (adimensional)
real phi_1D = 1.0; // porosidade do meio (adimensional)
real tempo = 60.0; // tempo de perfusão (segundos)
real dt; // passo de tempo
int t; // número de iterações
real r3 = 1.2;
real r = 0.03; // constante de reação do intravascular pro extravascular
real rf = 1.0; // constante de reação do extravascular pra fibrose.

/* INFARTO */
real r_problema = 0.025; //REFERÊNCIA: 0.025
real decaimento = 0.002;
real decaimento_infarto = 0.0009; //menor do que o decaimento nas outras regiões
real decaimento_fibrose = 0.0007; // referencia = 0.0007
real lambda = 0.25; // normal regions
real lambda_infarto = 1.0;
//real lambda_fibrose = 0.; // BIDOMÍNIO
 real lambda_fibrose = 0.5; // TRIDOMÍNIO. REFERENCIA: 0.50
real D11_infarto = 0.1e-7;
real D22_infarto = 0.5e-7;

/* TRIDOMÍNIO 1 */
//real r_problema = 0.015;
//real decaimento = 0.002;
//real decaimento_infarto = 0.0009; //menor do que o decaimento nas outras regiões
//real decaimento_fibrose = 0.0001; // referencia = 0.0007
//real lambda = 0.25; // normal regions
//real lambda_infarto = 1.0;
////real lambda_fibrose = 0.; // BIDOMÍNIO
// real lambda_fibrose = 1.0; // TRIDOMÍNIO. REFERENCIA: 0.50
//real D11_infarto = 0.1e-7;
//real D22_infarto = 0.5e-7;

/* ISQUEMIA */
//real r_problema = 0.03;
//real decaimento = 0.002;
//real decaimento_infarto = 0.002;
//real decaimento_fibrose = 0.0;
//real lambda = 0.25;
//// real lambda_infarto = 0.50; // isquemia
//real lambda_infarto = 0.25; //normal
//real lambda_fibrose = 0.0;
//real D11_infarto = 0.1e-3;
//real D22_infarto = 0.5e-3;

/* BIDOMÍNIO 1 */
//real r_problema = 0.025;
//real decaimento = 0.002;
//real decaimento_infarto = 0.0009;
//real decaimento_fibrose = 0.0007;
//real lambda = 0.25;
//// real lambda_infarto = 0.; // isquemia
//real lambda_infarto = 1.0; //normal
//real lambda_fibrose = 0.0;
//real D11_infarto = 0.1e-3;
//real D22_infarto = 0.5e-3;

/* BIDOMÍNIO 2 */
//real r_problema = 0.025;
//real decaimento = 0.002;
//real decaimento_infarto = 0.0000009;
//real decaimento_fibrose = 0.0007;
//real lambda = 0.25;
//// real lambda_infarto = 0.; // isquemia
//real lambda_infarto = 1.0; //normal
//real lambda_fibrose = 0.0;
//real D11_infarto = 0.1e-3;
//real D22_infarto = 0.5e-3;

real rins = 0.022;
real pi = 3.14159265;
real T_peak = 25.;
real K11 = 1.50;
real K22 = 0.15;
real D11 = 0.1e-3;
real D22 = 0.5e-3;
real K11_infarto = 1.50;
real K22_infarto = 0.15;
real sigma = 7.0;
int numero_de_fontes = 5;
real fluxo_pressao = 0.;

#define GET_IDX_3D(z, x, y, nz, nx, ny) (ny)*(x) + (nx)*(ny)*(z) + y
#define GET_IDX_2D(x, y, nx, ny) (ny)*(x) + y

int * leGeometria(int * tecido)
{
  int p=0;
  std::fstream fileIn;
  char nome[100];
  sprintf(nome, "vent%dx%dInfarto.dat",n, n);

  fileIn.open(nome, std::ios_base::in);
    for(int j=0;j<n;j++)
    {
        for(int i=0;i<n;i++)
        {
            fileIn>>tecido[GET_IDX_2D(j,i,n,n)];
            if (tecido[GET_IDX_2D(j,i,n,n)] == 128 || tecido[GET_IDX_2D(j,i,n,n)] == 50)
                p++;
        }
    }
    fileIn.close();
    return tecido;
}

void save_vtk(char * nome, int nx, int ny, real * v, int vswap)
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
    arqvtk << i*h << " ";
    arqvtk << "\n";

    arqvtk << "Y_COORDINATES " << ny << " double \n";
    for(int j=0; j<ny; j++)
    arqvtk << j*h << " ";
    arqvtk << "\n";

    arqvtk << "Z_COORDINATES 1 double \n";
    arqvtk << "0\n";

    arqvtk << "POINT_DATA " << nx*ny << "\n";
    arqvtk << "FIELD FieldData 1\n";
//   arqvtk << "Pressure(kPa) 1 " << nx*ny << " double \n";
    arqvtk << "C 1 " << nx*ny << " double \n";

    for(int j=0; j<ny; j++)
        for(int i=0; i<nx; i++)
            {
                arqvtk << v[GET_IDX_3D(vswap,j,i,2,nx,ny)] << " ";
            }

//    arqvtk << "VECTORS vectors double \n";
//    for(int i = 0; i < nx; i++)
//        for(int j=0; j<ny; j++)
//            {
//                arqvtk << vxGET_IDX_2D(i,j,nx,ny) << " "<< vyGET_IDX_2D(i,j,nx,ny) << " "<< 0.0 << "    ";
//            }

    arqvtk << "\n";
    arqvtk.close();
}

void save_press(char * nome, int nx, int ny, real * v)
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
    arqvtk << i*h << " ";
    arqvtk << "\n";

    arqvtk << "Y_COORDINATES " << ny << " double \n";
    for(int j=0; j<ny; j++)
    arqvtk << j*h << " ";
    arqvtk << "\n";

    arqvtk << "Z_COORDINATES 1 double \n";
    arqvtk << "0\n";

    arqvtk << "POINT_DATA " << nx*ny << "\n";
    arqvtk << "FIELD FieldData 1\n";
//   arqvtk << "Pressure(kPa) 1 " << nx*ny << " real \n";
    arqvtk << "Pressure 1 " << nx*ny << " double \n";

    for(int j=0; j<ny; j++)
        for(int i=0; i<nx; i++)
            {
                arqvtk << v[GET_IDX_2D(j,i,nx,ny)] << " ";
            }

//    arqvtk << "VECTORS vectors real \n";
//    for(int i = 0; i < nx; i++)
//        for(int j=0; j<ny; j++)
//            {
//                arqvtk << vxGET_IDX_2D(i,j,nx,ny) << " "<< vyGET_IDX_2D(i,j,nx,ny) << " "<< 0.0 << "    ";
//            }

    arqvtk << "\n";
    arqvtk.close();
}


int main (int argc, char *argv[])
{
    n  = atoi(argv[1]);
    int i, j, m = 0;
    char terminal_char[50];
    sprintf(terminal_char, "terminal_char.vtk");
    
    real * terminal = new real[n*n];
    int * tecido = new int[n*n];
    tecido = leGeometria(tecido);

    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if(tecido[GET_IDX_2D(i,j,n,n)] == 128 || tecido[GET_IDX_2D(i,j,n,n)] == 50)
            {
                terminal[GET_IDX_2D(i,j,n,n)] = 0.;
            }

            else
                terminal[GET_IDX_2D(i,j,n,n)] = -1.;
        }
    }
    
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if(tecido[GET_IDX_2D(i,j,n,n)] == 128 || tecido[GET_IDX_2D(i,j,n,n)] == 50)
            {
                    terminal[GET_IDX_2D(i,j,n,n)] = 10.;
                    printf("%d %.3f %.3f %.3f \n", m, j*h, i*h, 0.0);
                    m++;
            }
            j++;
        }
        i++;
    }
    save_press(terminal_char, n, n, terminal);
    
    return 0;
}


























