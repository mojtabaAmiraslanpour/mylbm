#include <iostream>
#include "src/grid.H"
#include "src/problem.H"

// GENERAL FLOW CONSTANTS
const int q_ = 9;
int lx_ = 400; // number of cells in x-direction
int ly_ = 100; // number of cells in y-direction
double pos_x = lx_/5; // position of the cylinder; (exact
double pos_y = ly_/2+3; // y-symmetry is avoided)
double radius = ly_/10+1; // radius of the cylinder
double uMax = 0.1; // maximum velocity of Poiseuille inflow
double Re = 100; // Reynolds number
double nu = uMax * 2. * radius / Re; // kinematic viscosity
double omega = 1. / (3 * nu+1./2.); // relaxation parameter, assuming deltaT = 1 so C_s = 1/3 p.117
int maxT = 400000; // total number of iterations
int tPlot = 50; // cycles

// Using D2Q9
double w[q_]  = {4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36}; // Weights
int cx[q_] = {0,   1,  0, -1,  0,    1,  -1,  -1,   1};
int cy[q_] = {0,   0,  1,  0, -1,    1,   1,  -1,  -1};

int L = ly_ - 2;
//double ux = 4. * uMax / (L*L) * (L*y - y*y);

int main(int, char**) {

    // Creating the mesh
    grid mesh_(lx_, ly_);

    // Defining the problem
    problem cylinder(mesh_, q_);

    // Initial conditions with macroscopic values
    double rho = 1.0;
    double Ux[ly_] = {};
    //double y_phys = ly_ - 0.5;
    
    for(int i =0; i<ly_; i++){
        Ux[i] = 4. * uMax / (L*L) * (L*i - i*i);
    }
    double Uy = 0.0;
    double test[5] = {};

    // initialization of the population with macroscopic variables
    cylinder.initialize(rho, Ux, Uy);
};