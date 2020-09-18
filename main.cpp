#include <iostream>
#include "grid.H"
#include "problem.H"

// GENERAL FLOW CONSTANTS
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
double w[9]  = {4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36}; // Weights
int cx[9] = {0,   1,  0, -1,  0,    1,  -1,  -1,   1};
int cy[9] = {0,   0,  1,  0, -1,    1,   1,  -1,  -1};

int L = ly_ - 2;
//double ux = 4. * uMax / (L*L) * (L*y - y*y);

int main(int, char**) {

    grid mesh_(lx_,ly_);
    double f;
    problem cylinder(mesh_, f);

    cylinder.initialize();
}