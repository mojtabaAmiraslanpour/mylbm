#include <iostream>
#include "src/grid.h"
#include "src/problem.h"

// GENERAL FLOW CONSTANTS
const int q_ = 9;
const int lx_ = 400; // number of cells in x-direction
const int ly_ = 100; // number of cells in y-direction
double pos_x_ = lx_/5; // position of the cylinder; (exact
double pos_y_ = ly_/2+3; // y-symmetry is avoided)
double radius_ = ly_/10+1; // radius of the cylinder
double Re_ = 100; // Reynolds number
double uMax_ = 0.1; // maximum velocity of Poiseuille inflow
double nu_ = uMax_ * 2. * radius_ / Re_; // kinematic viscosity
double omega_ = 1. / (3 * nu_+1./2.); // relaxation parameter, assuming deltaT = 1 so C_s = 1/3 p.117
int maxT_ = 400000; // total number of iterations
int tPlot_ = 50; // cycles

// Using D2Q9
double w_[q_]  = {4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36}; // Weights
int cx_[q_] = {0,   1,  0, -1,  0,    1,  -1,  -1,   1};
int cy_[q_] = {0,   0,  1,  0, -1,    1,   1,  -1,  -1};

int main(int, char**) {
    // Creating the mesh
    grid mesh_(lx_, ly_);

    // Defining the problem
    problem cylinder(mesh_, q_, cx_, cy_, w_);

    // initialization of the population with macroscopic variables
    cylinder.initialize(uMax_);
};