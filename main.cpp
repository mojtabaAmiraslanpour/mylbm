#include <iostream>
#include "src/grid.h"
#include "src/problem.h"

// GENERAL FLOW CONSTANTS
const int q_ = 9;
const int lx_ = 40; // number of cells in x-direction
const int ly_ = 10; // number of cells in y-direction
double pos_x_ = lx_/5; // position of the cylinder; (exact
double pos_y_ = ly_/2+3; // y-symmetry is avoided)
double radius_ = ly_/10+1; // radius of the cylinder
double Re_ = 100; // Reynolds number
double uMax_ = 0.1; // maximum velocity of Poiseuille inflow
double nu_ = uMax_ * 2. * radius_ / Re_; // kinematic viscosity
double omega_ = 1. / (3 * nu_+1./2.); // relaxation parameter, assuming deltaT = 1 so C_s = 1/3 p.117
int iter_ = 400000; // total number of iterations
int tPlot_ = 50; // cycles

// Using D2Q9
double w_[q_]  = {4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36}; // Weights
int cx_[q_] = {0,   1,  0, -1,  0,    1,  -1,  -1,   1};
int cy_[q_] = {0,   0,  1,  0, -1,    1,   1,  -1,  -1};

int main(int, char**) {

    int L = ly_ - 2;

    // Initial conditions with macroscopic values -----------------------------
    // Defining rho
    double** rho_ = new double* [lx_];
    for (int i = 0; i < lx_; i++) {
        rho_[i] = new double[ly_];
        for (int j = 0; j < ly_; j++) {
            rho_[i][j] = 1.0;
        }
    }

    // Defining Ux
    //double y_phys = ly_ - 0.5;
    double** Ux_ = new double* [lx_];
    for (int i = 0; i < lx_; i++) {
        Ux_[i] = new double[ly_];
        for (int j = 0; j < ly_; j++) {
            Ux_[i][j] = 4. * uMax_ / (L * L) * (L * j - j * i); // Not used y_phys
        }
    }

    // Defining Uy
    double** Uy_ = new double* [lx_];
    for (int i = 0; i < lx_; i++) {
        Uy_[i] = new double[ly_];
        for (int j = 0; j < ly_; j++) {
            Uy_[i][j] = 0;
        }
    }

    // Defining population array f(q,lx,ly)
    double*** f_;
    f_ = new double** [q_];
    for (int i = 0; i < q_; i++) {
        f_[i] = new double* [lx_];
        for (int j = 0; j < lx_; j++) {
            f_[i][j] = new double[ly_];
            for (int k = 0; k < ly_; k++) {
                f_[i][j][k] = 0;
            }
        }
    }

    // Creating the mesh
    grid mesh_(lx_, ly_);

    // Defining the problem
    problem cylinder(mesh_, q_, cx_, cy_, w_);
    
    // initialization of the population with macroscopic variables
    cylinder.initialize(uMax_, f_, rho_, Ux_, Uy_);

    for (int i = 0; i < q_; i++) {
        for (int j = 0; j < lx_; j++) {
            for (int k = 0; k < ly_; k++) {
                if (f_[i][j][k] != 0) {
                    cout << f_[i][j][k] << endl;
                }
            }
        }
    }

    for (int iter = 0; iter < iter_; iter++) {
        
    }

    //delete[] array;
    delete[] f_;

    return 0;
};