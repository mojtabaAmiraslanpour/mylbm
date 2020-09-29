#include <iostream>
#include "src/grid.h"
#include "src/problem.h"
#include <math.h>

// GENERAL FLOW CONSTANTS
const int q_ = 9;
const int lx_ = 2; // number of cells in x-direction
const int ly_ = 6; // number of cells in y-direction
double pos_x_ = lx_/5; // position of the cylinder; (exact
double pos_y_ = ly_/2+3; // y-symmetry is avoided)
double radius_ = ly_/10+1; // radius of the cylinder
double Re_ = 100; // Reynolds number
double uMax_ = 0.1; // maximum velocity of Poiseuille inflow
double nu_ = uMax_ * 2. * radius_ / Re_; // kinematic viscosity
double omega_ = 1. / (3 * nu_+1./2.); // relaxation parameter, assuming deltaT = 1 so C_s = 1/3 p.117
int iter_ = 400; // total number of iterations
int tPlot_ = 50; // cycles

// Using D2Q9
double w_[q_]  = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
 1.0/36, 1.0/36, 1.0/36, 1.0/36}; // Weights
int cx_[q_] = {0,   1,  0, -1,  0,    1,  -1,  -1,   1};
int cy_[q_] = {0,   0,  1,  0, -1,    1,   1,  -1,  -1};


int main(int, char**) {

    // Initial conditions with macroscopic values
    // Defining rho
    double** rho_ = new double* [lx_];
    for (int i = 0; i < lx_; i++) {
        rho_[i] = new double[ly_];
        for (int j = 0; j < ly_; j++) {
            rho_[i][j] = 1.0;
        }
    }

    // Defining Ux
    double y_;
    int L = ly_-2;
    double** Ux_ = new double* [lx_];
    for (int i = 0; i < lx_; i++) {
        Ux_[i] = new double[ly_];
        for (int j = 0; j < ly_; j++) {
            y_ = j - 0.5;
            Ux_[i][j] = 4. * uMax_ / (L * L) * (L * y_ - y_ * y_);
            cout << i <<"   "<< j <<"   "<< Ux_[i][j] << endl;
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
    for (int k = 0; k < q_; k++) {
        f_[k] = new double* [lx_];
        for (int i = 0; i < lx_; i++) {
            f_[k][i] = new double[ly_];
            for (int j = 0; j < ly_; j++) {
                f_[k][i][j] = 0;
            }
        }
    }

    // Defining population array fStar(q,lx,ly)
    double*** fStar_;
    fStar_ = new double** [q_];
    for (int k = 0; k < q_; k++) {
        fStar_[k] = new double* [lx_];
        for (int i = 0; i < lx_; i++) {
            fStar_[k][i] = new double[ly_];
            for (int j = 0; j < ly_; j++) {
                fStar_[k][i][j] = 0;
            }
        }
    }

    // Defining population array fEq(q,lx,ly)
    double*** fEq_;
    fEq_ = new double** [q_];
    for (int k = 0; k < q_; k++) {
        fEq_[k] = new double* [lx_];
        for (int i = 0; i < lx_; i++) {
            fEq_[k][i] = new double[ly_];
            for (int j = 0; j < ly_; j++) {
                fEq_[k][i][j] = 0;
            }
        }
    }

    // Defining c * u
    double** cu_ = new double*[lx_];
    for (int i=0; i<lx_; i++){
        cu_[i] = new double[ly_];
        for (int j=0; j<ly_; j++){
            cu_[i][j] = 0;
        }
    }

    // Creating the mesh
    grid mesh_(lx_, ly_);

    // Defining the problem
    problem cylinder(mesh_, q_);
    
    // initialization of the population with macroscopic variables
    cylinder.initialize(uMax_, f_, rho_, Ux_, Uy_, w_, cx_, cy_);

    // Main loop
    for (int iter = 0; iter < iter_; iter++) {

        // Computing macroscopic variables with f_
        for (int i = 0; i < lx_; i++) {
            for (int j = 0; j < ly_; j++) {
                rho_[i][j] = f_[0][i][j] + f_[1][i][j] + f_[2][i][j] + 
                                f_[3][i][j] + f_[4][i][j] + f_[5][i][j] +
                                f_[6][i][j] + f_[7][i][j] + f_[8][i][j];

                Ux_[i][j] = ((f_[0][i][j] + f_[4][i][j] + f_[7][i][j]) -
                                (f_[2][i][j] + f_[5][i][j] + f_[6][i][j])) / rho_[i][j];

                Uy_[i][j] = ((f_[1][i][j] + f_[4][i][j] + f_[5][i][j]) -
                                (f_[3][i][j] + f_[6][i][j] + f_[7][i][j])) / rho_[i][j];
            }
        }

        // Computing cu and fEq in each step from previous macroscopic values
        for (int k = 0; k < q_; k++){
            for (int i = 0; i < lx_; i++){
                for (int j = 0; j < ly_; j++){
                    cu_[i][j] = 3.0*(cx_[k]*Ux_[i][j] + cy_[k]*Uy_[i][j]);
                    fEq_[k][i][j] = rho_[i][j] * w_[k] * (1.0 + cu_[i][j] + 
                    0.5 * pow(cu_[i][j],2.0) - 1.5 * (pow(Ux_[i][j], 2.0) + 
                    pow(Uy_[i][j], 2.0)));
                }
            }
        }

        // No BCs yet. Only collision and streaming for a disturbance in a large field.
        for (int k = 0; k < q_; k++) {
            for (int i = 0; i < lx_; i++) {
                for (int j = 0; j < ly_; j++) {
                    fStar_[k][i][j] = f_[k][i][j] - omega_*(f_[k][i][j] - fEq_[k][i][j]);
                }
            }
        }

        // Streaming for a disturbance in a large field.
        for (int k = 0; k < q_; k++) {
            for (int i = 1; i < lx_ - 1; i++) {
                for (int j = 1; j < ly_ - 1; j++) {
                    f_[k][i + cx_[k]][j + cy_[k]] = fStar_[k][i][j];
                }
            }
        }
    }  
        
    delete[] f_;
    delete[] fStar_;
    delete[] fEq_;
    delete[] rho_;
    delete[] Ux_;
    delete[] Uy_;
    delete[] cu_;

    return 0;
};