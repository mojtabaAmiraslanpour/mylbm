#include <iostream>
#include "src/grid.h"
#include "src/problem.h"
#include <math.h>
#include <vector>

// GENERAL FLOW CONSTANTS
const int q_ = 9;
const int lx_ = 400; // number of cells in x-direction
const int ly_ = 100; // number of cells in y-direction
double pos_x_ = lx_/5.; // position of the cylinder; (exact
double pos_y_ = ly_/2. + 3; // y-symmetry is avoided)
double radius_ = ly_/10. + 1; // radius of the cylinder
double Re_ = 100; // Reynolds number
double uMax_ = 0.1; // maximum velocity of Poiseuille inflow
double nu_ = uMax_ * 2. * radius_ / Re_; // kinematic viscosity
double omega_ = 1. / (3. * nu_ + 0.5); // relaxation parameter, assuming deltaT = 1 so C_s = 1/3 p.117
int iter_ = 2000; // total number of iterations
int tPlot_ = 10; // cycles

// Using D2Q9
double w_[q_]  = {4./9., 1./9., 1./9., 1./9., 1./9.,
 1./36., 1./36., 1./36., 1./36.}; // Weights
int cx_[q_] = {0,   1,  0, -1,  0,    1,  -1,  -1,   1};
int cy_[q_] = {0,   0,  1,  0, -1,    1,   1,  -1,  -1};


int main(int, char**) {

    #include "createArrays.h"

    // Creating the mesh
    grid mesh_(lx_, ly_);

    // Defining the problem
    problem cylinder(mesh_, q_);

    double r;
    vector< pair <int, int> > obst;
    // Defining the cylinder nodes
    for (int i = 1; i < lx_ + 1; i++) {
        for (int j = 1; j < ly_ + 1; j++) {
            r = sqrt((i - pos_x_) * (i - pos_x_) + (j - pos_y_) * (j - pos_y_));
            if (r <= radius_) {
                obst.push_back(make_pair(i, j));
            }
        }
    }

    int m, n;
    int sizeObst = sizeof(obst);
    //cout << sizeObst << endl;
    //for (int i = 0; i < sizeObst; i++) {
    //    cout << obst[i].first << " " << obst[i].second << endl;
    //}
    
    // initialization of the population with macroscopic variables
    cylinder.initialize(uMax_, f_, rho_, Ux_, Uy_, w_, cx_, cy_);

    cylinder.writeVTK(rho_, Ux_, Uy_, 0);

    // Main loop
    for (int iter = 1; iter < iter_ + 1; iter++) {

        cout << endl;
        cout << "iter: " << iter << endl;

        // Computing macroscopic variables with f_
        for (int i = 1; i < lx_ + 1; i++) {
            for (int j = 1; j < ly_ + 1; j++) {
                rho_[i][j] = f_[0][i][j] + f_[1][i][j] + f_[2][i][j] + 
                                f_[3][i][j] + f_[4][i][j] + f_[5][i][j] +
                                f_[6][i][j] + f_[7][i][j] + f_[8][i][j];

                Ux_[i][j] = ((f_[1][i][j] + f_[5][i][j] + f_[8][i][j]) -
                                (f_[3][i][j] + f_[6][i][j] + f_[7][i][j])) / rho_[i][j];

                Uy_[i][j] = ((f_[2][i][j] + f_[5][i][j] + f_[6][i][j]) -
                                (f_[4][i][j] + f_[7][i][j] + f_[8][i][j])) / rho_[i][j];
            }
        }

        //Imposing BCs
        // f_[0][i + 0][j + 0] = fStar_[0][i][j];
        // f_[1][i + 1][j + 0] = fStar_[1][i][j];
        // f_[2][i + 0][j + 1] = fStar_[2][i][j];
        // f_[3][i - 1][j + 0] = fStar_[3][i][j];
        // f_[4][i + 0][j - 1] = fStar_[4][i][j];
        // f_[5][i + 1][j + 1] = fStar_[5][i][j];
        // f_[6][i - 1][j + 1] = fStar_[6][i][j];
        // f_[7][i - 1][j - 1] = fStar_[7][i][j];
        // f_[8][i + 1][j - 1] = fStar_[8][i][j];
        // In streaming we are sending populations to the neighbour nodes
        // so they can be used later in the collision step. So for the boundary nodes
        // you have to send populations to the neighbouring cells as normal. But there
        // is note that boundary nodes themselves wont have all the populations in the 
        // next collision step. Actually they are sending some populations beyond the
        // boundary. These populations will bounce back to themseleves which will be used
        // in the next collision step.

        // Setting Macroscopic BCs for inlet outlet
        for (int j = 2; j < ly_; j++) {
            // Setting U & rho for inlet fixed velocity
            y_ = j - 1.5;
            Ux_[1][j] = 4. * uMax_ / (L * L) * (L * y_ - y_ * y_);
            Uy_[1][j] = 0;
            rho_[1][j] = 1.0 / (1 - Ux_[1][j])
                * ((f_[3][1][j] + f_[6][1][j] + f_[7][1][j])
                    + (f_[0][1][j] + f_[2][1][j] + f_[3][1][j]
                        + f_[4][1][j] + f_[6][1][j] + f_[7][1][j]));

            // Setting U & rho for outlet fixed pressure
            rho_[lx_][j] = 1.0;
            Ux_[lx_][j] = -1.0 + 1.0 / rho_[lx_][j]
                * ((f_[0][lx_][j] + f_[2][lx_][j] + f_[4][lx_][j])
                    + 2 * (f_[1][lx_][j] + f_[5][lx_][j] + f_[8][lx_][j]));
            Uy_[lx_][j] = 0;
        }

        // Inlet vel BC
        for (int j = 2; j < ly_; j++) {
            f_[1][1][j] = f_[3][1][j] + (2.0 / 3.0) * rho_[1][j] * Ux_[1][j]; // inlet vel. bounce back
            f_[8][1][j] = f_[6][1][j] + 0.5 * (f_[2][1][j] - f_[4][1][j])
                + (1.0 / 6.0) * rho_[1][j] * Ux_[1][j] - 0.5 * rho_[1][j] * Uy_[1][j]; // inlet vel. bounce back
            f_[5][1][j] = f_[7][1][j] + 0.5 * (f_[4][1][j] - f_[2][1][j])
                + (1.0 / 6.0) * rho_[1][j] * Ux_[1][j] + 0.5 * rho_[1][j] * Uy_[1][j]; // inlet vel. bounce back
        }

        // Outlet zeroGrad vel
        for (int j = 2; j < ly_; j++) {
            f_[3][lx_][j] = f_[1][lx_][j] - 2.0 / 3.0 * rho_[lx_][j] * Ux_[lx_][j]; // Zou/He Pressure Boundary
            f_[7][lx_][j] = f_[5][lx_][j] + 0.5 * (f_[2][lx_][j] - f_[4][lx_][j])
                - (1.0 / 6.0) * rho_[lx_][j] * Ux_[lx_][j] - 0.5 * rho_[lx_][j] * Uy_[lx_][j]; // Zou/He Pressure Boundary
            f_[6][lx_][j] = f_[8][lx_][j] + 0.5 * (f_[4][lx_][j] - f_[2][lx_][j])
                - (1.0 / 6.0) * rho_[lx_][j] * Ux_[lx_][j] + 0.5 * rho_[lx_][j] * Uy_[lx_][j]; // Zou/He Pressure Boundary
        }

        // Collision
        for (int k = 0; k < q_; k++) {
            for (int i = 1; i < lx_ + 1; i++) {
                for (int j = 1; j < ly_ + 1; j++) {

                    cu_[i][j] = 3.0*(cx_[k]*Ux_[i][j] + cy_[k]*Uy_[i][j]);

                    fEq_[k][i][j] = rho_[i][j] * w_[k] * (1.0 + cu_[i][j] + 
                    0.5 * pow(cu_[i][j],2.0) - 1.5 * (pow(Ux_[i][j], 2.0) + 
                    pow(Uy_[i][j], 2.0)));

                    fStar_[k][i][j] = f_[k][i][j] - omega_*(f_[k][i][j] - fEq_[k][i][j]);
                }
            }
        }

        // Lower wall boundary bounce back streaming
        for (int i = 2; i < lx_; i++) {
            fStar_[2][i][1] = f_[4][i][1]; // bounce back
            fStar_[5][i][1] = f_[7][i][1]; // bounce back
            fStar_[6][i][1] = f_[8][i][1]; // bounce back
        }

        // Upper wall boundary bounce back streaming
        for (int i = 2; i < lx_; i++) {
            fStar_[4][i][ly_] = f_[2][i][ly_]; // bounce back
            fStar_[7][i][ly_] = f_[5][i][ly_]; // bounce back
            fStar_[8][i][ly_] = f_[6][i][ly_]; // bounce back
        }

        // Inlet-Bottom corner
        fStar_[6][1][1] = -0.5 * (f_[0][1][1] + 2 * (f_[1][1][1] + f_[2][1][1] + f_[5][1][1]));
        fStar_[8][1][1] = f_[6][1][1];
        fStar_[1][1][1] = f_[3][1][1];
        fStar_[2][1][1] = f_[4][1][1];
        fStar_[5][1][1] = f_[7][1][1];

        // Inlet-Top corner
        fStar_[5][1][ly_] = -0.5 * (f_[0][1][ly_] + 2 * (f_[1][1][ly_] + f_[4][1][ly_] + f_[8][1][ly_]));
        fStar_[7][1][ly_] = f_[5][1][ly_];
        fStar_[1][1][ly_] = f_[3][1][ly_];
        fStar_[4][1][ly_] = f_[2][1][ly_];
        fStar_[8][1][ly_] = f_[6][1][ly_];

        // Outlet-Bottom corner
        fStar_[5][1][ly_] = -0.5 * (f_[0][1][ly_] + 2 * (f_[2][1][ly_] + f_[3][1][ly_] + f_[6][1][ly_]));
        fStar_[7][1][ly_] = f_[5][1][ly_];
        fStar_[3][1][ly_] = f_[1][1][ly_];
        fStar_[2][1][ly_] = f_[4][1][ly_];
        fStar_[6][1][ly_] = f_[8][1][ly_];

        // // Outlet-Top corner
        fStar_[6][1][ly_] = -0.5 * (f_[0][1][ly_] + 2 * (f_[3][1][ly_] + f_[4][1][ly_] + f_[7][1][ly_]));
        fStar_[8][1][ly_] = f_[6][1][ly_];
        fStar_[3][1][ly_] = f_[1][1][ly_];
        fStar_[4][1][ly_] = f_[2][1][ly_];
        fStar_[7][1][ly_] = f_[5][1][ly_];

        // Cylinder bounce back
        for (int i = 0; i < sizeObst; i++) {
            m = obst[i].first;
            n = obst[i].second;
            fStar_[1][m][n] = f_[3][m][n];
            fStar_[2][m][n] = f_[4][m][n];
            fStar_[3][m][n] = f_[1][m][n];
            fStar_[4][m][n] = f_[2][m][n];
            fStar_[5][m][n] = f_[7][m][n];
            fStar_[6][m][n] = f_[8][m][n];
            fStar_[7][m][n] = f_[5][m][n];
            fStar_[8][m][n] = f_[6][m][n];
        }

        // // Check values
        // //for (int k = 0; k < q_; k++) {
        //     //cout << "k = " << k << endl;
        //     for (int j = 2; j < ly_; j++) {
        //         //for (int i = 0; i < lx_; i++) {
        //             cout << fStar_[4][1][j] << " ";
        //         //}
        //         cout << "\n";
        //     }
        // //}

        // Streaming of the fluid nodes
        for (int k = 0; k < q_; k++) {
            for (int i = 1; i < lx_ + 1; i++) {
                for (int j = 1; j < ly_ + 1; j++) {
                    f_[k][i + cx_[k]][j + cy_[k]] = fStar_[k][i][j];
                }
            }
        }

        if (iter % tPlot_ == 0) {
            cylinder.writeVTK(rho_, Ux_, Uy_, iter);
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