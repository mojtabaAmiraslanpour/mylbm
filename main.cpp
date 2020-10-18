#include <iostream>
#include "src/grid.h"
#include "src/problem.h"
#include <math.h>

// GENERAL FLOW CONSTANTS
const int q_ = 9;
const int lx_ = 3; // number of cells in x-direction
const int ly_ = 6; // number of cells in y-direction
double pos_x_ = lx_/5; // position of the cylinder; (exact
double pos_y_ = ly_/2+3; // y-symmetry is avoided)
double radius_ = ly_/10+1; // radius of the cylinder
double Re_ = 100; // Reynolds number
double uMax_ = 0.1; // maximum velocity of Poiseuille inflow
double nu_ = uMax_ * 2. * radius_ / Re_; // kinematic viscosity
double omega_ = 1. / (3 * nu_+1./2.); // relaxation parameter, assuming deltaT = 1 so C_s = 1/3 p.117
int iter_ = 2; // total number of iterations
int tPlot_ = 10; // cycles

// Using D2Q9
double w_[q_]  = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0}; // Weights
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
            //cout << i <<"   "<< j <<"   "<< Ux_[i][j] << endl;
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

        cout << endl;
        cout << "iter: " << iter << endl;

        // Computing macroscopic variables with f_
        for (int i = 0; i < lx_; i++) {
            for (int j = 0; j < ly_; j++) {
                rho_[i][j] = f_[0][i][j] + f_[1][i][j] + f_[2][i][j] + 
                                f_[3][i][j] + f_[4][i][j] + f_[5][i][j] +
                                f_[6][i][j] + f_[7][i][j] + f_[8][i][j];

                Ux_[i][j] = ((f_[1][i][j] + f_[5][i][j] + f_[8][i][j]) -
                                (f_[3][i][j] + f_[6][i][j] + f_[7][i][j])) / rho_[i][j];

                Uy_[i][j] = ((f_[2][i][j] + f_[5][i][j] + f_[6][i][j]) -
                                (f_[4][i][j] + f_[7][i][j] + f_[8][i][j])) / rho_[i][j];
            }
        }


        //BCs
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

        // // Lower wall boundary bounce back streaming
        // for (int i = 1; i < lx_ - 1; i++) {
        //     f_[0][i][0] = fStar_[0][i][0]; // form the node itself
        //     f_[1][i + 1][0] = fStar_[1][i][0]; // to neighbour boundary node
        //     f_[2][i][1] = fStar_[2][i][0]; // Upper fluid node
        //     f_[3][i - 1][0] = fStar_[3][i][0]; // to neighbour boundary node
        //     f_[2][i][0] = fStar_[4][i][0]; // bounce back
        //     f_[5][i + 1][1] = fStar_[5][i][0]; // Upper fluid node
        //     f_[6][i - 1][1] = fStar_[6][i][0]; // Upper fluid node
        //     f_[5][i][0] = fStar_[7][i][0]; // bounce back
        //     f_[6][i][0] = fStar_[8][i][0]; // bounce back
        // }
        
        // // Upper wall boundary bounce back streaming
        // for (int i = 1; i < lx_ - 1; i++) {
        //     f_[0][i][ly_ - 1] = fStar_[0][i][ly_ - 1]; // form the node itself
        //     f_[1][i + 1][ly_ - 1] = fStar_[1][i][ly_ - 1]; // to neighbour boundary node
        //     f_[4][i][ly_ - 1] = fStar_[2][i][ly_ - 1]; // bounce back
        //     f_[3][i - 1][ly_ - 1] = fStar_[3][i][ly_ - 1]; // to neighbour boundary node
        //     f_[4][i][ly_ - 2] = fStar_[4][i][ly_ - 1]; // Lower fluid node
        //     f_[7][i][ly_ - 1] = fStar_[5][i][ly_ - 1]; // bounce back
        //     f_[8][i][ly_ - 1] = fStar_[6][i][ly_ - 1]; // bounce back
        //     f_[7][i - 1][ly_ - 2] = fStar_[7][i][ly_ - 1]; // Lower fluid node
        //     f_[8][i + 1][ly_ - 2] = fStar_[8][i][ly_ - 1]; // Lower fluid node
        // }

        
        for (int j = 1; j < ly_ - 1; j++) {
            // Setting U & rho for inlet
            y_ = j - 0.5;
            Ux_[0][j] = 4. * uMax_ / (L * L) * (L * y_ - y_ * y_);
            Uy_[0][j] = 0;
            rho_[0][j] = 1.0 / (1 - Ux_[0][j])
                       * ((f_[3][0][j] + f_[6][0][j] + f_[7][0][j])
                       +  (f_[0][0][j] + f_[2][0][j] + f_[3][0][j]
                       +   f_[4][0][j] + f_[6][0][j] + f_[7][0][j]));
            
            rho_[lx_ - 1][j] = 1.0;
            Ux_[lx_ - 1][j] = -1.0 + 1.0 / rho_[lx_ - 1][j]
                            * ((f_[0][lx_ - 1][j] + f_[2][lx_ - 1][j] + f_[4][lx_ - 1][j])
                            + 2 * (f_[1][lx_ - 1][j] + f_[5][lx_ - 1][j] + f_[8][lx_ - 1][j]));
            Uy_[lx_ - 1][j] = 0;
        }



        // Inlet vel BC
        for (int j = 1; j < ly_ - 1; j++) {
            f_[1][0][j] = f_[3][0][j] + (2.0/3.0) * rho_[0][j] * Ux_[0][j]; // inlet vel. bounce back
            f_[8][0][j] = f_[6][0][j] + 0.5 * (f_[2][0][j] - f_[4][0][j])
                        + (1.0/6.0) * rho_[0][j] * Ux_[0][j] - 0.5 * rho_[0][j] * Uy_[0][j]; // inlet vel. bounce back
            f_[5][0][j] = f_[7][0][j] + 0.5 * (f_[4][0][j] - f_[2][0][j])
                        + (1.0/6.0) * rho_[0][j] * Ux_[0][j] + 0.5 * rho_[0][j] * Uy_[0][j]; // inlet vel. bounce back
        }

        // Outlet zeroGrad vel
        for (int j = 1; j < ly_ - 1; j++) {
            f_[3][lx_ - 1][j] = f_[1][lx_ - 1][j] - 2.0/3.0 * rho_[lx_ - 1][j] * Ux_[lx_ - 1][j]; // Zou/He Pressure Boundary
            f_[7][lx_ - 1][j] = f_[5][lx_ - 1][j] + 0.5 * (f_[2][lx_ - 1][j] - f_[4][lx_ - 1][j])
                        - (1.0/6.0) * rho_[lx_ - 1][j] * Ux_[lx_ - 1][j] - 0.5 * rho_[lx_ - 1][j] * Uy_[lx_ - 1][j]; // Zou/He Pressure Boundary
            f_[6][lx_ - 1][j] = f_[8][lx_ - 1][j] + 0.5 * (f_[4][lx_ - 1][j] - f_[2][lx_ - 1][j])
                        - (1.0/6.0) * rho_[lx_ - 1][j] * Ux_[lx_ - 1][j] + 0.5 * rho_[lx_ - 1][j] * Uy_[lx_ - 1][j]; // Zou/He Pressure Boundary
        }

        // // Inlet-Bottom corner
        // f_[0][0][0] = fStar_[0][0][0]; // form the node itself
        // f_[1][0][0] = fStar_[3][0][0] - 2 * w_[3] * rho_[1][0] * (cx_[3] * Ux_[0][0] + cy_[3] * Uy_[0][0]); // inlet vel. bounce back
        // f_[2][0][0] = fStar_[4][0][0]; // Bounce back from lower wall
        // f_[8][0][0] = fStar_[6][0][0] - 2 * w_[6] * rho_[1][0] * (cx_[6] * Ux_[0][0] + cy_[6] * Uy_[0][0]); // inlet vel. bounce back
        // f_[5][0][0] = fStar_[7][0][0]; // Bounce back from lower wall
        // f_[6][0][0] = fStar_[8][0][0]; // Bounce back from lower wall

        // // Inlet-Top corner
        // f_[0][0][ly_ - 1] = fStar_[0][0][ly_ - 1]; // form the node itself
        // f_[1][0][ly_ - 1] = fStar_[3][0][ly_ - 1] - 2 * w_[3] * rho_[1][ly_ - 1] * (cx_[3] * Ux_[0][ly_ - 1] + cy_[3] * Uy_[0][ly_ - 1]); // inlet vel. bounce back
        // f_[4][0][ly_ - 1] = fStar_[2][0][ly_ - 1]; // Bounce back from upper wall
        // f_[5][0][ly_ - 1] = fStar_[7][0][ly_ - 1] - 2 * w_[7] * rho_[1][ly_ - 1] * (cx_[7] * Ux_[0][ly_ - 1] + cy_[7] * Uy_[0][ly_ - 1]); // inlet vel. bounce back
        // f_[7][0][ly_ - 1] = fStar_[5][0][ly_ - 1]; // Bounce back from upper wall
        // f_[8][0][ly_ - 1] = fStar_[6][0][ly_ - 1]; // Bounce back from upper wall

        // // Outlet-Bottom corner
        // f_[0][lx_ - 1][0] = fStar_[0][lx_ - 1][0]; // form the node itself
        // f_[2][lx_ - 1][0] = fStar_[4][lx_ - 1][0]; // Bounce back from lower wall
        // f_[3][lx_ - 1][0] = fStar_[3][lx_ - 2][0]; // 1st order extrapolation from the inner fluid node, A A Mohammad P.119
        // f_[5][lx_ - 1][0] = fStar_[7][lx_ - 1][0]; // Bounce back from lower wall
        // f_[6][lx_ - 1][0] = fStar_[8][lx_ - 1][0]; // Bounce back from lower wall
        // f_[7][lx_ - 1][0] = fStar_[7][lx_ - 2][0]; // 1st order extrapolation from the inner fluid node, A A Mohammad P.119

        // // Outlet-Top corner
        // f_[0][lx_ - 1][ly_ - 1] = fStar_[0][lx_ - 1][ly_ - 1]; // form the node itself
        // f_[3][lx_ - 1][ly_ - 1] = fStar_[3][lx_ - 2][ly_ - 1]; // 1st order extrapolation from the inner fluid node, A A Mohammad P.119
        // f_[4][lx_ - 1][ly_ - 1] = fStar_[2][lx_ - 1][ly_ - 1]; // Bounce back from upper wall
        // f_[6][lx_ - 1][ly_ - 1] = fStar_[6][lx_ - 2][ly_ - 1]; // 1st order extrapolation from the inner fluid node, A A Mohammad P.119
        // f_[7][lx_ - 1][ly_ - 1] = fStar_[5][lx_ - 1][ly_ - 1]; // Bounce back from upper wall
        // f_[8][lx_ - 1][ly_ - 1] = fStar_[6][lx_ - 1][ly_ - 1]; // Bounce back from upper wall
        
        // Computing cu and fEq in each step & Collision
        for (int k = 0; k < q_; k++) {
            for (int i = 0; i < lx_; i++) {
                for (int j = 0; j < ly_; j++) {
                    cu_[i][j] = 3.0*(cx_[k]*Ux_[i][j] + cy_[k]*Uy_[i][j]);
                    fEq_[k][i][j] = rho_[i][j] * w_[k] * (1.0 + cu_[i][j] + 
                    0.5 * pow(cu_[i][j],2.0) - 1.5 * (pow(Ux_[i][j], 2.0) + 
                    pow(Uy_[i][j], 2.0)));
                    fStar_[k][i][j] = f_[k][i][j] - omega_*(f_[k][i][j] - fEq_[k][i][j]);
                }
            }
        }

            // Check values
            //for (int k = 0; k < q_; k++) {
                //cout << "k = " << k << endl;
                for (int j = 1; j < ly_ - 1; j++) {
                    //for (int i = 0; i < lx_; i++) {
                        cout << fStar_[4][lx_ - 3][j] << " ";
                    //}
                    cout << "\n";
                }
            //}

        // Streaming of the fluid nodes
        for (int k = 0; k < q_; k++) {
            for (int i = 1; i < lx_ - 1; i++) {
                for (int j = 1; j < ly_ - 1; j++) {
                    f_[k][i + cx_[k]][j + cy_[k]] = fStar_[k][i][j];
                }
            }
        }

        // Inlet boundary streaming
        for (int j = 1; j < ly_ - 1; j++) {
            f_[0][0][j] = f_[0][0][j]; // form the node itself
            f_[1][1][j] = f_[1][0][j]; // right fluid node
            f_[2][0][j + 1] = f_[2][0][j]; // to neighbour boundary node
            f_[4][0][j - 1] = f_[4][0][j]; // to neighbour boundary node
            f_[5][1][j + 1] = f_[5][0][j]; // right fluid node
            f_[8][1][j - 1] = f_[8][0][j]; // right fluid node
        }

        // Outlet boundary streaming
        for (int j = 1; j < ly_ - 1; j++) {
            f_[0][lx_ - 1][j] = f_[0][lx_ - 1][j]; // form the node itself
            f_[2][lx_ - 1][j + 1] = f_[2][lx_ - 1][j]; // to neighbour boundary node
            f_[3][lx_ - 2][j] = f_[3][lx_ - 1][j]; // left fluid node
            f_[4][lx_ - 1][j - 1] = f_[4][lx_ - 1][j]; // to neighbour boundary node
            f_[6][lx_ - 2][j + 1] = f_[6][lx_ - 1][j]; // left fluid node
            f_[7][lx_ - 2][j - 1] = f_[7][lx_ - 1][j]; // left fluid node
        }

        // if (iter % tPlot_ == 0) {
        //     cylinder.writeVTK(rho_, Ux_, Uy_, iter);
        // }

        // See if populated
        // if (iter = iter_) {
        //     int p;
        //     for (int k = 0; k < q_; k++) {
        //         cout << "k = " << k << endl;

        //         for (int j = 0; j < ly_; j++) {

        //             for (int i = 0; i < lx_; i++) {

        //                 if (f_[k][i][j] != 0) {
        //                     p = 1;
        //                 }
        //                 else
        //                 {
        //                     p = 0;
        //                 }

        //                 cout << p << "," << f_[k][i][j] << " ";
        //             }

        //             cout << "\n";
        //         }
        //         cout << "\n";
        //     }
        // } 
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