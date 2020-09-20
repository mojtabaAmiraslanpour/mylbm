#include "problem.h"
#include "grid.h"
#include <math.h>

problem::problem(const grid& mesh__, int q__,const int *cx__,const int *cy__, const double *w__):
lx(mesh__.lx),
ly(mesh__.ly),
q(q__),
cx(cx__),
cy(cy__),
w(w__)
{
    cout << "This is problem constructor.\n";

    // Defining population array f(q,lx,ly)
    f = new double**[q];
    for(int i =0; i<q; i++){
        f[i] = new double*[lx];
        for(int j =0; j<lx; j++){
            f[i][j] = new double[ly];
            for(int k = 0; k<ly;k++){
                f[i][j][k] = 0;
            }
        }
    }
};

problem::~problem() { 
    /*
    for(int i = 0; i < q; i++){
        delete[] f[i];  
        for(int j = 0; j < lx; j++){
            delete[] f[i][j];
        }
    }*/
    delete[] f; 
    }

int problem::initialize(const double uMax__)
{
    cout << "This is the initializer.\n";

    int L = ly - 2;
    uMax = uMax__;

    // Initial conditions with macroscopic values -----------------------------

    // Defining rho
    double** rho = new double*[lx];
    for (int i=0; i<lx; i++){
        rho[i] = new double[ly];
        for (int j=0; j<ly; j++){
            rho[i][j] = 1.0;
        }
    }

    // Defining Ux
    //double y_phys = ly_ - 0.5;
    double** Ux = new double*[lx];
    for (int i=0; i<lx; i++){
        Ux[i] = new double[ly];
        for (int j=0; j<ly; j++){
            Ux[i][j] = 4. * uMax / (L * L) * (L * j - j * i); // Not used y_phys
        }
    }

    // Defining Uy
    double** Uy = new double*[lx];
    for (int i=0; i<lx; i++){
        Uy[i] = new double[ly];
        for (int j=0; j<ly; j++){
            Uy[i][j] = 0;
        }
    }

    // Defining c * u
    double** cu = new double*[lx];
    for (int i=0; i<lx; i++){
        cu[i] = new double[ly];
        for (int j=0; j<ly; j++){
            cu[i][j] = 0;
        }
    }

    // Initializing f
    for (int k=0; k<q; k++){
        for (int i=0; i<lx; i++){
            for (int j=0; j<ly; j++){
                cu[i][j] = 3.0*(cx[k]*Ux[i][j] + cy[k]*Uy[i][j]);
                f[k][i][j] = rho[i][j] * w[k] * (1.0 + cu[i][j] + 
                0.5 * pow(cu[i][j],2.0) - 1.5 * (pow(Ux[i][j], 2.0) + 
                pow(Uy[i][j], 2.0)));
            }
        }
    }

    // for(int i =0; i<100; i++){
    //     std::cout << Ux[i] << endl;
    // }

    return 0;
};