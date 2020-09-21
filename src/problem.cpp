#include "problem.h"
#include "grid.h"
#include <math.h>

// Constructor
problem::problem(const grid& mesh__, int q__,const int *cx__,const int *cy__, const double *w__):
lx(mesh__.lx),
ly(mesh__.ly),
q(q__),
cx(cx__),
cy(cy__),
w(w__)
{
    cout << "This is problem constructor.\n";

};

// Destructor
problem::~problem() {}

int problem::initialize(const double uMax__, double*** f__, double** rho__, double** Ux__, double** Uy__)
{
    cout << "This is the initializer.\n";

    // Defining c * u
    double** cu = new double*[lx];
    for (int i=0; i<lx; i++){
        cu[i] = new double[ly];
        for (int j=0; j<ly; j++){
            cu[i][j] = 0;
        }
    }

    // Initializing f (populating cu and f)
    for (int k=0; k<q; k++){
        for (int i=0; i<lx; i++){
            for (int j=0; j<ly; j++){
                cu[i][j] = 3.0*(cx[k]*Ux__[i][j] + cy[k]*Uy__[i][j]);
                //f__[k][i][j] = rho__[i][j] * w[k] * (1.0 + cu[i][j] + 
                //0.5 * pow(cu[i][j],2.0) - 1.5 * (pow(Ux__[i][j], 2.0) + 
                //pow(Uy__[i][j], 2.0)));
                f__[k][i][j] = 1.0;
            }
        }
    }

    return 0;
};