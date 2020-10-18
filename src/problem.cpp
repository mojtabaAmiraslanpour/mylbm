#include "problem.h"
#include <iostream>
#include <sstream>
using namespace std;


// Constructor
problem::problem(const grid& mesh__, int q__):
lx(mesh__.lx),
ly(mesh__.ly),
q(q__)
{
    //cout << "This is problem constructor.\n";
}

// Destructor
problem::~problem() {}

int problem::initialize(const double uMax__, double*** f__, double** rho__, double** Ux__, double** Uy__,
const double* w__, const int* cx__, const int* cy__)
{
    //cout << "This is the initializer.\n";

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
                cu[i][j] = 3.0*(cx__[k]*Ux__[i][j] + cy__[k]*Uy__[i][j]);
                f__[k][i][j] = rho__[i][j] * w__[k] * (1.0 + cu[i][j] + 
                0.5 * pow(cu[i][j],2.0) - 1.5 * (pow(Ux__[i][j], 2.0) + 
                pow(Uy__[i][j], 2.0)));

                //cout << f__[k][i][j] << endl;
            }
        }
    }

    return 0;
}

int problem::writeVTK(double** rho__, double** Ux__, double** Uy__, int iter__) {

    string iter_str;
    stringstream transfer;
    transfer << iter__;
    transfer >> iter_str;

    ofstream UxResults;
    UxResults.open("file_no_" + iter_str + ".vtk");
    UxResults
        << "# vtk DataFile Version 3.0\n"
        << "first dataset\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_POINTS\n";
    UxResults << "DIMENSIONS " << lx << " " << ly << " 1\n";
    UxResults << "ORIGIN 0 0 0\n"
        << "SPACING 1 1 0\n";
    UxResults << "POINT_DATA " << lx * ly << endl;
    UxResults << "SCALARS U(x) float\n"
        << "LOOKUP_TABLE default\n";

    for (int j = 0; j < ly; j++) {
        for (int i = 0; i < lx; i++) {
            UxResults << Ux__[i][j] << " ";
        }
        UxResults << endl;
    }

    UxResults.close();

    return 0;
}