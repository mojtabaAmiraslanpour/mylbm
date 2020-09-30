#ifndef problem_H
#define problem_H

#include <iostream>
#include "grid.h"
#include <fstream>
#include <math.h>

using namespace std;

class problem
{
private:
    double uMax;

public:
    const int q;
    const int lx;
    const int ly;

    problem(const grid& mesh__, int q__);
    ~problem();

    int initialize(const double uMax__, double*** f__, double** rho__, double** Ux__, double** Uy__,
const double* w__, const int* cx__, const int* cy__);

    int writeVTK(double** rho__, double** Ux__, double** Uy__, int iter__);

};

#endif