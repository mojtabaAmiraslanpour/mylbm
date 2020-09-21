#ifndef problem_H
#define problem_H

#include <iostream>
#include "grid.h"

using namespace std;

class problem
{
private:
    double uMax;

public:
    const int q;
    const int lx;
    const int ly;
    const int *cx;
    const int *cy;
    const double *w;

    problem(const grid& mesh__, int q__,const int *cx__,const int *cy__, const double *w__);
    ~problem();

    int initialize(const double uMax__, double*** f__, double** rho__, double** Ux__, double** Uy__);

};

#endif