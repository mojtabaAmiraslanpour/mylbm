#ifndef problem_H
#define problem_H

#include <iostream>
#include "grid.H"

using namespace std;

class problem
{
    public:

    const int q;
    const int lx;
    const int ly;

    // Declaration of the population function
    double*** f;

    problem(const grid& mesh__, int q__);
    ~problem();

    int initialize(const double rho__, double Ux__[100], const double Uy__);
};

#endif