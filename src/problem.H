#ifndef problem_H
#define problem_H

#include <iostream>
#include "grid.H"

using namespace std;

class problem
{
    public:

    problem(const grid& mesh__, double f__);

    int initialize();
};

#endif