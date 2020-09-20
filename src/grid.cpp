#include "grid.h"

grid::grid(const int lx__, const int ly__):
lx(lx__),
ly(ly__)
{
    int* x = new int[lx];
    int* y = new int[ly];

    for (int i = 0; i < lx; i++) {
        x[i] = i;
    }

    for (int i = 0; i < ly; i++) {
        y[i] = i;
    }

    cout << "This is grid constructor.\n";

    delete[] x;
    delete[] y;
}