#include "grid.H"

grid::grid(const int lx__, const int ly__):
lx(lx__),
ly(ly__)
{
    int x[lx];
    int y[ly];

    for (int i = 0; i <= lx; i++) {
        x[i] = i;
    }

    for (int i = 0; i <= ly; i++) {
        y[i] = i;
    }

    cout << "This is grid constructor.\n";
}