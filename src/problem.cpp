#include "problem.H"
#include "grid.H"

problem::problem(const grid& mesh__, int q__):
lx(mesh__.lx),
ly(mesh__.ly),
q(q__)
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

int problem::initialize(const double rho__, double Ux__[100], const double Uy__)
{
    cout << "This is the initializer.\n";
    for(int i =0; i<100; i++){
        std::cout << Ux__[i] << endl;
    }
    return 0;
};