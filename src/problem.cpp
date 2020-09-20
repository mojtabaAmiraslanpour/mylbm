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

int problem::initialize(const double uMax__)
{
    cout << "This is the initializer.\n";

    int L = ly - 2;
    uMax = uMax__;
    // Initial conditions with macroscopic values
    double rho = 1.0;
    double Ux[100] = {};
    //double y_phys = ly_ - 0.5;
    for (int i = 0; i < ly; i++) { // Here we have not used y_phys which we have one
    // additional member in the vector beyond zero.
        Ux[i] = 4. * uMax / (L * L) * (L * i - i * i);
    }
    double Uy = 0.0;

    for(int i =0; i<100; i++){
        std::cout << Ux[i] << endl;
    }
    return 0;
};