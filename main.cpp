#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <chrono>
#include <fenv.h>
#include <mpi.h>

#include "Matrix.hpp"
#include "Numeric.hpp"

template <typename T>
void saveResults(std::string resultName, std::vector<T> data)
{
    std::ofstream f;
	f.open("results/" + resultName + ".txt", std::ios::out);

    for (int i = 0; i < data.size(); i++)
    {
        f << data[i] << "\n";
    }

    f << std::endl;
    f.close();    
}

int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    CHKERRQ(PetscInitialize(&argc, &argv, (char*)0, 0));

    int n = 1000;
    int m = 32; // +-sqrt(1000)

    Matrix L(n, n);

    for (int i = 0; i < n; i++)
    {
        L[i][i] = 4.0;

        if (i+1 < n)  { L[i][i+1] = -1.0; }

        if (i-1 >= 0) { L[i][i-1] = -1.0; }
        
        if (i+m < n)  { L[i][i+m] = -1.0; }

        if (i-m >= 0) { L[i][i-m] = -1.0; }
    }

    SparseMatrix S(L);
    
    std::vector<double> b(n, 2.0);

    std::vector<std::complex<double>> resultEigen = eig(L);
    
    std::vector<double> resultLapack = solve(L, b);

    std::vector<double> resultUmfpack = solve(S, b);

    std::vector<double> resultPetsc = solvePetsc(S, b);

    saveResults("vlCisla", resultEigen);
    saveResults("Lapack", resultLapack);
    saveResults("Umfpack", resultUmfpack);
    saveResults("Petsc", resultPetsc);


    CHKERRQ(PetscFinalize());

    return 0;
}