#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <chrono>
#include <fenv.h>


#include "Matrix.hpp"

int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    int n = 1000;
    int m = 32; // +-sqrt(1000)

    Matrix L(n, n);

    for (int i = 0; i < n; i++)
    {
        L[i][i] = 4.0;

        if (i+1 < n)
        {
            L[i][i+1] = -1.0;
        }

        if (i-1 >= 0)
        {
            L[i][i-1] = -1.0;
        }
        
        if (i+m < n)
        {
            L[i][i+m] = -1.0;
        }

        if (i-m >= 0)
        {
            L[i][i-m] = -1.0;
        }
    }
    
    std::vector<double> b(n, 2.0);

    std::vector<std::complex<double>> resultEigen = eig(L);
    std::vector<double> result = solve(L, b);

    return 0;
}