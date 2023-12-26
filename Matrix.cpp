#include <iostream>
#include <algorithm>

#include <lapacke.h>

#include "Matrix.hpp"


std::vector<std::complex<double>> eig(Matrix m)
{
    std::vector<double> auxReal(m.rows());
    std::vector<double> auxImag(m.rows());
    std::vector<double> auxVl(m.rows());
    std::vector<double> auxVr(m.rows());
    
    int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'N', m.rows(), m.data(), m.lda(), auxReal.data(), auxImag.data(), auxVl.data(), m.rows(), auxVr.data(), m.rows());

    std::vector<std::complex<double>> out(auxReal.size());
    for (int i = 0; i < auxReal.size(); i++)
    {
        out[i] = std::complex<double>(auxReal[i], auxImag[i]);
    }
    
    return out;
}

std::vector<double> solve(Matrix A, std::vector<double> b)
{
    std::vector<int> ipiv(A.rows());

    int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, A.rows(), 1, A.data(), A.lda(), ipiv.data(), b.data(), 1);

    return b;
}