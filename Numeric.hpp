#include "Matrix.hpp"
#include "SparseMatrix.hpp"
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>


std::vector<std::complex<double>> eig(Matrix m);

std::vector<double> solve(Matrix A, std::vector<double> b);

std::vector<double> solve(const SparseMatrix& A, std::vector<double> b);

std::vector<double> solvePetsc(const SparseMatrix& M, std::vector<double> b);