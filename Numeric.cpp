#include <iostream>
#include <lapacke.h>
#include "suitesparse/umfpack.h"
#include "Numeric.hpp"

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


std::vector<double> solve(const SparseMatrix& A, std::vector<double> b)
{
    std::vector<double> out(b.size());

    const std::vector<Triplet>& triplets = A.getTriplets();

    int n = A.getN();
    int m = A.getM();
    int nnz = A.getNnz();

    std::vector<int> ii(nnz);
    std::vector<int> jj(nnz);
    std::vector<double> AA(nnz);

    for (int i = 0; i < triplets.size(); i++)
    {
        ii[i] = triplets[i].i;
        jj[i] = triplets[i].j;
        AA[i] = triplets[i].val;
    }

    std::vector<int> cptr = std::vector<int>(m + 1);
    std::vector<int> cidx = std::vector<int>(nnz);
    std::vector<double> val = std::vector<double>(nnz);

    umfpack_di_triplet_to_col(n, m, nnz, ii.data(), jj.data(), AA.data(), cptr.data(), cidx.data(), val.data(), NULL);

    void *Symbolic, *Numeric;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];

    umfpack_di_defaults(Control);

    Control[UMFPACK_PRL] = 5;
    //umfpack_di_report_control(Control);

    umfpack_di_symbolic(n, m, cptr.data(), cidx.data(), val.data(), &Symbolic, Control, Info);

    umfpack_di_numeric(cptr.data(), cidx.data(), val.data(), Symbolic, &Numeric, Control, Info);
    if (Info[UMFPACK_STATUS] != UMFPACK_OK)
    {
        std::cout << "Chyba pri numericke faktorizaci.\n";
        return std::vector<double>();
    }

    umfpack_di_free_symbolic(&Symbolic);

    umfpack_di_solve(UMFPACK_A, cptr.data(), cidx.data(), val.data(), out.data(), b.data(), Numeric, Control, Info);

    umfpack_di_free_numeric(&Numeric);

    //umfpack_di_report_info(Control, Info);

    return out;
}




std::vector<double> solvePetsc(const SparseMatrix& M, std::vector<double> b)
{
    const std::vector<Triplet>& triplets = M.getTriplets();

    int n = M.getM();

    Vec y;
    VecCreateSeqWithArray(PETSC_COMM_WORLD, PETSC_DECIDE, n, b.data(), &y);
    
    Mat A;    
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
    MatSetFromOptions(A);
    MatSetUp(A);

    PetscInt istart, iend;
    MatGetOwnershipRange(A, &istart, &iend);
    for (int i = 0; i < triplets.size(); i++)
    {
        MatSetValue(A, triplets[i].i, triplets[i].j, triplets[i].val, INSERT_VALUES);
    }
    
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    Vec x;  
    VecDuplicate(y, &x);

    KSP solver;
    KSPCreate(PETSC_COMM_WORLD, &solver);
    KSPSetOperators(solver, A, A);
  
    KSPSetType(solver, KSPCG);
    PC prec;
    KSPGetPC(solver, &prec);
    PCSetType(prec, PCJACOBI);

    KSPSetFromOptions(solver);
    KSPSetUp(solver);

    KSPSolve(solver, y, x);

    //VecView(x, PETSC_VIEWER_STDOUT_WORLD);

    double *outPtr;
    VecGetArray(x, &outPtr);
    std::vector<double> out(outPtr, outPtr + n);

    KSPDestroy(&solver);
    MatDestroy(&A);
    VecDestroy(&x);
    VecDestroy(&y);
    
    return out;
}