#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP

#include <complex>
#include <vector>
#include <memory>

#include "Matrix.hpp"
#include "Triplet.hpp"

class SparseMatrix
{
    public:
        SparseMatrix() = delete;
        SparseMatrix(int m, int n) : m_(m), n_(n), data_(std::vector<Triplet>()) {}
        SparseMatrix(const Matrix& M);

        virtual ~SparseMatrix() {}

        double operator() (int i, int j) const;
        double& operator()(int i, int j);

        void sort();

        void addNew(Triplet t);

        int getM() const;
        int getN() const;
        int getNnz() const;

        const std::vector<Triplet>& getTriplets() const;

    protected:
        int m_;
        int n_;
        std::vector<Triplet> data_;
};

#endif // SPARSEMATRIX_HPP