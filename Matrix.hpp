#ifndef Matrix_HPP
#define Matrix_HPP

#include <complex>
#include <vector>
#include <memory>

class Matrix
{
    public:

        Matrix() : m_(0), n_(0), data_(std::vector<double>()) {}
        Matrix(int rows, int cols) : m_(rows), n_(cols), data_(std::vector<double>(rows*cols)) {}

        virtual ~Matrix() {}

        double* operator[](int i)
        {
            return data_.data() + i*n_;
        }

        const double* operator[](int i) const
        {
            return data_.data() + i*lda();
        }

        double* data()
        {
            return data_.data();
        }

        const double* data() const
        {
            return data_.data();
        }

        int rows() const
        {
            return m_;
        }

        int colls() const
        {
            return n_;
        }

        int lda() const
        {
            return n_;
        }

        int size() const
        {
            return m_*n_;
        }

        /*void operator- (const Matrix& m);
        void operator==(const Matrix& m);
        void operator+=(const Matrix& m);
        void operator-=(const Matrix& m);*/

    protected:
        int m_;
        int n_;
        std::vector<double> data_;
};

std::vector<std::complex<double>> eig(Matrix m);
std::vector<double> solve(Matrix A, std::vector<double> b);

#endif // Matrix_HPP