#include <algorithm>
#include "SparseMatrix.hpp"

SparseMatrix::SparseMatrix(const Matrix& M) : m_(M.rows()), n_(M.colls()), data_(std::vector<Triplet>())
{
    for (int i = 0; i < M.rows(); i++)
    {
        for (int j = 0; j < M.colls(); j++)
        {
            if (M[i][j] != 0.0)
            {
                addNew(Triplet(i, j, M[i][j]));
            }
        }
    }

    sort();
}

double SparseMatrix::operator() (int i, int j) const
{
    auto it = std::find(data_.begin(), data_.end(), Triplet(i, j, 0.0));

    if (it != data_.end())
    {
        return it->val;
    }
    else
    {
        return 0.0;
    }
}

double& SparseMatrix::operator()(int i, int j)
{
    auto it = std::find(data_.begin(), data_.end(), Triplet(i, j, 0.0));

    if (it != data_.end())
    {
        return it->val;
    }
    else
    {
        data_.push_back(Triplet(i, j, 0.0));
        return data_.back().val;
    }
}

void SparseMatrix::sort()
{
    data_.erase(std::remove_if(data_.begin(), 
                               data_.end(),
                               [](const Triplet& t) { return t.val == 0.0; }),
                data_.end());

    std::sort(data_.begin(), data_.end());
}

void SparseMatrix::addNew(Triplet t)
{
    data_.push_back(t);
}

int SparseMatrix::getM() const
{
    return m_;
}

int SparseMatrix::getN() const
{
    return n_;
}

int SparseMatrix::getNnz() const
{
    return data_.size();
}

const std::vector<Triplet>& SparseMatrix::getTriplets() const
{
    return data_;
}