#pragma once

#include <Matrix.hpp>

namespace Paraprog::MatrixLib
{

    class SparseMatrix : public Matrix
    {
    public:
        SparseMatrix(const size_t m, const size_t n) : Matrix(m, n) {}

        const double operator()(const size_t i, const size_t j) const override = 0;

        double &operator()(const size_t i, const size_t j) override = 0;
    };

} // namespace Paraprog::MatrixLib