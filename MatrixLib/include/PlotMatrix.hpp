#pragma once

#include <Matrix.hpp>
#include <vector>

namespace Paraprog::MatrixLib
{
    class PlotMatrix : public Matrix
    {
    public:
        PlotMatrix(const size_t m, const size_t n) : Matrix(m, n)
        {
            m_values = std::vector<double>(m_colSize * m_rowSize);
        }

        const double operator()(const size_t i, const size_t j) const override;

        double &operator()(const size_t i, const size_t j) override;

        PlotMatrix operator+(const Matrix &right) const;

    private:
        std::vector<double> m_values;
    };
} // namespace Paraprog::MatrixLib