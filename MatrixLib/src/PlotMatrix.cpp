#include <PlotMatrix.hpp>

#include <stdexcept>

namespace Paraprog::MatrixLib
{
    const double PlotMatrix::operator()(const size_t i, const size_t j) const
    {
        return m_values.at(m_rowSize * i + j);
    }

    double &PlotMatrix::operator()(const size_t i, const size_t j)
    {
        return m_values.at(m_rowSize * i + j);
    }

    PlotMatrix PlotMatrix::operator+(const Matrix &right) const
    {
        if (m_colSize != right.GetColSize() || m_rowSize != right.GetRowSize())
        {
            throw std::out_of_range("Matrix sizes not equals");
        }
        PlotMatrix result(m_rowSize, m_colSize);
        for (size_t i = 0; i < m_rowSize; i++)
        {
            for (size_t j = 0; j < m_colSize; j++)
            {
                result(i, j) = this->operator()(i, j) + right(i, j);
            }
        }
        return result;
    }
}