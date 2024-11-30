#include <RowSparseMatrix.hpp>

#include <stdexcept>
#include <algorithm>
#include <limits>
#include <utility>

namespace Paraprog::MatrixLib
{

    const double RowSparseMatrix::operator()(const size_t i, const size_t j) const
    {
        for (int index = m_rowIndexes.at(i); index < m_rowIndexes.at(i + 1); index++)
        {
            if (m_colIndexes[index] == j)
            {
                return m_values[index];
            }
        }
        return 0.0;
    }

    double &RowSparseMatrix::operator()(const size_t i,
                                        const size_t j)
    {
        for (int index = m_rowIndexes.at(i); index < m_rowIndexes.at(i + 1); index++)
        {
            if (m_colIndexes[index] == j)
            {
                return m_values[index];
            }
        }
        throw std::out_of_range("Try to get element out of range or element is zero");
    }

    RowSparseMatrix RowSparseMatrix::operator+(const RowSparseMatrix &right) const
    {
        if (m_colSize != right.GetColSize() || m_rowSize != right.GetRowSize())
        {
            throw std::out_of_range("Matrix sizes not equals");
        }
        RowSparseMatrix result(m_rowSize, m_colSize);
        result.m_rowIndexes.reserve(m_rowSize + 1);
        result.m_rowIndexes.push_back(0);
        for (size_t i = 0; i < m_rowSize; i++)
        {
            size_t nonZeroElementsCount = 0;
            std::vector<size_t> nonZeroColns;
            for (size_t j = 0; j < m_colSize; j++)
            {
                const auto value = this->operator()(i, j) + right(i, j);
                if (std::abs(value) < std::numeric_limits<decltype(value)>::epsilon())
                {
                    continue;
                }
                nonZeroElementsCount++;
                nonZeroColns.push_back(j);
                result.m_values.push_back(value);
            }
            for (const auto &colIndex : nonZeroColns)
            {
                result.m_colIndexes.push_back(colIndex);
            }
            result.m_rowIndexes.push_back(result.m_rowIndexes.at(i) + nonZeroElementsCount);
        }
        return result;
    }

    RowSparseMatrix RowSparseMatrix::operator-(const RowSparseMatrix &right) const
    {
        if (m_colSize != right.GetColSize() || m_rowSize != right.GetRowSize())
        {
            throw std::out_of_range("Matrix sizes not equals");
        }
        RowSparseMatrix result(m_rowSize, m_colSize);
        result.m_rowIndexes.reserve(m_rowSize + 1);
        result.m_rowIndexes.push_back(0);
        for (size_t i = 0; i < m_rowSize; i++)
        {
            size_t nonZeroElementsCount = 0;
            std::vector<size_t> nonZeroColns;
            for (size_t j = 0; j < m_colSize; j++)
            {
                const auto value = this->operator()(i, j) - right(i, j);
                if (std::abs(value) < std::numeric_limits<decltype(value)>::epsilon())
                {
                    continue;
                }
                nonZeroElementsCount++;
                nonZeroColns.push_back(j);
                result.m_values.push_back(value);
            }
            for (const auto &colIndex : nonZeroColns)
            {
                result.m_colIndexes.push_back(colIndex);
            }
            result.m_rowIndexes.push_back(result.m_rowIndexes.at(i) + nonZeroElementsCount);
        }
        return result;
    }

    RowSparseMatrix RowSparseMatrix::operator*(const double value) const
    {
        RowSparseMatrix result(m_rowSize, m_colSize);
        result.LoadProfile(m_values, m_colIndexes, m_rowIndexes);
        for (size_t i = 0; i + 1 < m_rowIndexes.size(); i++)
        {
            for (int index = m_rowIndexes.at(i); index < m_rowIndexes.at(i + 1); index++)
            {
                result.m_values[index] = m_values[index] * value;
            }
        }
        return result;
    }

    void RowSparseMatrix::LoadProfile(const std::vector<double> &values, const std::vector<size_t> &colIndexes, const std::vector<size_t> &rowIndexes)
    {
        if (rowIndexes.size() != m_rowSize + 1 || *std::max_element(colIndexes.begin(), colIndexes.end()) >= m_colSize || values.size() != rowIndexes.back())
        {
            throw std::logic_error("Incorrect profile for matrix");
        }
        m_values = values;
        m_colIndexes = colIndexes;
        m_rowIndexes = rowIndexes;
    }

    void RowSparseMatrix::GetProfile(std::vector<double> &values, std::vector<size_t> &colIndexes, std::vector<size_t> &rowIndexes) const
    {
        values = std::vector<double>(m_values);
        colIndexes = std::vector<size_t>(m_colIndexes);
        rowIndexes = std::vector<size_t>(m_rowIndexes);
    }

    RowSparseMatrix RowSparseMatrix::ConvertMatrix(const Matrix &matrix)
    {
        RowSparseMatrix result(matrix.GetRowSize(), matrix.GetColSize());
        result.m_rowIndexes.reserve(matrix.GetRowSize() + 1);
        result.m_rowIndexes.push_back(0);
        for (size_t i = 0; i < matrix.GetRowSize(); i++)
        {
            size_t nonZeroElementsCount = 0;
            std::vector<size_t> nonZeroColns;
            for (size_t j = 0; j < matrix.GetColSize(); j++)
            {
                const auto value = matrix(i, j);
                if (std::abs(value) < std::numeric_limits<decltype(value)>::epsilon())
                {
                    continue;
                }
                nonZeroElementsCount++;
                nonZeroColns.push_back(j);
                result.m_values.push_back(value);
            }
            for (const auto &colIndex : nonZeroColns)
            {
                result.m_colIndexes.push_back(colIndex);
            }
            result.m_rowIndexes.push_back(result.m_rowIndexes.at(i) + nonZeroElementsCount);
        }
        return result;
    }

    RowSparseMatrix RowSparseMatrix::ConvertMatrix(const std::initializer_list<std::initializer_list<double>> &matrix)
    {
        RowSparseMatrix result(matrix.size(), matrix.begin()->size());
        result.m_rowIndexes.reserve(matrix.size() + 1);
        result.m_rowIndexes.push_back(0);
        for (auto [row, i] = std::pair<decltype(matrix.begin()), size_t>(matrix.begin(), 0); row != matrix.end(); row++, i++)
        {
            size_t nonZeroElementsCount = 0;
            std::vector<size_t> nonZeroColns;
            for (auto [elem, j] = std::pair<decltype(row->begin()), size_t>(row->begin(), 0); elem != row->end(); elem++, j++)
            {
                const auto value = *elem;
                if (std::abs(value) < std::numeric_limits<decltype(value)>::epsilon())
                {
                    continue;
                }
                nonZeroElementsCount++;
                nonZeroColns.push_back(j);
                result.m_values.push_back(value);
            }
            for (const auto &colIndex : nonZeroColns)
            {
                result.m_colIndexes.push_back(colIndex);
            }
            result.m_rowIndexes.push_back(result.m_rowIndexes.at(i) + nonZeroElementsCount);
        }
        return result;
    }
} // namespace Paraprog::MatrixLib