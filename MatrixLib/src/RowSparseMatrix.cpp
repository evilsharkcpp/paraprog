#include <RowSparseMatrix.hpp>

#include <stdexcept>
#include <algorithm>
#include <limits>
#include <utility>
#include <iostream>
#include <fstream>

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

    RowSparseMatrix RowSparseMatrix::operator*(const RowSparseMatrix &right) const
    {
        RowSparseMatrix result(m_rowSize, right.m_colSize);
        result.m_rowIndexes.resize(m_rowSize + 1);
        for (size_t i = 0; i < m_rowSize; i++)
        {
            size_t nonZeroElementsCount = 0;
            for (size_t j = 0; j < right.m_colSize; j++)
            {
                double value = 0;
                for (size_t k = m_rowIndexes[i]; k < m_rowIndexes[i + 1]; k++)
                {
                    for (size_t l = right.m_rowIndexes[m_colIndexes[k]]; l < right.m_rowIndexes[m_colIndexes[k] + 1]; l++)
                    {
                        if (right.m_colIndexes[l] == j)
                        {
                            value += m_values[k] * right.m_values[l];
                        }
                    }
                }
                if (std::abs(value) < std::numeric_limits<decltype(value)>::epsilon())
                {
                    continue;
                }
                nonZeroElementsCount++;
                result.m_values.push_back(value);
                result.m_colIndexes.push_back(j);
            }
            result.m_rowIndexes[i + 1] = result.m_rowIndexes[i] + nonZeroElementsCount;
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
        m_values = values;
        m_colIndexes = colIndexes;
        m_rowIndexes = rowIndexes;
    }

    void RowSparseMatrix::LoadProfile(std::vector<double> &&values, std::vector<size_t> &&colIndexes, std::vector<size_t> &&rowIndexes)
    {
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

    std::tuple<const std::vector<double> &, const std::vector<size_t> &, const std::vector<size_t> &> RowSparseMatrix::GetProfile() const
    {
        return std::tuple<const std::vector<double> &, const std::vector<size_t> &, const std::vector<size_t> &>(m_values, m_colIndexes, m_rowIndexes);
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

    RowSparseMatrixPtr RowSparseMatrix::LoadMatrix(const std::string &path)
    {
        std::ifstream in(path);
        if (!in.is_open())
        {
            return nullptr;
        }
        size_t m, n;
        in >> m >> n;
        const auto matrix = std::make_shared<RowSparseMatrix>(m, n);
        for (size_t i = 0; i < m + 1; i++)
        {
            size_t value = 0;
            in >> value;
            matrix->m_rowIndexes.push_back(value);
        }
        for (size_t i = 0; i < matrix->m_rowIndexes.back(); i++)
        {
            double value = 0.0;
            in >> value;
            matrix->m_values.push_back(value);
        }
        for (size_t i = 0; i < matrix->m_rowIndexes.back(); i++)
        {
            size_t value = 0;
            in >> value;
            matrix->m_colIndexes.push_back(value);
        }
        return matrix;
    }

    void RowSparseMatrix::SaveMatrix(const std::string &path)
    {
        std::ofstream out(path);
        if (!out.is_open())
        {
            return;
        }
        out << m_rowSize << " " << m_colSize << std::endl;
        for (size_t i = 0; i < m_rowIndexes.size(); i++)
        {
            out << m_rowIndexes[i] << " ";
        }
        out << std::endl;
        for (size_t i = 0; i < m_values.size(); i++)
        {
            out << m_values[i] << " ";
        }
        out << std::endl;
        for (size_t i = 0; i < m_colIndexes.size(); i++)
        {
            out << m_colIndexes[i] << " ";
        }
    }
} // namespace Paraprog::MatrixLib