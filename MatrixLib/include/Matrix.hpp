#pragma once

#include <cstddef>
#include <string>
#include <memory>

namespace Paraprog::MatrixLib
{

    class Matrix
    {
    public:
        Matrix(const size_t m, const size_t n) : m_rowSize(m), m_colSize(n) {}

        virtual const double operator()(const size_t i, const size_t j) const = 0;

        virtual double &operator()(const size_t i, const size_t j) = 0;

        const size_t GetRowSize() const
        {
            return m_rowSize;
        }

        const size_t GetColSize() const
        {
            return m_colSize;
        }

        std::string GetPrettyMatrix() const;

    protected:
        size_t m_rowSize = 0;
        size_t m_colSize = 0;
    };

    using MatrixPtr = std::shared_ptr<Matrix>;
} // namespace Paraprog::MatrixLib