#pragma once

#include <SparseMatrix.hpp>
#include <vector>

namespace Paraprog::MatrixLib
{

    class RowSparseMatrix : public SparseMatrix
    {
    public:
        RowSparseMatrix(const size_t m, const size_t n) : SparseMatrix(m, n) {}

        const double
        operator()(const size_t i, const size_t j) const override;

        double &operator()(const size_t i, const size_t j) override;

        RowSparseMatrix operator+(const RowSparseMatrix &right) const;

        RowSparseMatrix operator*(const RowSparseMatrix &right) const;

        RowSparseMatrix operator-(const RowSparseMatrix &right) const;

        RowSparseMatrix operator*(const double value) const;

        static RowSparseMatrix ConvertMatrix(const Matrix &matrix);

        static RowSparseMatrix ConvertMatrix(const std::initializer_list<std::initializer_list<double>> &matrix);

        static std::shared_ptr<RowSparseMatrix> LoadMatrix(const std::string &path);

        void SaveMatrix(const std::string &path);

        void LoadProfile(const std::vector<double> &values,
                         const std::vector<size_t> &colIndexes,
                         const std::vector<size_t> &rowIndexes);

        void LoadProfile(std::vector<double> &&values,
                         std::vector<size_t> &&colIndexes,
                         std::vector<size_t> &&rowIndexes);

        void GetProfile(std::vector<double> &values,
                        std::vector<size_t> &colIndexes,
                        std::vector<size_t> &rowIndexes) const;

        std::tuple<const std::vector<double> &, const std::vector<size_t> &, const std::vector<size_t> &> GetProfile() const;

    private:
        std::vector<double> m_values;
        std::vector<size_t> m_colIndexes;
        std::vector<size_t> m_rowIndexes;
    };

    using RowSparseMatrixPtr = std::shared_ptr<RowSparseMatrix>;

} // namespace Paraprog::MatrixLib