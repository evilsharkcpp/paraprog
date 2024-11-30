#pragma once

#include <RowSparseMatrix.hpp>

namespace Paraprog::MatrixLib
{
    class MPITools
    {
    public:
        static RowSparseMatrixPtr AddRowSparseMatrices(const RowSparseMatrixPtr &left, const RowSparseMatrixPtr &right);

        static RowSparseMatrixPtr MultRowSparseMatrices(const RowSparseMatrixPtr &left, const RowSparseMatrixPtr &right);
    };
}