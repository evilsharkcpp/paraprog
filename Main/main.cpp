#include <RowSparseMatrix.hpp>
#include <PlotMatrix.hpp>
#include <MPITools.hpp>

#include <iostream>
#include <memory>

#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace Paraprog::MatrixLib;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (size > 0)
    {
        RowSparseMatrixPtr left;
        RowSparseMatrixPtr right;
        if (rank == 0)
        {
            left = std::make_shared<RowSparseMatrix>(RowSparseMatrix::ConvertMatrix({{1, 2, 3}, {1, 2, 3}, {1, 2, 3}}));
            right = std::make_shared<RowSparseMatrix>(RowSparseMatrix::ConvertMatrix({{1, 2, 3}, {1, 2, 3}, {1, 2, 3}}));
        }
        auto res = MPITools::MultRowSparseMatrices(left, right);
        MPI_Finalize();
        if (res != nullptr)
        {
            std::cout << res->GetPrettyMatrix();
        }
    }
    return 0;
}