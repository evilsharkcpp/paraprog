#include <RowSparseMatrix.hpp>
#include <MPITools.hpp>

#include <iostream>
#include <memory>
#include <chrono>
#include <iomanip>

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
    if (size > 1)
    {
        RowSparseMatrixPtr left;
        RowSparseMatrixPtr right;
        if (rank == 0)
        {
            left = std::make_shared<RowSparseMatrix>(RowSparseMatrix::ConvertMatrix({{1, 2, 3}, {1, 2, 3}, {1, 2, 3}}));
            right = std::make_shared<RowSparseMatrix>(RowSparseMatrix::ConvertMatrix({{1, 2, 3}, {1, 2, 3}, {1, 2, 3}}));
        }
        auto begin = MPI_Wtime();
        auto res = MPITools::MultRowSparseMatrices(left, right);
        auto end = MPI_Wtime();
        MPI_Finalize();
        if (res != nullptr)
        {
            std::cout << res->GetPrettyMatrix();
        }
        if (rank == 0)
        {
            std::cout << "Time: " << std::fixed << std::setprecision(3) << end - begin << "s" << std::endl;
        }
    }
    else
    {
        auto left = std::make_shared<RowSparseMatrix>(RowSparseMatrix::ConvertMatrix({{1, 2, 3}, {1, 2, 3}, {1, 2, 3}}));
        auto right = std::make_shared<RowSparseMatrix>(RowSparseMatrix::ConvertMatrix({{1, 2, 3}, {1, 2, 3}, {1, 2, 3}}));
        auto begin = std::chrono::high_resolution_clock::now();
        auto res = *left * *right;
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << res.GetPrettyMatrix();
        std::cout << "Time: " << std::fixed << std::setprecision(5) << std::chrono::duration<double>(end - begin).count() << "s" << std::endl;
    }
    return 0;
}