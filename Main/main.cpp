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
    if (argc != 5)
    {
        std::cout << "For mpi: mpiexec -n numThread ./Main f1.txt f2.txt add/mult result.txt" << std::endl;
        return -1;
    }
    if (size > 1)
    {
        RowSparseMatrixPtr left;
        RowSparseMatrixPtr right;
        if (rank == 0)
        {
            left = RowSparseMatrix::LoadMatrix(argv[1]);
            right = RowSparseMatrix::LoadMatrix(argv[2]);
        }
        auto begin = MPI_Wtime();
        RowSparseMatrixPtr res;
        if (argv[3] == "add")
        {
            res = MPITools::AddRowSparseMatrices(left, right);
        }
        else if (argv[3] == "mult")
        {
            res = MPITools::MultRowSparseMatrices(left, right);
        }

        auto end = MPI_Wtime();
        MPI_Finalize();
        if (rank == 0 && res != nullptr)
        {
            res->SaveMatrix(argv[4]);
            std::cout << "Time: " << std::fixed << std::setprecision(3) << end - begin << "s" << std::endl;
        }
    }
    else
    {
        auto left = RowSparseMatrix::LoadMatrix(argv[1]);
        auto right = RowSparseMatrix::LoadMatrix(argv[2]);
        auto begin = std::chrono::high_resolution_clock::now();
        auto res = *left * *right;
        auto end = std::chrono::high_resolution_clock::now();
        res.SaveMatrix(argv[4]);
        std::cout << "Time: " << std::fixed << std::setprecision(5) << std::chrono::duration<double>(end - begin).count() << "s" << std::endl;
    }
    return 0;
}