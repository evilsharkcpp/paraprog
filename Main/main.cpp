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
        if (rank == 0)
        {
            std::cout << "For mpi: mpiexec -n numThread ./Main f1.txt f2.txt add/mult result.txt" << std::endl;
        }
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
            if (left == nullptr || right == nullptr)
            {
                std::cout << "Not found input files" << std::endl;
                return -1;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        auto begin = MPI_Wtime();
        RowSparseMatrixPtr res;
        if (std::string(argv[3]) == "add")
        {
            res = MPITools::AddRowSparseMatrices(left, right);
        }
        else if (std::string(argv[3]) == "mult")
        {
            res = MPITools::MultRowSparseMatrices(left, right);
        }
        else
        {
            if (rank == 0)
            {
                std::cout << "For mpi: mpiexec -n numThread ./Main f1.txt f2.txt add/mult result.txt" << std::endl;
            }
            return -1;
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
        RowSparseMatrixPtr res;
        if (std::string(argv[3]) == "add")
        {
            res = std::make_shared<RowSparseMatrix>(*left + *right);
        }
        else if (std::string(argv[3]) == "mult")
        {
            res = std::make_shared<RowSparseMatrix>(*left * *right);
        }
        else
        {
            return -1;
        }
        auto end = std::chrono::high_resolution_clock::now();
        res->SaveMatrix(argv[4]);
        std::cout << "Time: " << std::fixed << std::setprecision(5) << std::chrono::duration<double>(end - begin).count() << "s" << std::endl;
    }
    return 0;
}