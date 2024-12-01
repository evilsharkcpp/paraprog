#include "MPITools.hpp"

#include <mpi.h>

#include <iostream>
#include <limits>

namespace Paraprog::MatrixLib
{

    static void SendMatrixPart(const size_t threadNum, const std::vector<double> &value, const std::vector<size_t> &row, const std::vector<size_t> &col)
    {
        const size_t vSize = value.size();
        const size_t rowSize = row.size();
        const size_t colSize = col.size();
        const int size = 3 * sizeof(size_t) + vSize * sizeof(double) + rowSize * sizeof(size_t) + colSize * sizeof(size_t);
        char pack[size];
        int bufPos = 0;
        MPI_Pack(&vSize, 1, MPI_UNSIGNED_LONG, &pack, size, &bufPos, MPI_COMM_WORLD);
        MPI_Pack(&rowSize, 1, MPI_UNSIGNED_LONG, &pack, size, &bufPos, MPI_COMM_WORLD);
        MPI_Pack(&colSize, 1, MPI_UNSIGNED_LONG, &pack, size, &bufPos, MPI_COMM_WORLD);
        MPI_Pack(value.data(), vSize, MPI_DOUBLE, &pack, size, &bufPos, MPI_COMM_WORLD);
        MPI_Pack(row.data(), rowSize, MPI_UNSIGNED_LONG, &pack, size, &bufPos, MPI_COMM_WORLD);
        MPI_Pack(col.data(), colSize, MPI_UNSIGNED_LONG, &pack, size, &bufPos, MPI_COMM_WORLD);
        MPI_Send(&size, 1, MPI_INT, threadNum, 2, MPI_COMM_WORLD);
        MPI_Send(pack, bufPos, MPI_BYTE, threadNum, 1, MPI_COMM_WORLD);
    }

    static std::tuple<const std::vector<double>, const std::vector<size_t>, const std::vector<size_t>> GetMatrixProfile(const size_t source)
    {
        size_t vSize;
        size_t rowSize;
        size_t colSize;
        MPI_Status status;
        int count;
        MPI_Recv(&count, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
        char pack[count];
        MPI_Recv(pack, count, MPI_BYTE, source, 1, MPI_COMM_WORLD, &status);
        int bufPos = 0;
        MPI_Unpack(pack, count, &bufPos, &vSize, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
        MPI_Unpack(pack, count, &bufPos, &rowSize, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
        MPI_Unpack(pack, count, &bufPos, &colSize, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
        std::vector<double> valuesPart(vSize);
        std::vector<size_t> rowPart(rowSize);
        std::vector<size_t> colPart(colSize);
        MPI_Unpack(pack, count, &bufPos, valuesPart.data(), vSize, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack(pack, count, &bufPos, rowPart.data(), rowSize, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
        MPI_Unpack(pack, count, &bufPos, colPart.data(), colSize, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

        return std::tuple<const std::vector<double>, const std::vector<size_t>, const std::vector<size_t>>(valuesPart, rowPart, colPart);
    }

    RowSparseMatrixPtr MPITools::AddRowSparseMatrices(const RowSparseMatrixPtr &left, const RowSparseMatrixPtr &right)
    {
        int rank, size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        size_t maxCount;
        size_t rowSize;
        RowSparseMatrixPtr result = nullptr;
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        if (rank == 0)
        {
            if (left->GetRowSize() != right->GetRowSize() || left->GetColSize() != right->GetColSize())
            {
                return nullptr;
            }
            rowSize = left->GetRowSize();
            result = std::make_shared<RowSparseMatrix>(left->GetRowSize(), left->GetColSize());
            rows.resize(left->GetRowSize() + 1);
        }
        MPI_Bcast(&rowSize, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        for (size_t rowIndex = 0; rowIndex < rowSize; rowIndex += (size - 1))
        {
            for (size_t threadNum = 1; threadNum < size; threadNum++)
            {
                if (rank == 0)
                {
                    size_t i = rowIndex + threadNum - 1;
                    const auto &[lvalues, lcol, lrow] = left->GetProfile();
                    const auto &[rvalues, rcol, rrow] = right->GetProfile();
                    size_t colSize = left->GetColSize();
                    MPI_Send(&colSize, 1, MPI_UNSIGNED_LONG, threadNum, 0, MPI_COMM_WORLD);
                    const size_t lElementsInRow = lrow[i + 1] - lrow[i];
                    auto lvaluesPart = std::vector<double>(lvalues.begin() + lrow[i], lvalues.begin() + lrow[i] + lElementsInRow);
                    auto lRowPart = std::vector<size_t>{0, lElementsInRow};
                    auto lColPart = std::vector<size_t>(lcol.begin() + lrow[i], lcol.begin() + lrow[i] + lElementsInRow);
                    SendMatrixPart(threadNum, lvaluesPart, lRowPart, lColPart);

                    const size_t rElementsInRow = rrow[i + 1] - rrow[i];
                    auto rvaluesPart = std::vector<double>(rvalues.begin() + rrow[i], rvalues.begin() + rrow[i] + rElementsInRow);
                    auto rRowPart = std::vector<size_t>{0, rElementsInRow};
                    auto rColPart = std::vector<size_t>(rcol.begin() + rrow[i], rcol.begin() + rrow[i] + rElementsInRow);
                    SendMatrixPart(threadNum, rvaluesPart, rRowPart, rColPart);
                }
                if (threadNum == rank)
                {
                    size_t colSize;
                    MPI_Status status;
                    MPI_Recv(&colSize, 1, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    const auto &[lvaluesPart, lRowPart, lColPart] = GetMatrixProfile(0);
                    const auto &[rvaluesPart, rRowPart, rColPart] = GetMatrixProfile(0);
                    RowSparseMatrix left(1, colSize);
                    left.LoadProfile(std::move(lvaluesPart), std::move(lColPart), std::move(lRowPart));
                    RowSparseMatrix right(1, colSize);
                    right.LoadProfile(std::move(rvaluesPart), std::move(rColPart), std::move(rRowPart));
                    auto result = left + right;
                    const auto &[values, col, row] = result.GetProfile();
                    SendMatrixPart(0, values, row, col);
                }
            }
            std::cout << "Data was send\n";
            if (rank == 0)
            {
                for (size_t threadNum = 1; threadNum < size; threadNum++)
                {
                    size_t i = rowIndex + threadNum - 1;
                    const auto &[valuesPart, rowsPart, colsPart] = GetMatrixProfile(threadNum);
                    rows[i + 1] = rows[i] + rowsPart.back();
                    values.insert(values.end(), valuesPart.begin(), valuesPart.end());
                    cols.insert(cols.end(), colsPart.begin(), colsPart.end());
                }
                std::cout << "Data was collected\n";
            }
        }
        if (rank == 0)
        {
            result->LoadProfile(std::move(values), std::move(cols), std::move(rows));
        }
        return result;
    }

    RowSparseMatrixPtr MPITools::MultRowSparseMatrices(const RowSparseMatrixPtr &left, const RowSparseMatrixPtr &right)
    {
        int rank, size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        size_t maxCount;
        size_t rowSize;
        size_t colSize;
        RowSparseMatrixPtr result = nullptr;
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        if (rank == 0)
        {
            if (left->GetRowSize() != right->GetRowSize() || left->GetColSize() != right->GetColSize())
            {
                return nullptr;
            }
            rowSize = left->GetRowSize();
            colSize = right->GetColSize();
            result = std::make_shared<RowSparseMatrix>(left->GetRowSize(), right->GetColSize());
            rows.resize(left->GetRowSize() + 1);
        }
        MPI_Bcast(&rowSize, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&colSize, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        for (size_t rowIndex = 0; rowIndex < rowSize; rowIndex++)
        {
            size_t nonZeroElementsCount = 0;
            for (size_t colIndex = 0; colIndex < colSize; colIndex += (size - 1))
            {
                for (size_t threadNum = 1; threadNum < size; threadNum++)
                {
                    if (rank == 0)
                    {
                        size_t j = colIndex + threadNum - 1;
                        const auto &[lvalues, lcol, lrow] = left->GetProfile();
                        const auto &[rvalues, rcol, rrow] = right->GetProfile();
                        size_t lcolSize = left->GetColSize();
                        MPI_Send(&lcolSize, 1, MPI_UNSIGNED_LONG, threadNum, 0, MPI_COMM_WORLD);
                        const size_t lElementsInRow = lrow[rowIndex + 1] - lrow[rowIndex];
                        auto lvaluesPart = std::vector<double>(lvalues.begin() + lrow[rowIndex], lvalues.begin() + lrow[rowIndex] + lElementsInRow);
                        auto lRowPart = std::vector<size_t>{0, lElementsInRow};
                        auto lColPart = std::vector<size_t>(lcol.begin() + lrow[rowIndex], lcol.begin() + lrow[rowIndex] + lElementsInRow);
                        SendMatrixPart(threadNum, lvaluesPart, lRowPart, lColPart);
                        auto rvaluesPart = std::vector<double>();
                        auto rRowPart = std::vector<size_t>(lcolSize + 1, 0);
                        auto rColPart = std::vector<size_t>();
                        for (size_t l = 0; l < lcolSize; l++)
                        {
                            for (size_t k = rrow[l]; k < rrow[l + 1]; k++)
                            {
                                if (rcol[k] == j)
                                {
                                    rvaluesPart.push_back(rvalues[k]);
                                    rRowPart[l + 1] = rRowPart[l] + 1;
                                    rColPart.push_back(0);
                                    break;
                                }
                            }
                        }
                        SendMatrixPart(threadNum, rvaluesPart, rRowPart, rColPart);
                    }
                    if (threadNum == rank)
                    {
                        size_t lcolSize;
                        MPI_Status status;
                        MPI_Recv(&lcolSize, 1, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                        const auto &[lvaluesPart, lRowPart, lColPart] = GetMatrixProfile(0);
                        const auto &[rvaluesPart, rRowPart, rColPart] = GetMatrixProfile(0);
                        RowSparseMatrix left(1, lcolSize);
                        left.LoadProfile(std::move(lvaluesPart), std::move(lColPart), std::move(lRowPart));
                        RowSparseMatrix right(lcolSize, 1);
                        right.LoadProfile(std::move(rvaluesPart), std::move(rColPart), std::move(rRowPart));
                        auto result = left * right;
                        const auto &[values, col, row] = result.GetProfile();
                        SendMatrixPart(0, values, row, col);
                    }
                }
                std::cout << "Data was send\n";
                if (rank == 0)
                {
                    for (size_t threadNum = 1; threadNum < size; threadNum++)
                    {
                        size_t j = colIndex + threadNum - 1;
                        const auto &[valuesPart, rowsPart, colsPart] = GetMatrixProfile(threadNum);
                        auto val = (valuesPart.empty()) ? 0.0 : valuesPart.front();
                        if (std::abs(val) < std::numeric_limits<decltype(val)>::epsilon())
                        {
                            continue;
                        }
                        nonZeroElementsCount++;
                        values.push_back(val);
                        cols.push_back(j);
                    }
                    std::cout << "Data was collected\n";
                }
            }
            if (rank == 0)
            {
                rows[rowIndex + 1] = rows[rowIndex] + nonZeroElementsCount;
            }
        }
        if (rank == 0)
        {
            result->LoadProfile(std::move(values), std::move(cols), std::move(rows));
        }
        return result;
    }
}