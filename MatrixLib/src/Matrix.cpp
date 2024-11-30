#include "Matrix.hpp"

std::string Paraprog::MatrixLib::Matrix::GetPrettyMatrix() const
{
    std::string result;
    for (size_t i = 0; i < m_rowSize; i++)
    {
        result += "[";
        for (size_t j = 0; j < m_colSize; j++)
        {
            result += std::to_string(this->operator()(i, j)) + (j + 1 == m_colSize ? "" : ", ");
        }
        result += "]\n";
    }
    return result;
}