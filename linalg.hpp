#pragma once

#include "Matrix.hpp"

namespace linalg {

template <typename T>
Matrix<T> mult(Matrix<T> matrix1, Matrix<T> matrix2) {
    assert(matrix1.cols == matrix2.rows);
    std::vector<T> data;
    data.reserve(matrix1.rows * matrix2.cols);
    for (std::size_t i = 0; i < matrix1.rows; i++) {
        for (std::size_t j = 0; j < matrix2.cols; j++) {
            T sum = 0;
            for (std::size_t k = 0; k < matrix1.cols; k++) {
                sum += matrix1.get(i, k) * matrix2.get(k, j);
            }
            data.push_back(sum);
        }
    }
    return Matrix<T>(data, matrix1.rows, matrix2.cols);
}

template <typename T>
void rowReductionInplace(Matrix<T>& matrix) {
    std::size_t pivotRow = 0;
    std::size_t pivotCol = 0;
    while (pivotRow < matrix.rows && pivotCol < matrix.cols) {
        std::size_t maxRow = pivotRow;
        for (std::size_t i = pivotRow + 1; i < matrix.rows; i++) {
            if (std::abs(matrix.get(i, pivotCol)) > std::abs(matrix.get(maxRow, pivotCol))) {
                maxRow = i;
            }
        }
        if (matrix.get(maxRow, pivotCol) == 0) {
            pivotCol++;
            continue;
        }
        matrix.swapRowsInplace(pivotRow, maxRow);
        for (std::size_t i = pivotRow + 1; i < matrix.rows; i++) {
            T factor = matrix.get(i, pivotCol) / matrix.get(pivotRow, pivotCol);
            matrix.addRowInplace(i, pivotRow, -factor);
        }
        pivotRow++;
        pivotCol++;
    }
}

template <typename T>
Matrix<T> rowReduction(Matrix<T> matrix) {
    rowReductionInplace(matrix);
    return matrix;
}

template <typename T>
int rank(Matrix<T> matrix) {
    rowReductionInplace(matrix);
    std::size_t rank = 0;
    for (std::size_t i = 0; i < matrix.rows; i++) {
        for (std::size_t j = 0; j < matrix.cols; j++) {
            if (matrix.get(i, j) != 0) {
                rank++;
                break;
            }
        }
    }
    return int(rank);
}

} // namespace linalg