#pragma once

#include "Matrix.hpp"

namespace linalg {

template <typename T>
Matrix<T> mult(Matrix<T> matrix1, Matrix<T> matrix2) {
    assert(matrix1.cols == matrix2.rows);
    std::vector<T> data;
    data.reserve(matrix1.rows * matrix2.cols);
    for (int i = 0; i < matrix1.rows; i++) {
        for (int j = 0; j < matrix2.cols; j++) {
            T sum = 0;
            for (int k = 0; k < matrix1.cols; k++) {
                sum += matrix1.get(i, k) * matrix2.get(k, j);
            }
            data.push_back(sum);
        }
    }
    return Matrix<T>(data, matrix1.rows, matrix2.cols);
}

template <typename T>
void rowReductionInplace(Matrix<T>& matrix) {
    int pivotRow = 0;
    int pivotCol = 0;
    while (pivotRow < matrix.rows && pivotCol < matrix.cols) {
        int maxRow = pivotRow;
        for (int i = pivotRow + 1; i < matrix.rows; i++) {
            if (std::abs(matrix.get(i, pivotCol)) > std::abs(matrix.get(maxRow, pivotCol))) {
                maxRow = i;
            }
        }
        if (matrix.get(maxRow, pivotCol) == 0) {
            pivotCol++;
            continue;
        }
        matrix.swapRowsInplace(pivotRow, maxRow);
        for (int i = pivotRow + 1; i < matrix.rows; i++) {
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
    int rank = 0;
    for (int i = 0; i < matrix.rows; i++) {
        for (int j = 0; j < matrix.cols; j++) {
            if (matrix.get(i, j) != 0) {
                rank++;
                break;
            }
        }
    }
    return rank;
}

} // namespace linalg