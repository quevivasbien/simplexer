#pragma once

#include "Matrix.hpp"

template <typename T>
Matrix<T> canonicalTableau(
    const std::vector<T>& objective,
    const Matrix<T>& constraints,
    const std::vector<T>& b
) {
    return BlockMatrix<T>({
        {
            Matrix<T>::fromScalar(1),
            Vector(objective).transpose() * -1,
            Matrix<T>::zeros(1, constraints.rows),
            Matrix<T>::fromScalar(0)
        },
        {
            Matrix<T>::zeros(constraints.cols, 1),
            constraints,
            Matrix<T>::identity(constraints.rows),
            Vector(b)
        }
    }).matrix();
}

template <typename T>
requires Number<T>
class Tableau {
public:
    Tableau(
        const std::vector<T>& objective,
        const Matrix<T>& constraints,
        const std::vector<T>& b
    ) : tableau(canonicalTableau(objective, constraints, b)) {}

    std::pair<std::vector<int>, std::vector<int>> basicCols() const {
        std::vector<int> basicCols;
        std::vector<int> nonBasicCols;
        for (int i = 1; i < this->tableau.cols - 1; i++) {
            if (this->isBasicCol(i)) {
                basicCols.push_back(i);
            }
            else {
                nonBasicCols.push_back(i);
            }
        }
        return { basicCols, nonBasicCols };
    }

    const Matrix<T>& matrixView() const {
        return tableau;
    }
private:
    Matrix<T> tableau;

    bool isBasicCol(int col) const {
        ColView<T> colView(this->tableau, col);
        bool foundNonZero = false;
        for (int i = 1; i < colView.size(); i++) {
            T val = colView.get(i);
            if (val != 0) {
                if (foundNonZero || val != 1) {
                    return false;
                }
                foundNonZero = true;
            }
        }
        return true;
    }

    // int pivotCol() const {
    //     // find the best column to pivot on
    // }
};