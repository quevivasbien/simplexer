#pragma once

#include <algorithm>
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
    ) : tableau(canonicalTableau(objective, constraints, b)) {
        this->initBasicCols();
    }

    const Matrix<T>& matrixView() const {
        return tableau;
    }

    std::pair<std::vector<T>, T> solve() {
        while (!this->done()) {
            std::cout << this->matrixView().str() << std::endl;
            std::cout << "Basic columns: " << Vector(this->basicCols).str() << std::endl;
            std::size_t pcol = this->pivotCol();
            std::size_t prow = this->pivotRow(pcol);
            this->pivot(prow, pcol);
        }
        return this->solution();
    }
private:
    Matrix<T> tableau;
    std::vector<std::size_t> basicCols;
    
    std::vector<std::size_t> nonBasicCols() const {
        std::vector<bool> markers(this->tableau.cols, false);
        for (std::size_t i : this->basicCols) {
            markers[i] = true;
        }
        std::vector<std::size_t> nonBasics;
        for (std::size_t i = 1; i < this->tableau.cols - 1; i++) {
            if (!markers[i]) {
                nonBasics.push_back(i);
            }
        }
        return nonBasics;
    }

    void initBasicCols() {
        this->basicCols = {};
        for (std::size_t i = 1; i < this->tableau.cols - 1; i++) {
            if (this->isBasicCol(i)) {
                this->basicCols.push_back(i);
            }
        }
    }

    bool isBasicCol(std::size_t col) const {
        ColView<T> colView(this->tableau, col);
        bool foundNonZero = false;
        for (std::size_t i = 1; i < colView.size(); i++) {
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

    std::size_t pivotCol() const {
        // return the column to pivot on
        // this is selected as the nonbasic column with the most negative value in the first row
        const auto nbcols = this->nonBasicCols();
        std::cout << "Nonbasic columns: " << Vector(nbcols).str() << std::endl;
        std::size_t col = nbcols[0];
        T minVal = this->tableau.get(0, nbcols[0]);
        for (std::size_t i : nbcols) {
            T val = this->tableau.get(0, i);
            if (val < minVal) {
                minVal = val;
                col = i;
            }
        }
        return col;
    }

    std::size_t pivotRow(std::size_t pivotCol) const {
        // given a pivot column, return the row to pivot on
        // this is selected as the row with the smallest value of b_i / a_i,
        // where a_i is the value in the pivot column of row i
        std::size_t row = 1;
        T minVal = this->tableau.get(row, this->tableau.cols - 1) / this->tableau.get(row, pivotCol);
        for (std::size_t i = 2; i < this->tableau.rows; i++) {
            T val = this->tableau.get(i, this->tableau.cols - 1) / this->tableau.get(i, pivotCol);
            if (val < minVal) {
                minVal = val;
                row = i;
            }
        }
        return row;
    }

    void pivot(std::size_t pivotRow, std::size_t pivotCol) {
        // perform a pivot operation on the tableau
        // this is done by dividing the pivot row by the value in the pivot column
        // then subtracting the pivot row from all other rows, multiplied by the value in the pivot column
        // this is done in-place
        std::cout << "Pivoting on " << pivotRow << ", " << pivotCol << std::endl;
        T pivotVal = this->tableau.get(pivotRow, pivotCol);
        this->tableau.scaleRowInplace(pivotRow, 1 / pivotVal);
        for (std::size_t i = 0; i < this->tableau.rows; i++) {
            if (i == pivotRow) {
                continue;
            }
            const T scale = -this->tableau.get(i, pivotCol);
            this->tableau.addRowInplace(i, pivotRow, scale);
            std::cout << "After row " << i << ": " << std::endl;
            std::cout << this->matrixView().str() << std::endl;
        }
        // finally, add pivotCol to list of basic cols
        this->basicCols.push_back(pivotCol);
    }

    bool done() const {
        // return true if the tableau is done
        // this is the case if all values in the first row are nonnegative
        for (std::size_t i = 0; i < this->tableau.cols - 1; i++) {
            if (this->tableau.get(0, i) < 0) {
                return false;
            }
        }
        return true;
    }

    std::pair<std::vector<T>, T> solution() const {
        // return the solution to the linear program
        // assuming that the solve is done
        std::vector<T> sol(this->tableau.rows - 1, 0);
        for (std::size_t i : this->basicCols) {
            if (i > this->tableau.rows - 1) {
                continue;
            }
            // figure out which row this corresponds to
            std::size_t index = 0;
            for (std::size_t j = 1; j < this->tableau.rows; j++) {
                if (this->tableau.get(j, i) == 1) {
                    index = j;
                    break;
                }
            }
            sol[i - 1] = this->tableau.get(index, this->tableau.cols - 1);
        }
        const T objective = this->tableau.get(0, this->tableau.cols - 1);
        return { sol, objective };
    }
};