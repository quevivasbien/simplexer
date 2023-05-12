#pragma once

#include <algorithm>
#include <limits>
#include "Matrix.hpp"

namespace linprog {

template <typename T>
Matrix<T> canonicalTableau(
    const std::vector<T>& c,
    const Matrix<T>& A,
    const std::vector<T>& b
) {
    return BlockMatrix<T>({
        {
            Matrix<T>::fromScalar(1),
            Vector(c).transpose() * -1,
            Matrix<T>::zeros(1, A.rows),
            Matrix<T>::fromScalar(0)
        },
        {
            Matrix<T>::zeros(A.cols, 1),
            A,
            Matrix<T>::identity(A.rows),
            Vector(b)
        }
    }).matrix();
}

enum Status {
    OK,
    NO_SOLUTION,
};

template <typename T>
struct Solution {
    Status status;
    std::vector<T> maximizer;
    T maximum;
};

template <typename T>
requires Number<T>
class Tableau {
public:
    Tableau(
        const std::vector<T>& c,
        const Matrix<T>& A,
        const std::vector<T>& b
    ) : tableau(canonicalTableau(c, A, b)) {
        for (std::size_t i = A.cols + 1; i < A.cols + A.rows + 1; i++) {
            this->basicCols.push_back(i);
        }
    }

    const Matrix<T>& matrixView() const {
        return tableau;
    }

    Solution<T> solve() {
        while (!this->done()) {
            #ifdef DEBUG
            std::cout << "Tableau:\n" << this->tableau.str() << std::endl;
            std::cout << "basic columns: " << Vector(this->basicCols).str() << std::endl;
            #endif
            std::size_t pcol = this->pivotCol();
            std::size_t prow = this->pivotRow(pcol);
            if (prow == 0) {
                return Solution<T>(NO_SOLUTION, {}, 0);
            }
            this->pivot(prow, pcol);
        }
        return this->solution();
    }
private:
    Matrix<T> tableau;
    // there is one entry for each row of A (each constraint)
    // each entry gives the column index of the basic variable for that row
    std::vector<std::size_t> basicCols;

    std::vector<std::size_t> nonBasicCols() const {
        std::vector<bool> markers(this->tableau.cols - 1, false);
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

    std::size_t pivotCol() const {
        // return the column to pivot on
        // this is selected as the nonbasic column with the most negative value in the first row
        const auto nbcols = this->nonBasicCols();
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
        // where a_i is the value in the pivot column of row i, filtered to only positive values
        // returning 0 means no feasible solution
        std::size_t row = 0;
        T minVal = std::numeric_limits<T>::max();
        for (std::size_t i = 1; i < this->tableau.rows; i++) {
            T a = this->tableau.get(i, pivotCol);
            if (a <= 0) {
                continue;
            }
            T val = this->tableau.get(i, this->tableau.cols - 1) / a;
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
        #ifdef DEBUG
        std::cout << "Pivoting on row " << pivotRow << ", column " << pivotCol << std::endl;
        #endif
        T pivotVal = this->tableau.get(pivotRow, pivotCol);
        this->tableau.scaleRowInplace(pivotRow, 1 / pivotVal);
        for (std::size_t i = 0; i < this->tableau.rows; i++) {
            if (i == pivotRow) {
                continue;
            }
            const T scale = -this->tableau.get(i, pivotCol);
            this->tableau.addRowInplace(i, pivotRow, scale);
        }
        // set entering variable to basic
        this->basicCols[pivotRow - 1] = pivotCol;
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

    Solution<T> solution() const {
        // return the solution to the linear program
        // assuming that the solve is done
        std::vector<T> sol(this->tableau.rows - 1, 0);
        #ifdef DEBUG
        std::cout << "Tableau:\n" << this->tableau.str() << std::endl;
        std::cout << "basic columns: " << Vector(this->basicCols).str() << std::endl;
        #endif
        for (std::size_t i = 0; i < this->basicCols.size(); i++) {
            if (this->basicCols[i] > this->tableau.rows - 1) {
                continue;
            }
            sol[this->basicCols[i] - 1] = this->tableau.get(i + 1, this->tableau.cols - 1);
        }
        const T objective_max = this->tableau.get(0, this->tableau.cols - 1);
        return Solution(OK, sol, objective_max);
    }
};

} // namespace linprog