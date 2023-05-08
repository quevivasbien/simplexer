#include <iostream>
#include "Matrix.hpp"
#include "simplex.hpp"
#include "linalg.hpp"

int main() {
    Vector<double> objective({1, 2, 3});
    Matrix<double> constraints({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    Vector<double> b({1, 2, 3});

    BlockMatrix<double> tableau({
        { Matrix<double>::fromScalar(1), objective.transpose() * -1, Matrix<double>::fromScalar(0) },
        { Matrix<double>::zeros(constraints.cols, 1), constraints, b }
    });
    std::cout << tableau.str() << std::endl;
    std::cout << tableau.matrix().str() << std::endl;
    Tableau<double> t(objective.dataView(), constraints, b.dataView());
    std::cout << t.matrixView().str() << std::endl;

    auto [basicCols, nonBasicCols] = t.basicCols();
    std::cout << "Basic columns: ";
    for (auto i : basicCols) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
    std::cout << "Non-basic columns: ";
    for (auto i : nonBasicCols) {
        std::cout << i << " ";
    }
    std::cout << std::endl;

    return 0;
}