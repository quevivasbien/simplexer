#include <iostream>
#include "Matrix.hpp"
#include "simplex.hpp"
#include "linalg.hpp"

int main() {
    std::vector<double> objective({4, 1, 4});
    Matrix<double> constraints({{2, 1, 1}, {1, 2, 3}, {2, 2, 1}});
    std::vector<double> b({2, 4, 8});

    Tableau<double> t(objective, constraints, b);
    std::cout << t.matrixView().str() << std::endl;

    const auto [solution, maximum] = t.solve();
    std::cout << "Solution: " << Vector(solution).str() << std::endl;
    std::cout << "Maximum: " << maximum << std::endl;

    return 0;
}