#include <iostream>
#include "Matrix.hpp"
#include "linprog.hpp"

int main() {
    std::vector<double> objective({4, 1, 4});
    Matrix<double> constraints({{2, 1, 1}, {1, 2, 3}, {2, 2, 1}});
    std::vector<double> b({2, 4, 8});

    linprog::Tableau<double> t(objective, constraints, b);
    std::cout << t.matrixView().str() << std::endl;

    const auto s = t.solve();
    std::cout << "Solution: " << Vector(s.solution).str() << std::endl;
    std::cout << "Maximum: " << s.maximum << std::endl;

    return 0;
}