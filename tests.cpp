#include <iostream>
#include <vector>
#include "Matrix.hpp"
#include "linalg.hpp"
#include "linprog.hpp"
#include "utils.hpp"

const double RTOL = 1e-6;
const double ATOL = 1e-8;

#define IS_TRUE(x) { if (!(x)) { std::cout << __FUNCTION__ << " failed: " << " on line " << __LINE__ << std::endl; return 1; } }

int testLinalg() {
    Matrix<double> m1({{0, 1, 2}, {3, 4, 5}, {6, 7, 8}});
    linalg::rowReductionInplace(m1);
    IS_TRUE(m1.approxEqual(Matrix<double>({{6, 7, 8}, {0, 1, 2}, {0, 0, 0}}), RTOL, ATOL));

    Matrix<double> m2({{0, 1, 2, 3}, {4, 5, 6, 7}, {8, 9, 10, 11}});
    IS_TRUE(linalg::mult(m2, m2.transpose()).approxEqual(Matrix<double>({{14, 38, 62}, {38, 126, 214}, {62, 214, 366}}), RTOL, ATOL));
    return 0;
}

int testLinprog() {
    const auto s1 = linprog::Tableau<double>(
        std::vector<double>({4, 1, 4}),
        Matrix<double>({{2, 1, 1}, {1, 2, 3}, {2, 2, 1}}),
        std::vector<double>({2, 4, 8})
    ).solve();
    IS_TRUE(Vector(s1.maximizer).approxEqual(Vector<double>({0.4, 0, 1.2}), RTOL, ATOL));
    IS_TRUE(utils::approxEqual(s1.maximum, 6.4, RTOL, ATOL));

    const auto s2 = linprog::Tableau<double>(
        std::vector<double>({1, -4}),
        Matrix<double>({{-3, 1}, {1, 2}}),
        std::vector<double>({6, 4})
    ).solve();
    IS_TRUE(Vector(s2.maximizer).approxEqual(Vector<double>({4, 0}), RTOL, ATOL));
    IS_TRUE(utils::approxEqual(s2.maximum, 4., RTOL, ATOL));
    
    const auto s3 = linprog::Tableau<double>(
        std::vector<double>({5, 6, 4, 4}),
        Matrix<double>({{9, 9, 5, 7}, {9, 7, 3, 8}, {8, 6, 6, 3}, {9, 2, 2, 9}}),
        std::vector<double>({3, 4, 4, 4})
    ).solve();
    IS_TRUE(Vector(s3.maximizer).approxEqual(Vector<double>({0, 0, 0.6, 0}), RTOL, ATOL));
    IS_TRUE(utils::approxEqual(s3.maximum, 2.4, RTOL, ATOL));

    return 0;
}

int main() {
    int linalg = testLinalg();
    int linprog = testLinprog();
    int res = linalg || linprog;
    if (res == 0) {
        std::cout << "All tests passed!" << std::endl;
    }
    return res;
}