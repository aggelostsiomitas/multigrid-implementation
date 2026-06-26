#ifndef GAUSSSEIDEL_HPP
#define GAUSSSEIDEL_HPP

#include<vector>
#include "matrix.hpp"

/**
 * @brief Gauss–Seidel smoother.
 *
 * Performs \c GaussSeidel iterations to solve A x = b.
 *
 * @param A  System matrix.
 * @param DL Lower triangular part of A.
 * @param x  Solution vector (updated in place).
 * @param b  Right-hand side.
 * @param nu Number of iterations.
 */
void GaussSeidel(
    const Matrix& A,
    const Matrix& DL,
    Matrix& x,
    const Matrix& b,
    int nu);
#endif // GAUSSSEIDEL_HPP