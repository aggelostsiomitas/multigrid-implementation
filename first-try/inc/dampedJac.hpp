#ifndef DAMPEDJAC_HPP
#define DAMPEDJAC_HPP

#include<vector>
#include "matrix.hpp"

/**
 * @brief Damped Jacobi smoother.
 *
 * Implements the damped Jacobi iteration for A x = b.
 *
 * @param A      System matrix.
 * @param D_mat  Diagonal of A.
 * @param x      Solution vector (updated in place).
 * @param b      Right-hand side.
 * @param omega  Relaxation parameter.
 * @param nu     Number of iterations.
 */

Matrix  DAMPEDJAC(
    const Matrix& A,
    const Matrix& D_mat,
    Matrix& x ,
    const Matrix& b,
    double omega,
    int nu);

#endif // DAMPEDJAC_HPP
