#ifndef MULTIGRID_HPP
#define MULTIGRID_HPP

#include "matrix.hpp"
#include "dampedJac.hpp"
#include "GaussSeidel.hpp"

/**
 * @brief V cycle Multigrid
 * Performs Multigrid iteration
 * @param A  System matrix A
 * @param DL Lower triangular part of A.
 * @param x  Solution vector (updated in place).
 * @param b  Right-hand side.
 * @param R Resitriction matrix 
 * @param P Prologation matrix 
 * @param nu1 Number of iterations for resitriction
 * @param nu2 number of iterations for proloagation
 * @param omega omega variable for Damped Jacobi method
 * @param l The number of levels of the method
 * @param smoother A character for the selection of the smoother
 *//** @file Multigrid.hpp */

 std::vector<Matrix> V(
    const std::vector<Matrix>&A,
    const std::vector<Matrix>&DL,
    std::vector<Matrix>&x,
    std::vector<Matrix>&b,
    const std::vector<Matrix>&R,
    const std::vector<Matrix>&P,
    int nu1,
    int nu2,
    double omega ,
    int l,
    char smoother
);

#endif //MULTIGRID_HPP