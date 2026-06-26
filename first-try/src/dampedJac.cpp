#include "dampedJac.hpp"

Matrix  DAMPEDJAC(
    const Matrix& A,
    const Matrix& D_mat,
    Matrix& x ,
    const Matrix& b,
    double omega,
    int nu)
{

    Matrix D=D_mat.getDiagonal();

    for(int iter=0;iter<nu;iter++){
    //   x+=omega*(b-(A*x)).Div(D);
        Matrix Ax = A * x;
        Matrix residual = b - Ax;
        Matrix correction = residual.Div(D);
        x += omega * correction;
    }
    return x;
}