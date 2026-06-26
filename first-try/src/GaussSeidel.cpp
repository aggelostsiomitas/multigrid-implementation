#include "GaussSeidel.hpp"

void GaussSeidel(
    const Matrix& A,
    const Matrix& DL,
    Matrix& x,
    const Matrix& b,
    int nu)
{
    for(int iter=0;iter<nu;iter++){
        Matrix BAX=b-(A*x);
        x+=DL.SolveLower(BAX);
    }
}