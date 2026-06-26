#include "multigrid.hpp"

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
    char smoother)
{
    // std::cout << "\n=== V level " << l << " ===  norm_b = " << b[l].norm() 
    //           << "  norm_x_before = " << x[l].norm() << std::endl;

    //if i am on level 1 solve directyl
    if (l == 0) {
        x[0](0, 0) = b[0](0, 0) / A[0](0, 0);
        // std::cout << "  After coarse solve: norm_x = " << x[0].norm() << std::endl;
        return x;
    }
        //choose smoother 
        if(smoother=='g'||smoother=='G'){
            GaussSeidel(A[l],DL[l],x[l],b[l],nu1);
        }else if(smoother=='j'||smoother=='J'){
            DAMPEDJAC(A[l],DL[l],x[l],b[l],omega,nu1);
        }
// std::cout << "  After pre-smooth: norm_x = " << x[l].norm() << std::endl;
        

        //restriction 
       Matrix residual = b[l] - (A[l] * x[l]);
    //    std::cout << "  Residual norm = " << residual.norm() << std::endl;

        b[l-1] = R[l-1] * residual;
        // std::cout << "  After restriction: norm_b[" << l-1 << "] = " << b[l-1].norm() << std::endl;
        
        
        //initialize solution to lower level equal to zero
        x[l-1].fill(0.0);
        
        //recall function 
        x=V(A,DL,x,b,R,P,nu1,nu2,omega,l-1,smoother);
        
        //prologate 
        Matrix correction=P[l - 1] * x[l - 1];
        x[l] += correction;

        //smooth again 
        if(smoother=='g'||smoother=='G'){
            GaussSeidel(A[l],DL[l],x[l],b[l],nu2);
        }else if(smoother=='j'||smoother=='J'){
            DAMPEDJAC(A[l],DL[l],x[l],b[l],omega,nu2);
        }
    
    return x;
    }