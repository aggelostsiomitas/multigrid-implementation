#include<iostream>
#include "matrix.hpp"
#include "dampedJac.hpp"
#include "GaussSeidel.hpp"
#include "multigrid.hpp"
#include<iostream>
#include<iomanip>
int main(){
    //parameters 
    int d=1;
    int nu1=1;
    int nu2=2;
    double omega=(2.0*d)/(2.0*d+1.0);
    int Nmax=100;
    double tol=1e-8;
    int levels=5;
    char smoother='g';


    //necessary Matrices
    std::vector<int> n(levels);
    std::vector<Matrix> A(levels);
    std::vector<Matrix> DL(levels);
    std::vector<Matrix> x(levels);
    std::vector<Matrix> b(levels);
    std::vector<Matrix> R(levels);
    std::vector<Matrix> P(levels);
    
    for(int i = 0; i < levels; i++){
        n[i] = (1 << (i + 1)) - 1;
    }
    
    for(int i=0;i<levels;i++){
        double h=1.0/(n[i]+1);

        //x(i)=zeros(n(i),1)
        x[i]=Matrix(n[i],1,0.0);
        //b(i)=zeros(n(i),1)
        b[i]=Matrix(n[i],1,0.0);

        //construct matrix A tridiag
        A[i]=Matrix(n[i],n[i]);
        for(int r=0;r<n[i];r++){
            A[i](r,r)=2.0/(h*h);
            if (r>0){ A[i](r,r-1)=-1.0/(h*h);}
            if (r<n[i]-1){ A[i](r,r+1)=-1.0/(h*h);}
        }
        
        //construct L matrix 
        DL[i]=Matrix(n[i],n[i]);
        for(int r=0;r<n[i];r++){
            DL[i](r, r) = A[i](r, r);
            if (r > 0) DL[i](r, r - 1) = A[i](r, r - 1);
        }
        //construct restriction and prologation matrices 
        if(i<levels-1){
            
            R[i]=Matrix(n[i],n[i+1]);
            std::cout << "R[" << i << "] created: " << n[i] << " x " << n[i+1] << std::endl;
            for(int j=0;j<n[i];j++){
                int idx=2*j+1;
                R[i](j,idx-1)=0.25;
                R[i](j,idx)=0.5;
                R[i](j,idx+1)=0.25;
            }
            P[i]=2.0*R[i].Transpose();
        }else {
            // Optional: initialize last level to something safe
            R[i] = Matrix(n[i], n[i]);   // or leave default
            P[i] = Matrix(n[i], n[i]);
        }
    }
    //create rhs for the finest grid
    b[levels - 1].fill(1.0);
    
    //compute norm of highest lvel 
    double normb=b[levels-1].norm();

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Starting CSR-Sparse Multigrid V-Cycle...\n";
    std::cout << "Finest Level Matrix NNZ Param Value: " << A[levels - 1].getNNZ() << "\n\n";
    for(int iter=0;iter<Nmax;iter++){
        x=V(A,DL,x,b,R,P,nu1,nu2,omega,levels-1,smoother);
        Matrix res=b[levels-1]-(A[levels-1]*x[levels-1]);
        double normr=res.norm();
        
        if(normr<tol*normb){
            std::cout << "\nConverged after iterations: " << iter << '\n';
            break;
        }
    }
    return 0;
}
