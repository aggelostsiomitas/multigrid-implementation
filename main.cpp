#include<iostream>
#include <mkl.h>
#include<vector>


//sparse matrix struct
struct Matrix{
    int n;
    int m;
    std::vector<double>values;
    std::vector<int>col_indices;
    std::vector<int>row_ptr;


    Matrix& operator+=( const Matrix& other){
        for(int i=0;i<values.size();i++){
            this->values[i]+=other.values[i];
        }
        return *this;
    }

};

//zeros function
Matrix Zeros(int r,int c){
    Matrix A;
    A.n=r;
    A.m=c;

    A.values.assign(r,0.0);
    A.col_indices.assign(r,0);
    A.row_ptr.reserve(r+1);
    for (int i = 0; i <= r; ++i) {
        A.row_ptr.push_back(i);       
    }
    return A;
}

//tridiag function 
Matrix tridiag(int n ,double h){
    Matrix A;
    A.n=n;
    A.m=n;

    //number of non zero elements
   int nnz=(n>1)?(3*n-2):1;

   A.values.reserve(nnz);
   A.col_indices.reserve(nnz);
   A.row_ptr.reserve(n+1);

   double scalar=1.0/(h*h);

   A.row_ptr.push_back(0.0);

   for(int r=0;r<n;r++){
    //lower diagonal 
    if(r>0){
        A.values.push_back(-1.0*scalar);
        A.col_indices.push_back(r-1);
    }

    //main. diagonal 
    A.values.push_back(2.0*scalar);
    A.col_indices.push_back(r);

    //upper diagonal
    if(r<n-1){
       A.values.push_back(-1.0*scalar);
       A.col_indices.push_back(r+1);
    }
    A.row_ptr.push_back(A.values.size());
   }
   return A;

};

Matrix spalloc(int r,int c,int space){
    Matrix R;
    R.n=r;
    R.m=c;
    R.values.reserve(space);
    R.col_indices.reserve(space);
    R.row_ptr.reserve(r+1);

    return R;
};


Matrix diag(const Matrix& A,int n){
    Matrix D;
    D.n=A.n;
    D.m=1;

    //diagonal vector of A has n elements
    D.values.reserve(D.n);
    D.col_indices.reserve(D.n);
    D.row_ptr.reserve(D.n + 1);
    
    D.row_ptr.push_back(0);

    //search for diagonal elements
    for(int r=0;r<A.n;r++){
        bool found_diagonal=false;
        for(int idx=A.row_ptr[r];idx<A.row_ptr[r+1];idx++){
            //get elements A[i][i]
            if (A.col_indices[idx] == r){
                D.values.push_back(A.values[idx]);
                found_diagonal = true;
                break;
            }
        }
        //if i do not have diagonal put a zero
        if (!found_diagonal) {
            D.values.push_back(0.0);
        }
        D.col_indices.push_back(0);
        D.row_ptr.push_back(D.values.size());
    }
    return D;
};


void print(const Matrix&A){
    for(int r=0;r<A.n;r++){
        int current_idx = A.row_ptr[r];
        int row_end_idx = A.row_ptr[r + 1];
       for (int c = 0; c < A.m; c++){
         if (current_idx < row_end_idx && A.col_indices[current_idx] == c) {
                std::cout << A.values[current_idx] << "\t";
                current_idx++; // Προχωράμε στο επόμενο μη-μηδενικό στοιχείο της γραμμής
            } else {
                // Αν δεν υπάρχει αποθηκευμένο στοιχείο για αυτή τη στήλη, είναι μηδέν
                std::cout << 0.0 << "\t";
            }
       }
       std::cout << "\n";
    }
    std::cout << "\n";
}


double norm(const Matrix&A){
    double sum=0.0;
    for(auto val:A.values){
        sum+=val*val;
    }
    return std::sqrt(sum);

}

Matrix mult(const Matrix& A,const Matrix& x){
    Matrix Ax;
    Ax.n=A.n;
    Ax.m=1;

    Ax.values.reserve(Ax.n);
    Ax.col_indices.reserve(Ax.n);
    Ax.row_ptr.reserve(Ax.n+1);

    Ax.row_ptr.push_back(0);
    for(int r=0;r<A.n;r++){
        double row_sum=0.0;

        for(int idx=A.row_ptr[r];idx<A.row_ptr[r+1];idx++){
            int c=A.col_indices[idx];
            double val=A.values[idx];
            row_sum+=val*x.values[c];
        }
        Ax.values.push_back(row_sum);
        Ax.col_indices.push_back(0);
        Ax.row_ptr.push_back(Ax.values.size());
    }

    return Ax;
}


Matrix sub(const Matrix& x,const Matrix& y){
    Matrix result;
    result.n=x.n;
    result.m=1;

    int N=x.values.size();
    result.values.reserve(N);
    result.col_indices.reserve(N);
    result.row_ptr.reserve(x.n+1);

    result.row_ptr.push_back(0);

    for(int i=0;i<x.n;i++){
        double diff=x.values[i]-y.values[i];
        result.values.push_back(diff);
        result.col_indices.push_back(0);
        result.row_ptr.push_back(result.values.size());
    }

    return result;
}

void add( Matrix& x,const Matrix& y){
   
    for (int i = 0; i < x.n; i++) {
        x.values[i] += y.values[i]; 
    }
}


void divide(Matrix &x,const Matrix & divider){
    for(int i=0;i<x.n;i++){
        x.values[i]/=divider.values[i];
    }
}


void scalar_mult(Matrix& A,double scalar)
{
   for(int i=0;i<A.values.size();i++){
    A.values[i]*=scalar;
   }

}


Matrix Dampedjac(const Matrix&A,const Matrix&D, Matrix&x,const Matrix&b,double omega,int nu){

    for(int i=0;i<nu;i++){

        //A*x
        Matrix Ax=mult(A,x);

        //b-A*x
        Matrix bAx=sub(b,Ax);
        
        //omega*(b-A*x)
        scalar_mult(bAx,omega);
        
        //omega*(b-A*x)./D
        divide(bAx,D);
        //x=x+omega*(b-A*x)./D
        add(x,bAx);
    }

    return x;


}


Matrix TransposeAndScale(const Matrix&R,double scalar){
    Matrix P;

    //P(n,m)=R(m,n)
    P.n=R.m;
    P.m=R.n;


    return P;
};

//R^T*x->vector with column=1
Matrix mult_transpose_scaled(const Matrix& R , const Matrix &x){

    Matrix result;
    result.n=R.m;
    result.m=1;

    result.values.assign(result.n,0.0);
    result.col_indices.assign(result.n,0);
    result.row_ptr.reserve(result.n + 1);

    for (int i = 0; i <= result.n; ++i) {
        result.row_ptr.push_back(i);
    }

    for(int r=0;r<R.n;r++){
        double x_val=x.values[r];

        for(int idx=R.row_ptr[r];idx<R.row_ptr[r+1];idx++){
            int c=R.col_indices[idx];
            double R_val=R.values[idx];
            result.values[c]+=2.0*R_val*x_val;
        }
    }
    return result;
}



void V_cycle(
    const std::vector<Matrix>& A,
    const std::vector<Matrix>& D,
    std::vector<Matrix>& x,
    std::vector<Matrix>& b,
    const std::vector<Matrix>&R,
    int nu1,
    int nu2,
    double omega,
    int l)
{
    //if level=0 solve directly 
    if(l==0){
        x[0].values[0]=b[0].values[0]/A[0].values[0];
    }else {

        //smooth 
        x[l]=Dampedjac(A[l],D[l],x[l],b[l],omega,nu1);

        //compute residual 
        Matrix Ax=mult(A[l],x[l]);
        Matrix bAx=sub(b[l],Ax);

        //restrict 
        b[l-1]=mult(R[l-1],bAx);

        //initialize solution to lower level=0
        x[l-1]=Zeros(b[l-1].n,1);
        
        //recall function  for lower level 
        V_cycle(A, D, x, b,R, nu1, nu2, omega, l - 1);

        
        //perform prologation
        Matrix correction=mult_transpose_scaled(R[l-1],x[l-1]);
        add(x[l],correction);

        //smooth again
        x[l] = Dampedjac(A[l], D[l], x[l], b[l], omega, nu2);
    }
    
}
int main(){
    

    //parameters 
    int d=1;
    double omega = 2.0 * d / (2.0 * d + 1.0);
    int nu1=1;
    int nu2=2;
    int Nmax=50;
    double tol =1e-8;

    //number of levels
    int levels=5;

    //uknowns per level
    std::vector<int>n(levels);

    for(int i=0;i<levels;i++){
        n[i]=(1<<(i+1))-1;
    }
    
    //mesh sizes
    std::vector<double>h(levels);

    for(int i=0;i<levels;i++){
        h[i]=1.0/(n[i]+1);
    }
    

    //coefficient matrices 
    std::vector<Matrix>A(levels);
    std::vector<Matrix>D(levels);
    std::vector<Matrix>b(levels);
    std::vector<Matrix>x(levels);
    std::vector<Matrix>R(levels-1);

    for(int i=0;i<levels;i++){
        A[i]=tridiag(n[i],h[i]);
        x[i]=Zeros(n[i],1);
        b[i]=Zeros(n[i],1);
        D[i]=diag(A[i],A[i].n);
    }
    

    //create RHS for finest level 
    for(int i=0;i<b[levels-1].n;i++){
        b[levels-1].values[i]=1.0;
    }
   
    //construct the Restriction matrix
   for(int i=0;i<levels-1;i++){
    //reserve space 
    R[i]=spalloc(n[i],n[i+1],3*n[i]);

    R[i].row_ptr.push_back(0);

    for(int j=0;j<n[i];j++){

        //indexes [2j-1:2j+1]=[0.25,-0.5,0.25]
        double val[3]={0.25,0.5,0.25};
        int index[3]={2*j,2*j+1,2*j+2};
        for(int k=0;k<3;k++){
            if(index[k]>=0 &&index[k]<n[i+1]){
            R[i].values.push_back(val[k]);
            R[i].col_indices.push_back(index[k]);
            } 
        }
        R[i].row_ptr.push_back(R[i].values.size());
    }
}

// print(R[1]);
// print(R[1]);
// double normb=norm(R[1]);
// std::cout<<normb;

double normb=norm(b[levels-1]);
for(int i=0;i<Nmax;i++){
     V_cycle(A, D, x, b,R, nu1, nu2, omega, levels-1);
     
     Matrix Ax=mult(A[levels-1],x[levels-1]);
     Matrix bAx=sub(b[levels-1],Ax);
     double normr=norm(bAx);
     std::cout<< "iteariotn: "<<i<< "with norm "<<normr<<'\n';
     if(normr<tol*normb){
        std::cout<<"converged after "<<i <<" iterations\n";
        break;
     }
     
}
    return 0;
}