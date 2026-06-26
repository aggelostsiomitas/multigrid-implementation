#include "matrix.hpp"

//default constructor
Matrix::Matrix():rows(0),cols(0),nnz(0){row_indices.push_back(0);}

//primary constructor
Matrix::Matrix(size_t r,size_t c):rows(r),cols(c),nnz(0){
    row_indices.assign(r+1,0);
}

//second constructor with initial value 
Matrix::Matrix(size_t r,size_t c ,double initial_value):rows(r),cols(c),nnz(0){
    
    if(c==1 && initial_value!=0){
        values.assign(r,initial_value);
        col_indices.assign(r,0);
        row_indices.resize(r + 1);
        for(size_t i=0;i<=r;i++){
            row_indices[i]=i;
        }
        nnz=static_cast<unsigned int>(r);
    }
    else{
        row_indices.assign(r+1,0);
        nnz=0;
    }
}

//get rows function
size_t Matrix::getRows()const{
    return rows;
}

//get cols function
size_t Matrix::getCols()const{
    return cols;
}
//check if the Matrix is a column
bool Matrix::isVector()const{
    return cols==1;
}
//get non zero elements function
unsigned int  Matrix::getNNZ()const{
    return nnz;
}

//norm function
double Matrix::norm()const {
    double sum=0.0;
    for(double val:values){
        sum+=val*val;
    }
    return std::sqrt(sum);
}

void Matrix::fill(double value ){
    if (value==0){
        values.clear();
        col_indices.clear();
        std::fill(row_indices.begin(),row_indices.end(),0);
        nnz=0;
    }else{
        if(isVector()){
            values.assign(rows,value);
            col_indices.assign(rows,0);
            for(size_t i=0;i<=rows;i++){
                row_indices[i]=i;
            }
            nnz=static_cast<unsigned int>(rows);
        }else{
            std::fill(values.begin(), values.end(), value);
        }
    }
}

//operator () overloading 
double Matrix::operator()(size_t r,size_t c)const{
    if (r >= rows || c >= cols) return 0.0;
    size_t  start=row_indices[r];
    size_t  end=row_indices[r+1];
    for(size_t  i=start;i<end;i++){
        if(col_indices[i]==static_cast<int>(c)){return values[i];}
    }
    return 0.0;
}

double& Matrix::operator()(size_t r,size_t c){

   if (r >= rows || c >= cols) {
        std::cerr << "OUT OF RANGE: r=" << r << " c=" << c 
                  << " rows=" << rows << " cols=" << cols << std::endl;
        throw std::out_of_range("Matrix::operator()(r,c): index out of range");
    }

    //search for the element
   size_t  start=row_indices[r];
   size_t  end=row_indices[r+1];

   for(size_t  i=start;i<end;++i){
    if(col_indices[i]==static_cast<int>(c)){  return values[i];}
   }
 

   //if i have not found the element insert it 
   size_t  insert_pos = start;
   //initialize the searching for the correct position for the element to be added (sorted)
   while (insert_pos < end && col_indices[insert_pos] < static_cast<int>(c)) {
            ++insert_pos;
        }
    // add the new element
    values.insert(values.begin() + insert_pos, 0.0);
        col_indices.insert(col_indices.begin() + insert_pos, static_cast<int>(c));

        //fix the row_indices
        for (size_t i = r + 1; i <= rows; ++i) {
            row_indices[i]++;
        }
        
        //increase the nnz since i have a new element
        nnz++; 
        return values[insert_pos];

}


double& Matrix::operator()(size_t idx){
    if (!isVector())
        throw std::runtime_error("Matrix::operator()(idx) called on non-vector");
    
    return (*this)(idx, 0); 
}
const double& Matrix::operator()(size_t idx)const{
     if (!isVector())
        throw std::runtime_error("Matrix::operator()(idx) called on non-vector");
    if (idx >= values.size())
        throw std::out_of_range("Matrix::operator()(idx): index out of range");
    return values[idx];
}

Matrix& Matrix::operator*=(double scalar){
    for (auto& val:values){
        val*=scalar;
    }
    return *this;
}

Matrix& Matrix::operator+=(const Matrix&other){
    if(rows!= other.rows || cols!=other.cols){
        throw std::invalid_argument("Matrix dimensions must agree");
    }
    if (isVector() && other.isVector()) {
        for (size_t i = 0; i < rows; ++i) {
            (*this)(i, 0) += other(i, 0);   
        }  
    }else{
        for (size_t i = 0; i < values.size(); ++i) {
            values[i] += other.values[i];
        }
    }
    return *this;
}   

//+= operator overload 
Matrix& Matrix::operator-=(const Matrix&other){
    if(rows!= other.rows || cols!=other.cols){
        throw std::invalid_argument("Matrix dimensions must agree");
    }
    for(size_t i=0;i<values.size();i++){
        values[i]-=other.values[i];
    }
    return *this;
}  

Matrix Matrix::operator*(double scalar) const {
    Matrix result = *this;  // copy
    result *= scalar;
    return result;
}

//matrix multiplication
Matrix Matrix::operator*(const Matrix& other)const{
    //if the other is a Vector 
    if(other.isVector() && cols==other.rows){
        Matrix result(rows ,1,0.0);
        for(size_t i=0;i<rows;i++){
            double sum=0.0;
            int start=row_indices[i];
            int end=row_indices[i+1];
            for(int k=start;k<end;k++){
                sum+=values[k]*other(col_indices[k],0);
            }
            result(i,0)=sum;
        }
        return result;
    }
    // General Matrix-Matrix multiplication
    Matrix result(rows, other.cols);

    for (size_t i = 0; i < rows; ++i) {
        size_t start = row_indices[i];
        size_t end = row_indices[i + 1];

        for (size_t k = start; k < end; ++k) {
            int c = col_indices[k];
            double val = values[k];

            // Multiply this row by the c-th row of 'other'
            size_t other_start = other.row_indices[c];
            size_t other_end = other.row_indices[c + 1];

            for (size_t j = other_start; j < other_end; ++j) {
                int other_col = other.col_indices[j];
                result(i, other_col) += val * other.values[j];
            }
        }
    }
    return result;
}


//Transpose function
Matrix Matrix::Transpose()const{
   
    Matrix result(cols,rows);
    if (nnz == 0) return result;
    result.nnz=nnz;
    
    std::vector<size_t > count(cols,0);
    // for(size_t i=0;i<values.size();i++){count[col_indices[i]]++;}
    for (int c : col_indices) {
        count[c]++;
    }
    result.row_indices.resize(cols+1,0);
    for(size_t i=0;i<cols;i++){
        result.row_indices[i+1]=result.row_indices[i]+count[i];
    }
    
    // result.values.resize(values.size());
    // result.col_indices.resize(col_indices.size());
    result.values.resize(nnz);
    result.col_indices.resize(nnz);

    std::vector<size_t> current_pos = result.row_indices;
    for(size_t r=0;r<rows;r++){
        size_t  start=row_indices[r];
        size_t  end=row_indices[r+1];
        for(size_t  k=start;k<end;++k){
            int  c=col_indices[k];
            size_t  dest=current_pos[c]++;
            result.values[dest] = values[k];
            result.col_indices[dest] = static_cast<int>(r);
        }
    }
   return result;
}

// extract D matrix
Matrix Matrix::getDiagonal()const{
    Matrix diagVec(rows,1,0.0);
    for(size_t i=0;i<rows;i++){
        diagVec(i,0)=(*this)(i, i);
    }
    return diagVec;
}

Matrix Matrix::Div(const Matrix& other)const {
    if (!isVector() || !other.isVector() || rows != other.rows)
        throw std::runtime_error("Div requires two vectors of same size");

    Matrix result(rows, 1, 0.0);

    for (size_t i = 0; i < rows; ++i) {
        double val   = (*this)(i, 0);
        double denom = other(i, 0);
        if (std::abs(denom) < 1e-14) {
            result(i, 0) = 0.0;
        } else {
            result(i, 0) = val / denom;
        }
    }
    return result;
}

Matrix Matrix::SolveLower(const Matrix& rhs)const{
    Matrix y(rows,1,0.0);

    for (size_t i=0;i<rows;i++){
        int start=row_indices[i];
        int end=row_indices[i+1];
        double sum=0.0;
        double diag_val=0.0;

        for (size_t k=start;k<end;k++){
            int c=col_indices[k];
            if(c<static_cast<int>(i)){
                sum+=values[k]*y(c,0);
            }else if(c == static_cast<int>(i)){
                diag_val=values[k];
            }
        }
        y(i,0)=(rhs(i,0)-sum)/diag_val;
    }
    return y;
}


void Matrix::print(const std::string&name )const{
    if (!name.empty()) {
        std::cout << name << " (" << rows << " x " << cols << "):" << std::endl;
    } else {
        std::cout << "Matrix (" << rows << " x " << cols << "):" << std::endl;
    }

    for (size_t r = 0; r < rows; ++r) {
        std::cout << "Row " << r << ": ";
        size_t start = row_indices[r];
        size_t end = row_indices[r + 1];

        if (start == end) {
            std::cout << " [all zeros]" << std::endl;
            continue;
        }

        for (size_t k = start; k < end; ++k) {
            std::cout << "(col=" << col_indices[k] << ", val=" << values[k] << ") ";
        }
        std::cout << std::endl;
    }
    std::cout << "NNZ = " << nnz << std::endl << std::endl;
}

