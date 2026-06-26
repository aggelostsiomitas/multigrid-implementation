#ifndef MATRIX_HPP
#define MATRIX_HPP

#include<iostream>
#include<cstddef>
#include<vector>
#include<stdint.h>
#include<cmath>
#include<stdexcept>
/**
 * @brief Matrix class
 * \param rows the Number of rows of my matrix 
 * \param cols the Number of columns of my matrix 
 * \param values a vector that will store the values of my matrix
 * \param col_indices a vector that will store the column indices of the values in value vector
 * \param row_indices a vector that will store the row indices of the values in value vector
 * \param nnz a variable that will store the number of non zero elements
 */

 class Matrix{
private:
    size_t rows;
    size_t cols;
    std::vector<double>values;
    std::vector<int>col_indices;
    std::vector<size_t>row_indices;
    unsigned int nnz;
public:
    /**
     * \brief Default Constructor
     * The default Constructor of the class 
     */
    Matrix();
    /**
     * \brief Primary Constructor
     * The primary Constructor of the class 
     * \param r the number of rows of my matrix 
     * \param c the number of columns of my matrix 
     */
    Matrix(size_t r,size_t c);
    /**
     * \brief Second Constructor
     * This Constructor initialize a matrix with an initial value 
     * \param r the number of rows of my matrix 
     * \param c the number of columns of my matrix 
     * \param initial_value the initial value if my matrix values
     */
    Matrix(size_t r,size_t c,double initial_value);

    /**
     * \brief Helper getter functions to return the data of the matrix class
     */
    size_t getRows()const; 
    size_t getCols()const;
    bool isVector()const;
    unsigned int getNNZ()const;
    double norm() const;

    /**
     * \brief Fill function to assign values the the values vector
     * \param value the value i will assign
     */
    void fill(double value);

    /**
     * \brief Extract diagonal element of Matrix
     */
    Matrix getDiagonal()const;

    /**
     * \brief Divide with diagonal element of Matrix
     * \param other The Matrix with who i will perform the division
     */
    Matrix Div(const Matrix& other)const;

    /**
     * \brief Lower Triangualr solve
     * \param rhs The right hand side Matrix
     */
    Matrix SolveLower(const Matrix& rhs)const;


    /**
     * \brief Operator overloading
     * The necessary operator overloading for my Matrix class
     */
    double operator()(size_t r,size_t c)const;
    double& operator()(size_t r,size_t c);
    double& operator()(size_t idx);
    const double& operator()(size_t idx)const;

    /**
      * \brief multiplication operator overloading 
      * \param scalar The scalar va;ue with which i will multiply the values in my Matrix
      */
    Matrix operator*(double scalar) const;

     /**
      * \brief multiplication operator overloading 
      * \param scalar The scalar va;ue with which i will multiply the values in my Matrix
      */
     Matrix& operator*=(double scalar);

     /**
      * \brief in place addition operator overloading 
      * \param other The Matrix with who i will perform the division
      */
     Matrix& operator+=(const Matrix& other);
     /**
      * \brief in place substraction operator overloading
      * \param other The Matrix with who i will perform the division 
      */

     Matrix& operator-=(const Matrix& other);
     /**
      * \brief multiplication operator overloading 
      * \param other The Matrix with who i will perform the division
      */
     Matrix operator*(const Matrix& other)const;

     /**
      * \brief Transpose function
      * Compute the transpose equivalent of my matrix
      */
     Matrix Transpose()const;
     void print(const std::string&name="")const;
 }; 
  /**
      * \brief addition operator overloading 
      */
    inline Matrix operator+(Matrix a, const Matrix& b) { a += b; return a; }     
     /**
      * \brief substraction operator overloading 
      */
    inline Matrix operator-(Matrix a, const Matrix& b) { a -= b; return a; }
    
     /**
      * \brief multiplication operator overloading 
      */
    inline Matrix operator*(double scalar, Matrix m) { m *= scalar; return m; }

#endif 