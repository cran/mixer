//  matrix.h
// 
//  Copyright (C) 2008 Laboratoire Statistique & Genome
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or (at
//  your option) any later version.
// 
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
//

//! \file    matrix.h
//! \author  Gilles Grasseau
//! \brief   Basic tools for matrix
//! \version 0.1
//! \date    November 2008

# include <iostream>
# include <cstring>
# include <cmath>
# include <cfloat>
# include <string.h>
# include <math.h>
# include "util.h"

# ifndef _MATRIX_H_
# define _MATRIX_H_

static const char *not_same_dims   = "Matrix don't have the same dimensions";
static const char *not_dims_match  = "Matrix dimensions don't match";
static const char *not_scalar      = "Matrix can't be trasformed to scalar";
static const char *out_of_dims     = "Index out of matrix dimensions";
static const char *not_square_matrix = "Not a squared Matrix";
static const char *not_implemented = "Not yet implemented";

//! \cond
#define OPT
#define EQUAL_DIM( X ) ( (row_ == (X).row_) && (col_ == (X).col_) ) 

#define EQUAL_DIM_MAT( X, Y ) ( ((X).getRow() == (Y).getRow()) \
                             && ((X).getCol() == (Y).getCol()) )

#define EQUAL_COL_ROW( X ) ( (col_ == (X).row_) )  
//! \endcond

//! \class Matrix matrix.h
//! \brief Class used to make matrix computations
//!  The matrix elements are row-major ordered i.e. A(i,k) and A(i+1,k) are 
//!  contiguous in the memory. 
class Matrix {

  private :

  size_t row_; 
  size_t col_;

  struct s_data {
    // Reference number on matrix array matrix_
    int ref_;
    // Row ordered first
    double *matrix_;
  } *link_;
 
  struct s_data *newLink(size_t size_ ){

    struct s_data *link = new struct s_data;
    link -> matrix_ = new double [size_];
    link -> ref_    = 1;



    return(link);
  }

  struct s_data *newLink( double *ptr_ ){

    struct s_data *link = new struct s_data;
    link -> matrix_ = ptr_;
    link -> ref_    = 1;

    return(link);
  }

  void deleteLink( ){

    if ( (-- (link_-> ref_)) == 0 ) {

      delete [] (link_->matrix_);
      delete link_;
      link_ = 0;
    }

  }

  void fillData( double a ) {     
    for (size_t i=0; i < row_*col_; i++) {
      link_->matrix_[i] = a;
    }
  } 

  void fillData( const double *ptr ) {   

    memcpy( link_->matrix_, ptr, row_*col_*sizeof(double) );
  } 

  public :

    //! \brief Constructor - The matrix is not created.
    Matrix() {
      row_ = 0; col_ = 0;
      // To initialize link_
      link_ = newLink(1);
    }

    //! \brief Constructor - The matrix A(row,col) contains indeterminated
    //!  values.
    //! \param row Number of rows.
    //! \param col Number of lines.
    Matrix( const size_t row, const size_t col) {
      row_ = row; col_ = col;
      link_ = newLink( row_*col_ );
    }

    //! \brief Constructor - The matrix  A(row,col) is filled with a value.
    //! \param row Number of rows.
    //! \param col Number of lines.
    //! \param a   Value to fill the matrix with.
    Matrix( const size_t row, const size_t col, double a) {
      row_ = row; col_ = col;
      link_ = newLink( row_*col_ );
      fillData( a );
    }

    //! \brief Constructor -  The matrix  A(row,col) is filled with
    //!        the "ptr" pointer contain.
    //! \param row  Number of rows.
    //! \param col  Number of lines.
    //! \param ptr  Floating point array.
    //! \param copy Copy the  Floating point "ptr" array if true, 
    //!             else "ptr" is directly assign to the internal matrix data.
    //!             In "false" case do not use the pointer for other 
    //!              computations  
    Matrix( const size_t row, const size_t col, double *ptr, 
	    const bool copy=true ) {

      row_ = row; col_ = col;
      if ( copy ) {
	link_ = newLink(row_*col_) ;
	fillData( ptr );
      } else {
	link_ = newLink( ptr );
      }

    }
    //! \brief Constructor -  The returned matrix (rowcol,col) or 
    //!        (row, rowcol) is filled with the vector V
    //!        (matrix(1,col) or matrix(row,1).
    //! \param rowcol  Number of rows or columns.
    //! \param V       Vector (Matrix (1,col) or (row, 1)). 
    //! \param tranpose Tranpose parameter (default false), 
    //!             else "ptr" is directly assign to the internal matrix data.
    //!             else "ptr" is directly assign to the internal matrix data.
    //!             In "false" case do not use the pointer for other 
    //!              computations  
    Matrix( const size_t rowcol, const Matrix &V, 
	    const char *transpose=false ) {
      double *tmp_tab = 0;



      const double *tabV = V.getData();

      if( V.row_ == 1 ) {
	row_ = rowcol;
        col_ = V.col_;
	tmp_tab = new double [ row_ * col_ ];
	for (size_t j=0; j < col_; j++) {
	  for (size_t i=0; i < row_; i++) {
	    tmp_tab[j*row_+i] = tabV[j];
	  }
	}
	link_ = newLink( tmp_tab );
      } else if( V.col_ == 1 ) {
	row_ = V.row_;
        col_ = rowcol;
	tmp_tab = new double [ row_ * col_ ];
	for (size_t j=0; j < col_; j++) {
	  for (size_t i=0; i < row_; i++) {
	    tmp_tab[j*row_+i] = tabV[i];
	  }
	}
	link_ = newLink( tmp_tab );
      }

    }


    //! Copy constructor
    //! \param A The new matrix shared the same data with matrix A.
    //!        The new matrix is an alias of A.
    //!        Do not used this constuctor to duplicate A. 
    //! \code
    //!  B = Matrix(A);        // Forbiden to copy matrix
    //!  Matrix C = A;         // Forbiden to copy matrix
    //! \endcode
    //! See getCopy() to duplicate Matrix objects.
    Matrix( const Matrix& A ) {
      row_ = A.row_;
      col_ = A.col_;
      A.link_-> ref_++;
      link_ = A.link_;
    }

    //! Destructor
    ~Matrix() {

      deleteLink();

    }

    //! Assignement operator
    Matrix& operator = ( const Matrix& A) {

      row_ = A.row_;
      col_ = A.col_;
      deleteLink();

      link_ = A.link_;
      link_ -> ref_++;


      return(*this);
    }

    //! \cond  Symetric operator
    friend Matrix inline operator+(const double a, const Matrix& A );
    friend Matrix inline operator-(const double a, const Matrix& A );
    friend Matrix inline operator*(const double a, const Matrix& A );
    //! \endcond

    //! \brief Compute the maximum between the A matrix elements and the scalar 
    //! value a.
    //! \param A Matrix.
    //! \param a Scalar.
    //! \return Result Matrix.
    friend inline Matrix max( const Matrix& A, const double  a  );
    //! \cond  Symetric operator
    friend inline Matrix max( const double  a, const Matrix& A );
    //! \endcond

    //! \brief Compute the maximum element-by-element between the A matrix 
    //! and the B matrix elements. 
    //! \param A Matrix.
    //! \param B Matrix.
    //! \return Result Matrix.
    friend inline Matrix max( const Matrix& A, const Matrix& B );

    //! \brief Compute the minimum between the A matrix elements and the scalar 
    //! value a.
    friend inline Matrix min( const Matrix& A, const double  a );
    //! \cond  Symetric operator
    friend inline Matrix min( const double  a, const Matrix& A );
    //! \endcond

    //! \brief Compute the minimun element-by-element between the A matrix 
    //! and the B matrix elements. 
    //! \param A Matrix.
    //! \param B Matrix.
    //! \return Result Matrix.
    friend inline Matrix min( const Matrix& A, const Matrix& B );

    //! \brief Compute the sum of all A matrix elements. 
    //!        sum = sum( A(i,j), 0 < i < row, 0 < j < col) 
    //! \param A Matrix.
    //! \return The sum of all elements.

    //! \brief Compute the sum by row or by column of 
    //!        a matrix or the sum of all matrix elements. 
    //         Return a vector or scalar :
    //!          if dim = "col" (first dimension) 
    //!             v(1,k) = sum( A(i, k), 0 < i < row ).
    //!          if dim = "row" (second dimension) 
    //!             v(k,1) = sum( A(k, j), 0 < j < col ).
    //!          if dim = "all"
    //!            sum = sum( A(i,j), 0 < i < row, 0 < j < col) 
    //! \param A    Matrix.
    //! \param dim  Specify on which the summation is computed
    //!               dim="row", the sum is computed on rows 
    //!               dim="col", the sum is computed on columns 
    //!               dim="all", the sum is computed on all matrix element 
    //! 
    //! \return Return a sum vector. 
    friend inline Matrix sum( const Matrix& A, const char *dim="all" );

    //! \brief Compute the log() of A matrix elements. 
    //!        X(i,j) = log( A(i,j) ); 0 <= i < row, 0 <= j < col 
    //! \param A Matrix.
    //! \return log(A) Matrix.
    friend inline Matrix log( const Matrix& A );

    //! \brief Compute the exp() of A matrix elements. 
    //!        X(i,j) = exp( A(i,j) ); 0 <= i < row, 0 <= j < col 
    //! \param A Matrix.
    //! \return exp(A) Matrix.
    friend inline Matrix exp( const Matrix& A );

    //! \brief Compute the abs() of A matrix elements. 
    //!        X(i,j) = abs( A(i,j) ); 0 <= i < row, 0 <= j < col 
    //! \param A Matrix.
    //! \return abs(A) Matrix.
    friend Matrix inline abs( const Matrix& A );

    // Operators 

    //! \brief Type cast operator from a matrix (1,1) to a scalar.
    inline operator double(void) {
      double a=0;

      if ((row_ == 1) && (col_ == 1)) {
	a = link_->matrix_[0];
      }

      return( a );
    }

    //! \brief Type cast operator from a Matrix to const *Matrix.
    //inline operator const &Matrix(void) {
    //  return( a );
    //}

    //! \brief Indexing of the matrix elements
    //! \param i Row indice.
    //! \param j Column indice.
    //! \return Return A(i,j) value.
    double inline &operator()(size_t i, size_t j) { 
      return ( link_->matrix_[j*row_+i] ); 
    };

    //! \brief Linear indexing of the matrix array
    //! Remark: The matrix elements are row-major ordered.
    //! \param i Indice, 0 <= i < row*col .
    //! \return Return A(i) scalar.
    double inline &operator()(size_t i) { 
      return ( link_->matrix_[i] ); 
    };

    //!cond
    //! Functionnal operators for "const Matrix"
    double inline operator()(size_t i, size_t j) const { // For "const Matrix" 
      return ( link_->matrix_[j*row_+i] ); 
    };    
    double inline operator()(size_t i) const { // For "const Matrix" 
      return ( link_->matrix_[i] ); 
    };    
    // \endcond
    
    // \cond
    //! Unitary operator
    Matrix inline operator+( );
    Matrix inline operator-( );
    // \endcond

    //! \brief Addition operator of a Matrix with a scalar.
    //! All Matrix elements are added with "a" value.
    //! \param a Scalar to add.
    //! \return Matrix 
    Matrix inline operator+(const double a );
    //! \brief Subtraction operator of a Matrix with a scalar.
    //! All Matrix elements are substracted with "a" value.
    //! \param a Scalar to substract.
    //! \return Matrix 
    Matrix inline operator-(const double a );
    //! \brief Multiply operator of a Matrix with a scalar.
    //! All Matrix elements are multiplied with "a" value.
    //! \param a Scalar to multiply.
    //! \return Result Matrix 
    Matrix inline operator*(const double a );
    //! \brief Divide operator of a Matrix with a scalar.
    //! All Matrix elements are divided with "a" value.
    //! \param a Scalar to divide.
    //! \return Result Matrix 
    Matrix inline operator/(const double a );

    //! \brief Add element-by-element the two matrices.
    //! \param A Matrix to add.
    //! \return Result Matrix  
    Matrix inline operator+(const Matrix& A);
    //! \brief Substract element-by-element the two matrices.
    //! \param A Matrix to substract.
    //! \return Result Matrix  
    Matrix inline operator-(const Matrix& A);
    //! \brief Mutiply element-by-element the two matrices.
    //! \param A Matrix to mutiply.
    //! \return Result Matrix  
    Matrix inline operator|(const Matrix& A);
    //! \brief Matrices multiplication.
    //! \param A Matrix to multiply.
    //! \return Result Matrix  
    Matrix inline operator*(const Matrix& A);

    //! \brief Duplicate the Matrix.
    //! Perform a deep copy of the object
    //! \return Return a copy of the object.  
    Matrix getCopy();

    //! \brief Get the number of rows in Matrix.
    //! \return Return the row number.  
    size_t inline getRow( ) const { return ( row_ ); }; 
    //! \brief Get the number of rows in Matrix.
    //! \return Return the row number.  
    size_t inline getCol( ) const { return ( col_ ); }; 

    //! \brief Get the matrix array.
    //! The content of the pointer cannot be modified. 
    //! \return Return a pointer to the data array.  
    const double *getData() const { return( link_->matrix_ ) ; }

    //! \brief Copy the matrix array in the "ptr_dest" content.
    //! \param Pointer to an array of doubles. 
    //!        The ptr_dest MUST be allocated by the user.  
    void copyData(double *ptr_dest) const { 
      double *mat_ = link_->matrix_;

      for (size_t i=0; i < row_*col_; i++) {
	ptr_dest[i] = mat_[i];
      }
     }

    //! \brief Tranpose the matrix. 
    //! \param A Matrix to multiply.
    //! \return Return a new Matrix transposed.  
    Matrix inline t(void);

    //! \brief Get the row vector i.
    //! \param i Row index.
    //! \return Return the vector as Matrix object.  
    Matrix inline getRowVector( const size_t i ); 
    //! \brief Get the column vector j.
    //! The content of the pointer cannot be modified. 
    //! \param i Column index.
    //! \return Return the vector as Matrix object.  
    Matrix inline getColVector( const size_t j ); 
    //! \brief Get the upper triangular matrix.
    //! Return a Matrix (row, col):
    //! - if (diag = false):
    //!   - A(i,j), 0 <= i < row, i < j < col
    //!   - zero for other matrix elements.
    //! - else (diag = true):
    //!   - A(i,j), 0 <= i < row, i <= j < col
    //!   - zero for other matrix elements.
    //! \param diag Specify if the diagonal is included
    //!        (default diag=false).
    //! \return Return the vector as Matrix object.  
    Matrix inline getUpperTriang( const bool diag=false );
    //! \brief Get the lower triangular matrix.
    //! Return a Matrix (row, col)
    //! if (diag = false):
    //!  - A(i,j), 0 <= i < row, 0 <= j < i
    //!  - zero for other matrix elements.
    //! else (diag = true):
    //!  - A(i,j), 0 <= i < row, 0 <= j <= i
    //!  - zero for other matrix elements.
    //! \param diag Specify if the diagonal is included
    //!        (default diag=false).
    //! \return Return the vector as Matrix object.  
    Matrix inline getLowerTriang( const bool diag=false );

    //! \brief Conditionnal assignment
    //! \param test Logical test to evaluate.
    //! \param a    Value to test with matrix elements.
    //! \param tvalue If the test is true set the corresponding
    //!               matrice element to "tvalue".
    //! \param fvalue If the test is false set the corresponding
    //!               matrice element to "fvalue".
    //! \return Return the Matrix result.
    Matrix inline Compare( const char   *str_, const double a, 
		    const double t_val, const double f_val );
    //! \brief Conditionnal assignment
    //! \param test Logical test to evaluate. Only "==" available.
    //! \param B    Matrix. Evaluate the test element-by-element.
    //! \param tvalue If the test is true set the corresponding
    //!               matrice element to "tvalue".
    //! \param fvalue If the test is false set the corresponding
    //!               matrice element to "fvalue".
    //! \return Return the Matrix result.
    Matrix inline Compare( const char   *str_, const Matrix& A, 
		    const double t_val, const double f_val );
    //! \brief Remove the row i and one column j from the object matrix
    //!        if i (or j) is out of matrix bounds, nothing is done.
    //! \param i row index. 
    //! \param j row index.
    //! \return Return the Matrix result.
    Matrix inline RemoveRowCol( const size_t i, const size_t j);

    //! \brief Normalize the matrix according to the string "str_" value:
    //!        - "row" normalization : A(i,j) = A(i,j) / sum( A(i,j), j=1,..,col )
    //!        - "col" normalization : A(i,j) = A(i,j) / sum( A(i,j), i=1,..,row )
    //! \param A    Matrix to normalize.
    //! \param str_ Normalization mode ("row", "col").
    //! \return Return the Matrix result.
    Matrix inline Normalize( const Matrix& A, const char *str_ );

    //! \brief Sum by row (or column) the matrix excepted the jth colum 
    //!        (or ith row) terms:
    //!        - "row" summation : 
    //!             S(i,j) = sum( A(i,j), j=0,..,j-1,j+1,..,col-1 )
    //!        - "col" summation : 
    //!             S(i,j) = sum( A(i,j), i=0,..,i-1,i+1,..,row-1 )
    //! \param A    Initial matrix.
    //! \param str_ Summation mode ("row", "col").
    //! \return Return the vector result ( Matrix(row,1) for the "row" mode,
    //!          Matrix(1,col) for the "col" mode. 
    Matrix inline RestrictSum( const Matrix& A, const char *str_ );

    //! \brief Calculate the product :
    //!          Cij = sum( Aik*Bkj, k = 1,.., i-1,i+1,..,col)
    //!          A must be a squared matrix. 
    //! \param  A    Initial matrix.
    //! \param  B    Initial matrix.
    //! \return Return the restricted Matrix product.
    Matrix inline RestrictProd( const Matrix& A, const Matrix& B );

    //! \brief Display the matrix.
    //! \param name Matrix name.    
    void print( const char *name, size_t first_indexes=2 );

};


//   -------------------------------------
//
//   Define inlined functions and operators
//
//   -------------------------------------

//
//   Operators
//

Matrix inline Matrix::operator+( ) {

  double *tab = new double [ row_ * col_ ];
  double *mat_ = link_->matrix_;

  for (size_t i=0; i < row_*col_; i++) {
    tab[i] = mat_[i];
  }
  Matrix tmp(row_, col_, tab, false);

  return ( tmp );
}

Matrix inline Matrix::operator-( ) {

  double *tab = new double [ row_ * col_ ];
  double *mat_ = link_->matrix_;

  for (size_t i=0; i < row_*col_; i++) {
    tab[i] = - mat_[i];
  }
  Matrix tmp(row_, col_, tab, false);

  return ( tmp );
}

Matrix Matrix::operator+( const double a ) {

  double *tab = new double [ row_ * col_ ];
  double *mat_ = link_->matrix_;

  for (size_t i=0; i < row_*col_; i++) {
    tab[i] = mat_[i]+a;
  }
  Matrix tmp(row_, col_, tab, false);

  return ( tmp );
}

Matrix  Matrix::operator+(const Matrix& A) {
    
  Matrix tmp( row_,col_);

  double *tmp_tab = tmp.link_->matrix_; 
  double *_tab    = link_->matrix_;
  double *A_tab   = A.link_->matrix_;
 
  if( EQUAL_DIM( A ) ) {
    
    for (size_t i=0; i < row_*col_; i++) {
      tmp_tab[i] = _tab[i] +  A_tab[i];
    }
  }
  else {

  }
  return ( tmp );

}

Matrix Matrix::operator-( const double a ) {

  double *tab = new double [ row_ * col_ ];
  double *mat_ = link_->matrix_;

  for (size_t i=0; i < row_*col_; i++) {
    tab[i] = mat_[i] - a;
  }
  Matrix tmp(row_, col_, tab, false);

  return ( tmp );
}

Matrix inline Matrix::operator-(const Matrix& A) {
    
  Matrix tmp( row_,col_);

  double *tmp_tab = tmp.link_->matrix_; 
  double *_tab    = link_->matrix_;
  double *A_tab   = A.link_->matrix_;

  if( EQUAL_DIM( A ) ) {
    for (size_t i=0; i < row_*col_; i++) {
      tmp_tab[i] = _tab[i] -  A_tab[i];
    }
  }
  else {

  }
  return ( tmp );

}

Matrix Matrix::operator*( double a ) {

  double *tab = new double [ row_ * col_ ];
  double *mat_ = link_->matrix_;

  for (size_t i=0; i < row_*col_; i++) {
    tab[i] = mat_[i] * a;
  }
  Matrix tmp(row_, col_, tab, false );

  return ( tmp );
}

Matrix inline Matrix::operator*(const Matrix& A) {
    
  Matrix tmp(row_, A.col_) ; 
  double *tmp_tab = tmp.link_->matrix_; 
  double *_tab    = link_->matrix_;
  double *A_tab   = A.link_->matrix_;
 


  if( EQUAL_COL_ROW( A ) ) {

    size_t A_row = A.row_;

#ifndef OPT
    // Initial loops
    double sum = 0.0;
    for (size_t i=0; i < row_; i++) {
      for (size_t j=0; j < A.col_; j++) {
	sum = 0.0;
	for (size_t k=0; k < col_; k++) {
	  sum += _tab[k*row_+i] * A_tab[j*A_row + k];
	}
	tmp_tab[j*row_+i] = sum;
      }
    }
#else
    double *sum_tab = new double[row_]; 
    // Loops i and k have switched for performance reason.
    for (size_t j=0; j < A.col_; j++) {
      memset(sum_tab, 0, sizeof(double) * row_ );
      for (size_t k=0; k < col_; k++) {
	for (size_t i=0; i < row_; i++) {
	  sum_tab[i] += _tab[k*row_+i] * A_tab[j*A_row + k];
	}
      }
      for (size_t i=0; i < row_; i++) {
	tmp_tab[j*row_+i] = sum_tab[i];
      }
    }
    delete [] sum_tab;
#endif
  }
  else {

  }

  return ( tmp );
}

Matrix inline Matrix::operator/( double a ) {

  double *tab = new double [ row_ * col_ ];
  double *mat_ = link_->matrix_;

  for (size_t i=0; i < row_*col_; i++) {
    tab[i] = mat_[i]/a;
  }
  Matrix tmp(row_, col_, tab, false);

  return ( tmp );
}

Matrix inline Matrix::operator|(const Matrix& A) {
    
  Matrix tmp(row_, col_) ; 

  double *tmp_tab = tmp.link_->matrix_; 
  double *_tab    = link_->matrix_;
  double *A_tab   = A.link_->matrix_;


  if( EQUAL_DIM( A ) ) {
    for (size_t i=0; i < row_*col_; i++) {
      tmp_tab[i] = _tab[i] * A_tab[i]; 
    }
  }
  else {

  }


  return ( tmp );
}

//
//    Member functions
//    ----------------

Matrix inline Matrix::t() {
    
  size_t t_row = col_;
  size_t t_col = row_;
  Matrix tmp(t_row, t_col) ; 
  double *tmp_tab = tmp.link_->matrix_; 
  double *_tab    = link_->matrix_;



  for (size_t j=0; j < t_col; j++) {
    for (size_t i=0; i < t_row; i++) {
      tmp_tab[j*t_row+i] = _tab[i*row_+j];
    }
  }


  return ( tmp );
}

Matrix inline Matrix::getRowVector( const size_t i ) {

  const double* mat_ = getData();
  double *tmp_mat    = new double [col_];
  for (size_t j=0; j < col_; j++) {
    tmp_mat[j] = mat_[j*row_+i];
  }
  
  Matrix tmp( 1, col_, tmp_mat, false );
  return( tmp );
}

Matrix inline Matrix::getColVector( const size_t j ) {

  const double* mat_ = getData();
  double *tmp_mat    = new double [row_];

  memcpy ( tmp_mat, &mat_[j*row_], row_*sizeof(double));
  
  Matrix tmp( row_, 1, tmp_mat, false );
  return( tmp );
}

Matrix inline Matrix::getUpperTriang(const bool diag) {
  Matrix tmp = getCopy();
  double *mat_= tmp.link_->matrix_;
  size_t shift;

  if( diag ) 
    shift = 1;
  else
    shift = 0;

  for (size_t j=0; j < col_; j++) {
    for (size_t i=j+shift; i < row_; i++) {
      mat_[j*row_+i] = 0;
    }
  }
  return( tmp );
}

Matrix inline Matrix::getLowerTriang(const bool diag) {
  Matrix tmp = getCopy();
  double *mat_= tmp.link_->matrix_;
  size_t shift;
  if( diag ) 
    shift = 1;
  else
    shift = 0;
  for (size_t j=shift; j < col_; j++) {
    for (size_t i=0; i < j-1; i++) {
      mat_[j*row_+i] = 0;
    } 
  } 
  return( tmp );
}

Matrix inline Matrix::Compare( const char *str_,   const Matrix& A, 
			const double t_val, const double f_val ) {

  std::string str(str_);
  const double *mat_  = getData();
  const double *mat_A = A.getData();
  double *tmp_mat = new double[row_*col_];

  if( EQUAL_DIM( A ) ) {
    if ( str.compare("==") == 0){
      for (size_t i=0; i < row_*col_; i++) {
	if ( mat_[i] ==  mat_A[i]) 
	  tmp_mat[i] = t_val;
	else
	  tmp_mat[i] = f_val;
      }
    }
    else {

    }
  } else {

  }
  Matrix tmp( col_, row_, tmp_mat,  false);
  return ( tmp );

}

Matrix inline Matrix::Compare( const char *str_,   const double test_val, 
			const double t_val, const double f_val ) {

  std::string str(str_);
  const double *mat_  = getData();
  double *tmp_mat = new double[row_*col_];

  if ( str.compare("==") == 0){
    for (size_t i=0; i < row_*col_; i++) {
      if ( mat_[i] == test_val ) 
	tmp_mat[i] = t_val;
      else
	tmp_mat[i] = f_val;
    }
  }
  else {

  }
  Matrix tmp( col_, row_, tmp_mat, false);
  return ( tmp );

}

Matrix inline Matrix::RemoveRowCol( const size_t ii, const size_t jj) {
  size_t nrow, ncol;
  size_t ix, jx;

  if( (ii >= 0) && (ii < row_) ) {
    nrow = row_ - 1;
    ix = ii;
  } else {
    nrow = row_;
    ix = row_;
  }
  if( (jj >= 0) && (jj < col_) ) {
    ncol = col_ - 1;
    jx   = jj;
  } else {
    ncol = col_;
    jx   = col_;
  }

  const double *mat  = getData();
  double *tmp_mat = new double[nrow*ncol];


  // First block
  //
  for (size_t j=0; j < jx; j++) {
    for (size_t i=0; i < ix; i++) {
      tmp_mat[j*nrow+i] = mat[j*row_+i];
    }
  }    

  // Second block
  //
  for (size_t j=jx+1; j < col_; j++) {
    for (size_t i=0; i < ix; i++) {
      tmp_mat[(j-1)*nrow+i] = mat[j*row_+i];
    }
  }    

  // Third block
  //
  for (size_t j=0; j < jx; j++) {
    for (size_t i=ix+1; i < row_; i++) {
      tmp_mat[(j)*nrow+(i-1)] = mat[j*row_+i];
    }
  }    

  // Fourth block
  //
  for (size_t j=jx+1; j < col_; j++) {
    for (size_t i=ix+1; i < row_; i++) {
      tmp_mat[(j-1)*nrow+(i-1)] = mat[j*row_+i];
    }
  }    

  Matrix tmp(nrow, ncol, tmp_mat, false);
  return ( tmp );

}

//
//   Friend functions and operators
// 
Matrix inline operator+(const double a, const Matrix& A ) {

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getData();

  for (size_t i=0; i < row*col; i++) {
    tab[i] = a + mat_[i];
  }
  Matrix tmp(row, col, tab, false);


  return ( tmp );
}  

Matrix inline operator-(const double a, const Matrix& A ) {

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getData();

  for (size_t i=0; i < row*col; i++) {
    tab[i] = a - mat_[i];
  }
  Matrix tmp(row, col, tab, false);


  return( tmp );
}

Matrix inline operator*(const double a, const Matrix& A ) {

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getData();

  for (size_t i=0; i < row*col; i++) {
    tab[i] = a * mat_[i];
  }
  Matrix tmp(row, col, tab, false);


  return( tmp );
}

Matrix inline max( const double x, Matrix &A) {


  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getData();

  for (size_t i=0; i < row*col; i++) {
    tab[i] = MAX( mat_[i], x );
  }
  Matrix tmp(row, col, tab, false);



  return ( tmp );
}

Matrix inline max( const Matrix &A, const double x) {



  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getData();

  for (size_t i=0; i < row*col; i++) {
    tab[i] = MAX( mat_[i], x );
  }
  Matrix tmp(row, col, tab, false);



  return ( tmp );
}

Matrix inline max( const Matrix &A, Matrix &B ) {



  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *matA_ = A.getData();
  const double *matB_ = B.getData();

  if( EQUAL_DIM_MAT( A, B ) ) {
    for (size_t i=0; i < row*col; i++) {
      tab[i] = MAX( matA_[i], matB_[i]);
    }
  } else {

  }

  Matrix tmp(row, col, tab, false );

  return ( tmp );
}


Matrix inline min(  const double x, const Matrix &A) {


  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getData();

  for (size_t i=0; i < row*col; i++) {
    tab[i] = MIN( mat_[i], x );
  }
  Matrix tmp(row, col, tab, false);



  return ( tmp );
}

Matrix inline min( const Matrix &A, const double x) {



  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getData();

  for (size_t i=0; i < row*col; i++) {
    tab[i] = MIN( mat_[i], x );
  }
  Matrix tmp(row, col, tab, false);



  return ( tmp );
}

Matrix inline min( const Matrix &A, const Matrix &B ) {



  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *matA_ = A.getData();
  const double *matB_ = B.getData();

  if( EQUAL_DIM_MAT( A, B ) ) {
    for (size_t i=0; i < row*col; i++) {
      tab[i] = MIN( matA_[i], matB_[i]);
    }
  } else {

  }

  Matrix tmp(row, col, tab, false);



  return ( tmp );
}

Matrix inline sum( const Matrix &A, const char *dim  ) {

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  std::string  str(dim);

  size_t row_;
  size_t col_;
  double *tab = 0;

  const double *mat = A.getData();

  if ( (str.compare("row") != 0) && (str.compare("col") != 0) && 
       (str.compare("all") != 0)) {

  }
  if ( str.compare("col") == 0 ) {
    tab = new double [ col ];
    for (size_t j=0; j < col; j++) {
      tab[j] = 0.0;
      for (size_t i=0; i < row; i++) {
	tab[j] += mat[j*row+i];
      }
    }  
    row_ = 1;
    col_ = col;
  } else if ( str.compare("row") == 0 ) {
    tab = new double [ row ];
    for (size_t i=0; i < row; i++) {
      tab[i] = 0.0;
      for (size_t j=0; j < col; j++) {
	tab[i] += mat[j*row+i];
      }
    }  
    row_ = row;
    col_ = 1;
  } else {
    // sum all terms
    tab = new double [ 1 ];
    tab[ 0 ] = 0.0;
    for (size_t i=0; i < row*col; i++) {
      tab[0] +=  mat[i];
    }
    row_ = 1;
    col_ = 1;
  }

  Matrix tmp(row_, col_, tab, false);

  return( tmp );
}

Matrix inline log( const Matrix &A ) {

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getData();
  for (size_t i=0; i < row*col; i++) {
    tab[i] = std::log( mat_[i] );
  }
  Matrix tmp(row, col, tab, false);

  return( tmp );
}

Matrix inline exp( const Matrix &A ) {

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getData();
  for (size_t i=0; i < row*col; i++) {
    tab[i] = std::exp( mat_[i] );
  }
  Matrix tmp(row, col, tab, false);

  return( tmp );
}

Matrix inline abs( const Matrix &A ) {

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getData();
  for (size_t i=0; i < row*col; i++) {
    tab[i] = fabs( mat_[i] );
  }
  Matrix tmp(row, col, tab, false);

  return( tmp );
}

Matrix inline Normalize( const Matrix &A, const char *str_ ) {

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  const double *matA = A.getData();
  double *tab = new double [ row * col ];

  std::string str(str_);
  double sum   = 0.0;

  if ( str.compare("row") == 0) {
    for (size_t i=0; i < row; i++) {
      sum = 0.0;
      for (size_t j=0; j < col; j++) {
	sum = sum + matA[j*row+i];
      }
      if (sum < DBL_MIN ) {    

      }
      for (size_t j=0; j < col; j++) {
	tab[j*row+i] = matA[j*row+i] / sum;
      }
    }
  } else if( str.compare("col") == 0 ) {
    for (size_t j=0; j < col; j++) {
      sum = 0.0;
      for (size_t i=0; i < row; i++) {
	sum = sum + matA[j*row+i];
      }
      if (sum < DBL_MIN ) {    

      }
      for (size_t i=0; i < row; i++) {
	tab[j*row+i] = matA[j*row+i] / sum;
      }
    }
  } else {

  }

  Matrix tmp(row, col, tab, false);
  return( tmp );
}

Matrix inline RestrictSum( const Matrix& A, const char *str_ ) {

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  const double *matA = A.getData();

  std::string str(str_);
  double *tab = new double [ row * col ];   
  memset( tab, 0, sizeof(double) * row * col);


  if ( str.compare("col") == 0) {
    // i row index not considered in the sum
    for (size_t i=0; i < row; i++) {
      for (size_t j=0; j < col; j++) {
	for (size_t k=0; k < i; k++) {
	  tab[j*row+i] =  tab[j*row+i] + matA[j*row+k];
	}
	for (size_t k=i+1; k < row; k++) {
	  tab[j*row+i] =  tab[j*row+i] + matA[j*row+k];
	}
      }
    }
  } else if ( str.compare("row") == 0) {
    // j column index not considered in the sum
    for (size_t i=0; i < row; i++) {
      for (size_t j=0; j < col; j++) {
	for (size_t k=0; k < j; k++) {
	  tab[j*row+i] =  tab[j*row+i] + matA[k*row+i];
	}
	for (size_t k=j+1; k < col; k++) {
	  tab[j*row+i] =  tab[j*row+i] + matA[k*row+i];
	}
      }
    }
  } else {

  }

  Matrix tmp(row, col, tab, false);
 
  return( tmp );
}

Matrix inline RestrictProd(const Matrix& B, const Matrix& A) {
    
  // Computation of B*A
  size_t row_ = B.getRow();
  size_t col_ = B.getCol();
  const double *_tab    =  B.getData();
  const double *A_tab   =  A.getData();
  double *tmp_tab =  new double [ row_ * A.getCol() ];
 


  if( col_ == A.getRow() )  {

    // ???
    if( EQUAL_DIM_MAT( A, A ) ) {

      size_t A_row = A.getRow();
      size_t A_col = A.getCol();

      // Initial loops
      double sum = 0.0;
      for (size_t i=0; i < row_; i++) {
	for (size_t j=0; j < A_col; j++) {
	  sum = 0.0;
	  for (size_t k=0; k < i; k++) {
	    sum += _tab[k*row_+i] * A_tab[j*A_row + k];
	  }
	  for (size_t k=i+1; k < col_; k++) {
	    sum += _tab[k*row_+i] * A_tab[j*A_row + k];
	  }
	  tmp_tab[j*row_+i] = sum;
	}
      }
    } else {

    }
  } else {

  }

  Matrix tmp(row_,  A.getCol(), tmp_tab, false);
  return ( tmp );
}
  
#endif
