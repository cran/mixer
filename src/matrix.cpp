//  matrix.cpp
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

//! \file    matrix.cpp
//! \author  Gilles Grasseau
//! \brief   Basic tools for matrix
//! \version 0.1
//! \date    June 2008

//! \cond
#include "matrix.h"

//
//    Member functions
//    ----------------

Matrix Matrix::getCopy() {
    
  double *tab = new double[row_*col_];
  memcpy ( tab, link_->matrix_, row_*col_*sizeof(double));

  Matrix tmp(row_, col_, tab, false) ; 

  return ( tmp );
}

void Matrix::print( const char *name, size_t first_indexes ) { 
 
  double *mat_ = link_->matrix_;
  
  
  size_t max_row = MIN( row_, first_indexes);
  size_t max_col = MIN( col_, first_indexes);
  
  for (size_t i=0; i < max_row; i++) {

    
    for (size_t j=0; j < max_col; j++) {

    }
    if( col_ > max_col ) {

    }

  }
  if( row_ > max_row ) {

    
    for (size_t j=0; j < max_col; j++) {

    }
    if( col_ > max_col ) {

    }

  }
}
//! \endcond


