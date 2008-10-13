/* BandMemMatrix.h
 *
 * Copyright (C) 2006 Laboratoire Statistique & Génome
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*!
 *   \file
 *     \author Vincent Miele
 *       \brief Band Memory Matrix class
 *
 *          13/03/2006
 *          
*/


#ifndef ERMG_BANDMEMMATRIX_H
#define ERMG_BANDMEMMATRIX_H
/*!
  \class BandMemMatrix libermg/BandMemMatrix.h
  \brief BandMem Matrix class
*/


#include <vector>
#include <algorithm>

namespace ermg {

  template <typename Content> 
    class BandMemMatrix
    { 
    protected:
      typedef std::vector<Content> Band;

      //! matrix as a vector of Content
      std::vector<Content> _matrix;
      //! the 1st dimension
      int _isize;
      //! the 2st dimension
      long _jsize;
      //! begin of _matrix
      typename Band::iterator _itbeg;
      //! end of _matrix
      typename Band::iterator _itend;

    public:
      typedef typename Band::iterator BandCursor;

      //! default
      BandMemMatrix(){
	_isize=0; _jsize=0;
      }

      //! constructor for isize lines, jsize columswith the default value def_val
      BandMemMatrix(int isize, int jsize, double def_val=Content())
	{
	  _isize = isize;      
	  _jsize = jsize;     
	  _matrix.assign(isize*jsize, def_val);
	  _itbeg = _matrix.begin();
	  _itend = _matrix.end();
	}

      //! copy constructor
      /*
	/param bmm BandMemMatrix to copy
      */
      BandMemMatrix(const BandMemMatrix& bmm)
	{   
	  _isize = bmm._isize;      
	  _jsize = bmm._jsize;     
	  _matrix = bmm._matrix;
	  _itbeg = _matrix.begin();
	  _itend = _matrix.end();
	}
  
      //! copy operator
      /*
	/param bmm BandMemMatrix to copy
      */
      BandMemMatrix& operator=(const BandMemMatrix& bmm)
	{    
	  _isize = bmm._isize;      
	  _jsize = bmm._jsize;     
	  _matrix = bmm._matrix;
	  _itbeg = _matrix.begin();
	  _itend = _matrix.end();
	  return *this;    
	}
  
      //! destructor
      ~BandMemMatrix(){}


      //! access to the iline-th line, jcol-th column
      /*
	/param iline concerned line number
	/param jcol concerned column number
      */
      Content& value(int iline, int jcol){
	return _matrix[iline*_jsize+jcol];
      }

      //! access to iline-th line, first column
      /*
	/param iline concerned line number
      */  
      BandCursor access(int iline){
	return _itbeg+(iline*_jsize);
      }
  
      //! access to iline-th line, jcol-th column
      /*
	/param iline concerned line number
      */
      BandCursor access(int iline, int jcol){
	return _itbeg+(iline*_jsize)+jcol;
      }

      //! assign val to all elements
      void assign(Content val){
	_matrix.assign(_isize*_jsize, val);
	_itbeg = _matrix.begin();
	_itend = _matrix.end();
      }   

      //! assign val to all elements of matrix isizeXjsize 
      void assign(int isize, int jsize, Content val){
	_isize = isize; _jsize = jsize;
	_matrix.assign(isize*jsize, val);
	_itbeg = _matrix.begin();
	_itend = _matrix.end();
      }
  
      //! fill with val 
      void fill(Content val){
	std::fill(_matrix.begin(), _matrix.end(), val); 
      }

      //! Band begin
      BandCursor begin(){
	return _itbeg;
      }

      //! Band end
      BandCursor end(){
	return _itend;
      }


      //! return the 1st dimension
      int nLines() const{
	return _isize;
      }
      //! return the 2st dimension
      int nCols() const{
	return _jsize;
      }  
    };

}
#endif


