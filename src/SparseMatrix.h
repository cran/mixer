/* SparseMatrix.h
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
 *       \brief Sparse Matrix class
 *
 *          13/03/2006
 *          
 */

#ifndef ERMG_SPARSEMATRIX_H
#define ERMG_SPARSEMATRIX_H
/*!
  \class SparseMatrix libermg/Sparsematrix.h
  \brief Sparse Matrix
*/

#include <cstdlib>

#include <string>
#include <vector>
#include <deque>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sstream>

extern "C" { 
  void dnaupd_( int& ido, char* bmat, int& n, char* which, int& nev, double& tol, 
		double* resid, int& ncv, double* v, int& ldv, int* iparam,
		int* ipntr, double* workd, double* workl, int& lworkl, int& info );
  void dneupd_( int& rvec, char* howmany, int* select, double* dr, double* di, 
		double*z, int& ldz, double&sigmar, double& sigmai, double* workev, char*bmat, 
		int& n, char* which, int& nev, double& tol, double* resid, 
		int& ncv, double* v, int& ldv, int* iparam, int* ipntr, double* workd, double* workl, 
		int& lworkl, int& info );
}



namespace ermg {
  
  template < typename Content, typename Container = std::vector<Content> >
  class SparseMatrix
  {
  protected:
    //! matrix
    std::vector< Container > _matrix;
    //! position in each line of the first non null element after the diagonal
    std::vector<int> _first_position_after_diagonal;
    //! checks ifthe diagonal is null or not
    std::deque<bool> _full_diagonal;
    //! 2nd eigen vector
    std::vector<double> _eigv2;

    void laplacianMatrixProduct(double* vecin, double* vecout)
    {  
      //std::cerr<<"[";	
      for (int i=0; i<int(_matrix.size()); i++){
	double sum = 0.;
	int s = _matrix[i].size();	  
	typename Container::iterator it_l = _matrix[i].begin();
	int pos = _first_position_after_diagonal[i];
	// avoiding matrix[i][i]
	if (fullDiagonal(i)){
	  pos--;   
	} 

	int ind1=-1, ind2;
	for (int nbe=0; nbe<pos; nbe++){
	  sum -= vecin[*it_l]; // -1 if connected
	  ind2 = *it_l;
	  //for (int j=ind1+1; j<ind2; j++){
	  //std::cerr<<"0 ";
	  //}
	  //std::cerr<<"-1 ";
	  ind1 = ind2;
	  it_l++;
	}
	//for (int j=ind1+1; j<i; j++){
	//std::cerr<<"0 ";
	//}
	if (fullDiagonal(i)){
	  sum += sqrt( static_cast<double> (s-1) )*vecin[i]; // degree
	  //sum += (s-1)*vecin[i]; // degree
	  //std::cerr<<s-1<<" ";
	  pos++;
	  it_l++;
	}
	else{
	  // GG Nov-2009
	  // sum += sqrt( s )*vecin[i]; // degree
	  sum += sqrt( static_cast<double> (s) )*vecin[i]; // degree
	  //sum += s*vecin[i]; // degree
	}
	    
	ind1=i;
	for (int nbe=pos; nbe<s; nbe++){
	  sum -= vecin[*it_l]; // -1 if connected
	  ind2 = *it_l;
	  //for (int j=ind1+1; j<ind2; j++){
	  //std::cerr<<"0 ";
	  //}
	  //std::cerr<<"-1 ";
	  ind1 = ind2;
	  it_l++;
	}
	//for (int j=ind1+1; j<_matrix.size(); j++){
	//std::cerr<<"0 ";
	//}
	vecout[i] = sum;
	//if (i<_matrix.size()-1)
	//std::cerr<<",";
	//else
	//std::cerr<<"]";
	//std::cerr<<std::endl;
      }
    }
 
  public:
    //! default constructor
    SparseMatrix(){}

    //! constructor with nline empty lines
    SparseMatrix(int nline)
      : _matrix(nline, Container())
    {}
 
    //! constructor from a vector
    /*
      /param spmvec sparse matrix vector containing the coordinates of the non null elements [numline,numcol]
    */
    SparseMatrix(const std::vector< Container >& spmmat )
    {          
      for (int numline=0; numline<int(spmmat.size()); numline++){
	int numcol;
	Container tmp;
	bool seeking_first_position_after_diagonal = true;
	bool nonulldiag = false;
	if (spmmat[numline].size()!=0){
	  for (unsigned int c=0; c<spmmat[numline].size(); c++){
	    numcol=spmmat[numline][c];
	    tmp.push_back(numcol);
	    if (seeking_first_position_after_diagonal){
	      if (numcol>numline){	      
		_first_position_after_diagonal.push_back( tmp.size()-1 );
		seeking_first_position_after_diagonal = false;
	      }
	      if (numcol==numline)
		nonulldiag = true;
	    }
	  }
	}
	_full_diagonal.push_back(nonulldiag);
	_matrix.push_back(tmp);
	if (seeking_first_position_after_diagonal)
	  _first_position_after_diagonal.push_back( tmp.size() );
      }
    }

 
    //! copy constructor
    /*
      /param spm SparseMatrix to copy
    */
    SparseMatrix(const SparseMatrix& spm)
    {
      _matrix = spm._matrix;
      _first_position_after_diagonal = spm._first_position_after_diagonal;
    }
      
 
    //! copy operator
    /*
      /param spm SparseMatrix to copy
    */
    SparseMatrix& operator=(const SparseMatrix& spm)
    {
      _matrix = spm._matrix;
      _first_position_after_diagonal = spm._first_position_after_diagonal;
      return *this;
    }

    //! submatrix copy constructor
    /*
      /param spm SparseMatrix copy into a submatrix
    */
    SparseMatrix(const SparseMatrix& spm, const std::vector<int> concerned_vertices)
    {
      _matrix.assign(concerned_vertices.size(), Container());
	
      for (int i_cv=0; i_cv<int(concerned_vertices.size()); i_cv++){	
	int numline=concerned_vertices[i_cv];
	int numcol;	  
	
	Container tmp;
	Container& line = const_cast<SparseMatrix&>(spm)[numline];	  
	int j_cv = 0;
	std::vector<int>::const_iterator it_cv = concerned_vertices.begin();
	typename Container::iterator it_l = line.begin();
	
	bool seeking_first_position_after_diagonal = true;
	bool nonulldiag = false;
	while (it_l != line.end()){
	  bool stop = false;
	  while (!stop){
	    if (it_cv!=concerned_vertices.end()){
	      if (*it_cv<*it_l){
		it_cv++;
		j_cv++;
	      }
	      else
		stop = true;
	    }
	    else 
	      stop = true;
	  }
	  if (it_cv==concerned_vertices.end())
	    it_l = line.end();
	  else{
	    if (*it_cv == *it_l){
	      numcol = *it_cv; 
	      tmp.push_back(j_cv);// vertices number in the submatrix
	      if (seeking_first_position_after_diagonal){
		if (numcol>numline){
		  _first_position_after_diagonal.push_back(tmp.size()-1);
		  seeking_first_position_after_diagonal = false;
		}
		if (numcol==numline)
		  nonulldiag = true;	      }
	    }
	    it_l++;
	  }
	}	
	_full_diagonal.push_back(nonulldiag);
	//_matrix.push_back(tmp);
	_matrix[i_cv] = tmp;
	if (seeking_first_position_after_diagonal)
	  _first_position_after_diagonal.push_back(tmp.size());
      }
    }
      

  //! destructor
    ~SparseMatrix(){}


  //! access to the iline-th line
  /*
    /param iline concerned line number
  */
  Container& operator[] (int iline){
    return _matrix[iline];
  }

  //! returns the column of the first non null element after the diagonal
  int firstPositionAfterDiagonal(int iline){
    return _first_position_after_diagonal[iline];
  }
 
  //! returns false if the diagonal is null
  int fullDiagonal(int iline){
    return _full_diagonal[iline];
  }

  //! return the number of lines
  int nLines() const{
    return _matrix.size();
  }

  //! print
  void stderrPrint() const{
    std::cerr<<"SparseMatrix:"<<std::endl;
    for (int i=0; i<int(_matrix.size()); i++){
      //std::cerr<<_matrix[i].size()<<" => ";
      for (int j=0; j<int(_matrix[i].size()); j++){
	std::cerr<<_matrix[i][j]<<" ";
      }
      std::cerr<<std::endl;
    }
      exit(1);
  }


  //! computes and return second eigen-vector
  std::vector<double>& secondEigenVector()
  {
    if (_eigv2.empty()){	  
      int n = _matrix.size();
      //std::cerr<<"size: "<<n<<std::endl;
      int ido = 0; 
      char* bmat = new char[1];
      bmat[0] = 'I';
      char* which = new char[2] ;
      which[0]='L'; which[1]='M'; 
      int nev = 2; // <n-1
      double tol = 0; 
      double* resid = new double[n];
      for (int i=0; i<n; i++)
	resid[i]=0;
      int ncv = 2*nev+2; // HARD TO CHOOSE
      double** v = new double*[ncv]; // 
      // array of size n*ncv stored in column
      v[0] = new double[ncv*n];
      //std::cerr<<"Allocation v of size :"<<ncv*n*8<<std::endl;
      for(int i=1;i<ncv;i++)
	v[i]=v[i-1]+n;    
      int ldv = n; 
      int* iparam = new int[11];
      int maxitr = 300;
      iparam[0] = 1;
      iparam[2] = maxitr;
      iparam[6] = 1;
      int* ipntr = new int[14]; 
      double* workd = new double[3*n];
      double* workl = new double[3*ncv*ncv+6*ncv];
      int lworkl = 3*int(pow(double(ncv), 2))+6*ncv; 
      int info = 0;
	  
      bool stop = false;
      while (!stop){
	dnaupd_( ido, bmat, n, which, nev, tol, 
		 //resid, ncv,  &v[0][0],  ldv, iparam, 
		 resid, ncv,  v[0],  ldv, iparam,
		 ipntr, workd,  workl, lworkl, info);
	if ((ido==-1)||(ido==1)){
	  // WARNING: -1 because c++
	  laplacianMatrixProduct(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
	} 
	else
	  stop = true;
      }
	  
      if (info<0){
	//std::cerr<<"\ndnaupd_: problem "<<info<<std::endl;
	_eigv2 = std::vector<double>();
      }
      else{
	int rvec = 1;
	char* howmany = new char[1]; 
	howmany[0] = 'A';  
	int* select = new int[ncv];  
	double* workev = new double[3*ncv];
	double* dr = new double[ncv];
	double* di = new double[ncv];
	double sigmar, sigmai;
	int ierr;
	    
	//dneupd_( rvec, howmany, select, dr, di, &v[0][0], ldv, 
	dneupd_( rvec, howmany, select, dr, di, v[0], ldv, 
		 sigmar, sigmai, workev, bmat, n, which, nev, tol, 
		 //resid, ncv, &v[0][0], ldv, iparam, ipntr, workd, workl, 
		 resid, ncv, v[0], ldv, iparam, ipntr, workd, workl,
		 lworkl, ierr ); 
#ifdef VERBOSE
	std::cerr<<"Ierr:"<<ierr<<std::endl;
	std::cerr<<"Eigen values:"<<std::endl;
	for (int k=0; k<nev; k++)
	  std::cerr<<dr[k]<<" ";
	std::cerr<<std::endl;
	    
	// fortran columns of v, i.e. c++ lines of v contain eigenvectors
	for (int k=0; k<nev; k++){
	  for (int j=0; j<n; j++)
	    std::cerr<<v[k][j]<<" ";
	  std::cerr<<std::endl;
	}
#endif
    
	_eigv2.reserve(n);
	for (int j=0; j<n; j++)
	  _eigv2.push_back(v[1][j]);
	    
	delete[] howmany;
	delete[] select;
	delete[] workev;
	delete[] dr;
	delete[] di;     
      }
	  
      delete[] bmat;
      delete[] which;
      delete[] resid;
      delete[] v[0];
      delete[] v;
      delete[] iparam;
      delete[] ipntr;
      delete[] workd;
      delete[] workl; 
    }	
    return _eigv2;
  }
};
  
}
#endif
