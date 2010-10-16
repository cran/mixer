/* Barycenters.h
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
 *       \brief Erdös Reni Mixture for Graphs
 *
 *          13/03/2006
 *          
*/

#ifndef ERMG_BARYCENTERS_H
#define ERMG_BARYCENTERS_H
/*!
  \class Barycenters libermg/Barycenters.h
  \brief Barycenters
*/

#include <SparseMatrix.h>
#include <BandMemMatrix.h>
#include <cmath>
#include <cstdlib>

namespace ermg {

const double PRECISION=1e-10;

class Barycenters
{
 protected:
  //! number of classes
  int _nbclass;
  //! number of vertices
  int _n;

  //! graph as a sparse matrix
  SparseMatrix< int, std::vector<int> > & _Xadj;
  //! iterator points to Xadj[i]
  std::vector<int>::iterator _it_Xadj_i;

  //! barycenters at iteration k+1
  BandMemMatrix<double> *_G1;
  //! iterator points to _G1
  BandMemMatrix<double>::BandCursor _it_G1_cj;
  //! auxiliary iterator points to _G1
  BandMemMatrix<double>::BandCursor _it_G1_c1j;
  //! auxiliary iterator points to _G1
  BandMemMatrix<double>::BandCursor _it_G1_c2j;

  //! barycenters at iteration k+1h
  BandMemMatrix<double> *_G2;
  //! iterator points to _G2
  BandMemMatrix<double>::BandCursor _it_G2_cj;
  

 public:
  //! constructor
  Barycenters( SparseMatrix< int, std::vector<int> > & Xadj ) 
    : _Xadj(Xadj) 
    {_G1=NULL; _G2=NULL;}
  
  //! destructor
  ~Barycenters(){
    if (_G1)
      delete _G1;
    if (_G1)
      delete _G2;
    _G1=NULL; _G2=NULL; 
  };

  //! switch G2 to G1
  void switchG2toG1();
  //! memory init
  void init(int nbclass, int n);
  //! initiallize from selected vertices
  void initFromSelected( const std::vector<int>& selected );
  //! computes the barycenters from class
  void initFromClass(const std::vector<int>& cl, const std::vector<int>& cardinal_class);
  //! normalizes and check convergency
  bool normalizeAndConvergency( const std::vector<int>& cardinal_class,
				bool check_convergency );
  //! computes the kmeans distance vertex i and barycenter c
  double kmeansDistance(int i, int c);
  //! update the barycenter min_c
  void kmeansUpdate(int i, int min_c);
  //! cah distance between classes c1 and c2
  double cahDistance(int c1, int c2);
  //! update the barycenter min_c1 from itself and min_c2
  void cahUpdate(int min_c1, int cardinal_min_c1, int min_c2, int cardinal_min_c2);
};

}
#endif
