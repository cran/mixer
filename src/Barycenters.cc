/* Barycenters.cc
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

#include <Barycenters.h>
using namespace std;


namespace ermg {

  void Barycenters::init(int nbclass, int n)
  {
    _nbclass = nbclass;
    _n = n;
    if (_G1)
      delete _G1;
    if (_G2)
      delete _G2;
    _G1 = new BandMemMatrix<double>(_nbclass, _n); 
    _G2 = new BandMemMatrix<double>(_nbclass, _n); 
  }

  void  Barycenters::initFromSelected( const vector<int>& selected )
  {
    int s, ind1, ind2;
    _it_G1_cj = (*_G1).access(0);
    for (int c=0; c<_nbclass; c++){
      s = _Xadj[ selected[c] ].size();
      if (s>0){  
	_it_Xadj_i = _Xadj[ selected[c] ].begin();
	ind1=-1;
	for (int nbe=0; nbe<s; nbe++){
	  ind2 = *_it_Xadj_i;
	  // all the null elements between
	  // test if barycenters for partial _Xadj
	  if (ind2>_n){
	    ind2=_n;
	  }
	  for (int j=ind1+1; j<ind2; j++){
	    //*_it_G1_cj = 0;
	    (*_G1).value(c, j) = 0;
	    _it_G1_cj++;
	  }
	  
	  // test if barycenters for partial _Xadj
	  if (ind2<_n){
	    // non null element j=ind2
	    *_it_G1_cj = 1;
	    _it_G1_cj++;
	    //(*_G1).value(c, ind2) = 1;
	  }
	  ind1 = ind2;
	  _it_Xadj_i++;
	} 
	// last null elements
	for (int j=ind1+1; j<_n; j++){
	  *_it_G1_cj = 0;
	  _it_G1_cj++;
	  //(*_G1).value(c, j) = 0;
	} 
      }
      else{
	for (int j=0; j<_n; j++){
	  *_it_G1_cj = 0;
          _it_G1_cj++;
          //(*_G1).value(c, j) = 0;
	}
      }
    }
    _G2->fill( 0. );
  }

  void Barycenters::initFromClass( const vector<int>& cl,
                                   const vector<int>& cardinal_class )
  {
    int s, ind1, ind2;

    _G2->fill( 0. );
    for (int i=0; i<_n; i++){
      s = _Xadj[i].size();

      if (s>0){
        _it_Xadj_i = _Xadj[i].begin();
        ind1=-1;
        for (int nbe=0; nbe<s; nbe++){
          ind2 = *_it_Xadj_i;

          // all the null elements between
          //(*_G2).value(_class[i], j) += 0;

          // non null element j=ind2
          (*_G2).value(cl[i], ind2) += 1;
          ind1 = ind2;
          _it_Xadj_i++;

          // last null elements
          //(*_G2).value(_class[i], j) += 0;
        }
      }
    }
    this->normalizeAndConvergency(cardinal_class, false);
    this->switchG2toG1();
  }

  void Barycenters::kmeansUpdate(int i, int min_c)
  {
    int s = _Xadj[i].size(), ind1, ind2;

    if (s>0){
      _it_Xadj_i = _Xadj[i].begin();
      ind1=0; 
    
      for (int nbe=0; nbe<s; nbe++){
	ind2 = *_it_Xadj_i;
	// non null element j=ind2
	(*_G2).value(min_c,ind2) += 1.;
	_it_Xadj_i++;	  
      }
    } 
  }

  double Barycenters::kmeansDistance(int i, int c)
  {
    int s = _Xadj[i].size(), ind1, ind2; 
    double score_ic = 0;

    if (s>0){
      _it_Xadj_i = _Xadj[i].begin();
      _it_G1_cj = (*_G1).access( c );    
      ind1=-1; 
    
      for (int nbe=0; nbe<s; nbe++){
	ind2 = *_it_Xadj_i;
	// all the null elements between
	for (int j=ind1+1; j<ind2; j++){
	  double tmp = -*_it_G1_cj;
	  _it_G1_cj++;
	  //double tmp = -(*_G1).value(c,j);
	  score_ic += tmp*tmp;
	}
      
	// non null element j=ind2
	double tmp = (1-*_it_G1_cj);
	_it_G1_cj++;
	//double tmp = (1-(*_G1).value(c,ind2));
	score_ic += tmp*tmp;
	ind1 = ind2;
	_it_Xadj_i++;	  
      }
      // last null elements
      for (int j=ind1+1; j<_n; j++){
	double tmp = *_it_G1_cj;
	_it_G1_cj++;
	//double tmp = -(*_G1).value(c,j);
	score_ic += tmp*tmp;
      }
    }
  
    return score_ic;
  }

  bool Barycenters::normalizeAndConvergency( const vector<int>& cardinal_class, 
					     bool check_convergency )
  {
    _it_G1_cj = (*_G1).access(0);
    _it_G2_cj = (*_G2).access(0);
    bool convergency = true;
    for (int c=0; c<_nbclass; c++){
      if (cardinal_class[c]==0){
	for (int j=0; j<_n; j++){
	  *_it_G2_cj = 0;
	  if (check_convergency)
	    if ( fabs(*_it_G1_cj - *_it_G2_cj)> PRECISION )
	      convergency = false;
	  _it_G1_cj++;
	  _it_G2_cj++;
	}
      }
      else{
	for (int j=0; j<_n; j++){
	  //(*_G2).value(c,j) /= _cardinal_class[c];
	  //if (check_convergency)
	  //if (((*_G1).value(c,j) != (*_G2).value(c,j))) convergency = false;
	  *_it_G2_cj /= cardinal_class[c];
	  if (check_convergency)
	    if ( fabs(*_it_G1_cj - *_it_G2_cj)> PRECISION )
	      convergency = false;
	  _it_G1_cj++;
	  _it_G2_cj++;
	}
      }
    }
    return convergency;
  }

  void Barycenters::switchG2toG1()
  {
    BandMemMatrix<double>* tmp = _G1;
    _G1 = _G2;
    _G2 = tmp;
    _G2->fill( 0. );
  }


  void Barycenters::cahUpdate( int min_c1, int cardinal_min_c1, 
			       int min_c2, int cardinal_min_c2 )
  {
    int card = cardinal_min_c1+cardinal_min_c2;
    _it_G1_c1j = (*_G1).access(min_c1); 
    _it_G1_c2j = (*_G1).access(min_c2);
    for (int j=0; j<_n; j++){
      *_it_G1_c1j = (*_it_G1_c1j*cardinal_min_c1 + *_it_G1_c2j*cardinal_min_c2 )/(card); 
      _it_G1_c1j++; _it_G1_c2j++;
      //(*_G1).value(min_c1,j) = ( (*_G1).value(min_c1,j)*cardinal_min_c1 + (*_G1).value(min_c2,j)*cardinal_min_c2 )/(card);
    }
  }

  double Barycenters::cahDistance( int c1, int c2 )
  {
    double score_c1c2 = 0;
    _it_G1_c1j = (*_G1).access(c1);
    _it_G1_c2j = (*_G1).access(c2); 
    for (int j=0; j<_n; j++){
      double tmp =  (*_it_G1_c1j-*_it_G1_c2j);
      _it_G1_c1j++; _it_G1_c2j++;
      //double tmp = (*_G1).value(c1,j) - (*_G1).value(c2,j);
      score_c1c2 += tmp*tmp;
    }
    return score_c1c2;
  }

}
