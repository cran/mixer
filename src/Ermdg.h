/* Ermdg.h
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
 *       \brief Erdös Reni Mixture for Directed Graphs
 *
 *          13/10/2006
 *          
*/

#ifndef ERMG_ERMDG_H
#define ERMG_ERMDG_H
/*!
  \class Ermdg src/Ermdg.h
  \brief Erdös Reni Mixture for Directed Graphs class
*/

#include <Ermg.h>

namespace ermg {

  class Ermdg : public Ermg
  {   
  protected:
    //! graph as a sparse matrix = _Xadj transposed
    SparseMatrix< int, std::vector<int> > _tXadj;

    //! Barycenter for _Xadj 
    Barycenters& _Gplus;
    //! Barycenter for _tXadj 
    Barycenters _Gminus;

    //! forbidden default constructor
    Ermdg();
    
    //! initialize the Barycenters object
    virtual void initBarycenters(int nbclass, int n){
      _Gplus.init(nbclass, n);
      _Gminus.init(nbclass, n);
    }
    //! computes the barycenters from class
    virtual void initBarycentersFromClass(){
      _Gplus.initFromClass(_class, _cardinal_class);
      _Gminus.initFromClass(_class, _cardinal_class);    
    }
    //! computes the kmeans distance vertex i and barycenter c
    virtual double kmeansDistance(int i, int c){ 
      return (_Gplus.kmeansDistance(i, c) + _Gminus.kmeansDistance(i, c));
    }
    //! update the barycenter min_c
    virtual void kmeansUpdateBarycenters(int i, int min_c){
      _Gplus.kmeansUpdate(i, min_c);
      _Gminus.kmeansUpdate(i, min_c);
    }
    //! normalizes and check convergency
    virtual bool kmeansNormalizeBarycentersAndConvergency( bool check_convergency ){  
      return ( (_Gplus.normalizeAndConvergency(_cardinal_class, check_convergency))
	       || (_Gminus.normalizeAndConvergency(_cardinal_class, check_convergency)) );
    }
    //! switch last barycenters to previous
    virtual void kmeansSwitchBarycenters(){
      _Gplus.switchG2toG1();
      _Gminus.switchG2toG1();
    }
    //! initiallize barycenters from selected vertices
    virtual void initBarycentersFromSelected(const std::vector<int>& selected){
      _Gplus.initFromSelected(selected);
      _Gminus.initFromSelected(selected);
    }

    //! cah distance between classes c1 and c2
    virtual double cahDistance(int c1, int c2){
      return (_Gplus.cahDistance(c1,  c2) + _Gminus.cahDistance(c1,  c2));
    }
    //! update the barycenter min_c1 from itself and min_c2
    virtual void cahUpdateBarycenters(int min_c1, int cardinal_min_c1, int min_c2, int cardinal_min_c2){
      _Gplus.cahUpdate(min_c1, cardinal_min_c1, min_c2, cardinal_min_c2);
      _Gminus.cahUpdate(min_c1, cardinal_min_c1, min_c2, cardinal_min_c2);
    }
  

  public:
    //! constructor from spm and its transposed matrix tspm
    Ermdg( const std::vector< std::vector<int> >& spm, const std::vector< std::vector<int> >& tspm, 
	   bool enable_loop,
	   int kmeans_nitermax,
	   int kmeans_nbclass, 
	   int cah_minnbclass, int cah_maxnbclass )
      : Ermg(spm,enable_loop,kmeans_nitermax,kmeans_nbclass,cah_minnbclass,cah_maxnbclass), 
      _tXadj( tspm ), _Gplus(_G), _Gminus(_tXadj)
      {
#ifdef VERBOSE
	_tXadj.stderrPrint();
#endif
      }
      
      //! destructor
      virtual ~Ermdg(){};

      //! returns the ICL
      virtual double ICL(){ 
	return this->completeLikelihood() - (_q-1)/2.0*log(double(_n)) - _q*_q/2.0*log(double(_n*_n));
      }
      //! returns the BIC
      virtual double BIC(){ 
	return this->incompleteLikelihoodApproximation() - (_q-1)/2.0*log(double(_n)) - _q*_q/2.0*log(double(_n*_n));
      }  

      friend class Emd;    
  };

}
#endif
