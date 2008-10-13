/* EmCore.h
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
 *       \brief EM algorith
 *
 *          10/01/2007
 *          
*/

#ifndef ERMG_EMCORE_H
#define ERMG_EMCORE_H
/*!
  \class Ermg libermg/EmCore.h
  \brief EM algorithm core
*/


#include <Ermg.h>

const int MININITNBV=200;

namespace ermg {

  class EmCore
  {
  protected:    
    //! model to estimate
    Ermg& _ermg;
    
    //! vertices number
    int _n;
    //! vertices number considered to intialize
    int _n_init;
    //! classes number
    int _q;

    //! epsilon for the convergency of the em
    double _em_eps;
    //! maximum number of iteration of the em
    int _em_nitermax;
    //! copy of the maximum number of iteration of the em
    int _back_em_nitermax;

    //! iterator points to Xadj[i]
    std::vector<int>::iterator _it_Xadj_i;
    //! pointer to transposed Tau at iteration k
    BandMemMatrix<double>* _t_Tau1;
    //! iterator points to _t_Tau1[l,j]
    BandMemMatrix<double>::BandCursor _it_t_Tau1_lj;
    //! auxiliary variable
    std::vector<double> _sumsTau_qi;


    //! alpha at previous em iteration
    std::vector<double> _Alpha_back;
    //! Pi at previous em iteration
    BandMemMatrix<double> _Pi_back;
    //! iterator points to _Pi[q,l]
    BandMemMatrix<double>::BandCursor _it_Pi_ql;
    //! iterator points to _Pi[q,l]
    BandMemMatrix<double>::BandCursor _it_Pi_back_ql;

    //! denominator for Pi calculus
    BandMemMatrix<double> _denom_Pi;
    //! iterator points to _denom_Pi[q,l]
    BandMemMatrix<double>::BandCursor _it_denom_Pi_ql;
  
    //! possible values of Beta_ijql when Xadj_ij=1
    BandMemMatrix<double> _val1;
    //! iterator points to _val1[q,l]
    BandMemMatrix<double>::BandCursor _it_val1_ql;
    //! possible values of Beta_ijql when Xadj_ij=0
    BandMemMatrix<double> _val2;
    //! iterator points to _val2[q,l]  
    BandMemMatrix<double>::BandCursor _it_val2_ql;

    //! init the beta_ijql
    void initBetaijql();
        //! Pi and denomPi to null
    void initPiDenomPi();

    //! initialization
    virtual bool init( );  
    //! e step
    virtual void eStep() = 0;
    //! first m step
    virtual void mFirstStep() = 0;
    //! m step
    virtual void mStep() = 0;
    //! m step and likelihood step
    virtual void mStepMaximization() = 0;
    //! m step maximization for Pi
    virtual void mStepMaximizationPi() = 0;
    //! Likelhood computation
    virtual void mLikelihood() = 0; 


    //! check the convergency on the alpha and Pi
    bool emConvergency();
    

  public:
    EmCore( Ermg& ermg,
	    int em_nitermax, double em_eps=PRECISION )
      : _ermg(ermg) 
      {
	_em_nitermax = em_nitermax;
	_em_eps = em_eps;
	_t_Tau1 = NULL;
	_n = _ermg._n;
	_q = -1;
	_n_init=_n;
      }
      
      virtual ~EmCore(){};
      
      //! performs the em step
      virtual bool run(int nbclass) = 0;

      // returns the number of vertices to be initialized
      virtual int requiredNumberOfInitializedVertices() const{
	return _n_init;
      }

      // sets the number of vertices to be initialized
      void setRequiredNumberOfInitializedVertices(int n_init){
	_n_init = n_init;
      }

      //! set minimal nb iteration
      virtual void setMinimalNiter(){
	_back_em_nitermax = _em_nitermax;
	_em_nitermax = 1;
      }
      //! unset minimal nb iteration
      virtual void unsetMinimalNiter(){      
	_em_nitermax = _back_em_nitermax;
      }

      //! returns the number of em iterations
      int nitermax() const{
	return _em_nitermax;
      }
      //! sets the number of em iterations
      void setNitermax(int em_nitermax){
      _em_nitermax = em_nitermax;
      }
  };
  
}
#endif
