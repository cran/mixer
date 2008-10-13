/* OCEm.h
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
 *       \brief Online Classification EM algorithm
 *
 *          10/01/2007
 *          
*/

#ifndef ERMG_OCEM_H
#define ERMG_OCEM_H
/*!
  \class OCEm libermg/OCEm.h
  \brief Online Classification EM algorithm
*/

#include <limits>
#include <EmCore.h>


namespace ermg {

  class OCEmT : public EmCore
  {
  protected:
    //! current index
    int _curr_i;
    //! current class[i] before estep change class
    int _curr_old_classi;

    //! iterator points to _class
    std::vector<int>::iterator _it_class;

    virtual void sufficientStatInit()=0;
    virtual void nqIncremente(int c)=0;
    virtual void nqDecremente(int c)=0;
    virtual void nqlIncremente(int c1, int c2)=0;
    virtual void nqlDecremente(int c1, int c2)=0;


    virtual bool init();

    //! e step
    virtual void eStep();
    //! E step iteration for vertex i
    virtual void eStepVertex(int i, bool stop_at_diagonal=false)=0;
    //! core calculus of E step iteration
    double eStepCore( SparseMatrix< int, std::vector<int> >& Xadj, int q, int n, bool stop_at_diagonal=false );
    //! core calculus of E step iteration for empty  Xadj[i]
     double eStepDefault(int q, int i, bool stop_at_diagonal=false );

    //! first m step
    virtual void mFirstStep();
    //! m step
    virtual void mStep();     

    //! prepare m step for vertex i by substraction for old class of i 
    void mStepPrepareVertex(int i);   
    //! m step for vertex i
    virtual void mStepVertex(int i, bool stop_at_diagonal); 
    //! m step and likelihood step
    virtual void mStepMaximization()=0;
    //! m step maximization for Pi
    virtual void mStepMaximizationPi()=0;

    //! computes the complete likelihood
    virtual void mLikelihood();    
    //! likelihood step for vertex i
    void mLikelihoodVertex(int i, bool stop_at_diagonal); 
    //! m step patch at end of em
    void emConcludeAndLikelihood();

    
  public:
    OCEmT( Ermg& ermg, int n_init, int em_nitermax=3 )
      :  EmCore(ermg, em_nitermax){
      _curr_i=0;
      _curr_old_classi=-1;  
      _n_init=n_init;
      _em_nitermax=em_nitermax;
    }
      
      OCEmT(const OCEmT& ocem)
	:  EmCore(ocem){
	_curr_i=ocem._curr_i;
	_curr_old_classi=ocem._curr_old_classi;
	_n_init=ocem._n_init;
	_em_nitermax=ocem._em_nitermax;
      }
	
      virtual ~OCEmT(){};
      
      virtual bool run(int nbclass);      
  };

  class OCEm : public OCEmT
  {
  private:
    //! n_q
    std::vector<int> _n_q;
    //! n_ql
    BandMemMatrix<int> _n_ql;
    //! iterator on _n_ql
    BandMemMatrix<int>::BandCursor _it_n_ql;


    //! E step iteration for vertex i
    virtual void eStepVertex(int i, bool stop_at_diagonal=false);
    //! m step maximization 
    virtual void mStepMaximization();
    //! m step maximization for Pi
    virtual void mStepMaximizationPi();

    virtual void sufficientStatInit(){
      _n_q.assign(_q, 0);
      _n_ql.assign(_q, _q, 0);
      for(int i=0; i<_n_init; i++)    
	nqIncremente( _ermg._class[i] );
    }
    virtual void nqIncremente(int c){
      _n_q[c]++;
    }
    virtual void nqDecremente(int c){
      _n_q[c]--;
    }
    virtual void nqlIncremente(int c1, int c2){
      _n_ql.value(c1, c2)++;
    }
    virtual void nqlDecremente(int c1, int c2){
      _n_ql.value(c1, c2)--;      
    }
  public:
    OCEm( Ermg& ermg, int n_init, int em_nitermax=3 )
      : OCEmT(ermg, n_init, em_nitermax) {}
      
      OCEm(const OCEm& ocem)
	: OCEmT(ocem) {}
	
	virtual ~OCEm(){};
  };

}
#endif
