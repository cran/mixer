/* SOCEm.h
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
 *       \brief Stochastic Online Classification EM algorithm
 *
 *          10/01/2007
 *          
*/

#ifndef ERMG_SOCEM_H
#define ERMG_SOCEM_H
/*!
  \class SOCEm libermg/SOCEm.h
  \brief Stochastic Online Classification EM algorithm
*/

#include <OCEm.h>


namespace ermg {

  class SOCEm : public OCEmT
  {
  protected:
    //! n_q
    std::vector<double> _n_q;
    //! n_ql
    BandMemMatrix<double> _n_ql;
    //! iterator on _n_ql
    BandMemMatrix<double>::BandCursor _it_n_ql;

    BandMemMatrix<double> _m_ql;
    //! iterator on _m_ql
    BandMemMatrix<double>::BandCursor _it_m_ql;

    //! ponderate of the updated values
    double _ponderate;
    void initPonderate(){
      _ponderate = 1.;
    }
    void updatePonderate(){
      if (_curr_i<_n)
	_ponderate = 1.;
      else
	_ponderate = 1./( 1 + pow(double(_curr_i)/_n, (2*log(10))/log(double(_em_nitermax))) );     
    }	
    
    //! e step iteration for vertex i
    virtual void eStepVertex(int i, bool stop_at_diagonal=false);

  
    //! m step for vertex i
    virtual void mStepVertex(int i, bool stop_at_diagonal); 
    //! m step m_ql computation
    void mStepMqlUpdate(int i, bool stop_at_diagonal=false);
    //! m step maximization 
    virtual void mStepMaximization();
    //! m step maximization for Pi
    virtual void mStepMaximizationPi();

    virtual void sufficientStatInit(){
      _n_q.assign(_q, 0.);
      _n_ql.assign(_q, _q, 0.);
      initPonderate();
      for(int i=0; i<_n_init; i++)    
	nqIncremente( _ermg._class[i] );
      _m_ql.assign(_q, _q, 0.);
    }
    virtual void nqIncremente(int c){
      if(_curr_i<=_n)
	_n_q[c]+=1.; 
      else
	_n_q[c]+=_ponderate;	
    }
    virtual void nqDecremente(int c){
      if(_curr_i<=_n)
      _n_q[c]-=1.;
      else
	_n_q[c]-=_ponderate;
    }
    virtual void nqlIncremente(int c1, int c2){
      if(_curr_i<=_n)
	_n_ql.value(c1, c2)+=1.;
      else
      _n_ql.value(c1, c2)+=_ponderate;   
    }
    virtual void nqlDecremente(int c1, int c2){
      if(_curr_i<=_n)
	_n_ql.value(c1, c2)-=1.;
      else  
	_n_ql.value(c1, c2)-=_ponderate;   
    }
    
  public:
    SOCEm(Ermg& ermg, int n_init, int em_nitermax=3)
      : OCEmT(ermg, n_init, em_nitermax) {
      _ponderate=1.;
    }
      
      SOCEm(const SOCEm& socem)
	: OCEmT(socem) {
	_ponderate=socem._ponderate;
      }
	
	virtual ~SOCEm(){};	
  };
  
}
#endif
