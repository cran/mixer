/* EmCore.cc
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

#include <EmCore.h>

using namespace std;


namespace ermg {

 bool EmCore::emConvergency()
  {
    for (int q=0; q<_ermg._q; q++)
      if ((_ermg._Alpha[q]-_Alpha_back[q])/_ermg._Alpha[q] > _em_eps){
	_Alpha_back = _ermg._Alpha;
	_Pi_back = _ermg._Pi;
	return false;
      }

    _it_Pi_ql = _ermg._Pi.access(0);
    _it_Pi_back_ql = _Pi_back.access(0);
    for (int q=0; q<_ermg._q; q++)
      for (int l=0; l<_ermg._q; l++)
	if ((*_it_Pi_ql-*_it_Pi_back_ql)/(*_it_Pi_ql) > _em_eps){
	  _Alpha_back = _ermg._Alpha;
	  _Pi_back = _ermg._Pi;
	  //       if ((_Pi.value(q,l)-_Pi_back.value(q,l))/_Pi.value(q,l) > _em_eps){
	  // 	_Alpha_back = _Alpha;
	  // 	_Pi_back = _Pi;
	  _it_Pi_ql++;
	  _it_Pi_back_ql++;
	  return false;
	}

    return true;
  } 

  bool EmCore::init( )
  {
    _Alpha_back.assign(_ermg._n, 0.);
    _Pi_back.assign(_ermg._q, _ermg._q, 0.);
    _denom_Pi.assign(_ermg._q, _ermg._q, 0.);
    
    _sumsTau_qi.assign( _ermg._n, 0);
    _val1.assign(_ermg._q, _ermg._q, 0.);
    _val2.assign(_ermg._q, _ermg._q, 0.);

    _t_Tau1 = _ermg._t_Tau;
    if (_t_Tau1)
      return true;
    else 
      return false;
  }
  

  void EmCore::initBetaijql()
  {  
    _it_Pi_ql = _ermg._Pi.access(0);
    _it_val1_ql = _val1.access(0);
    _it_val2_ql = _val2.access(0);

    for (int q=0; q<_ermg._q; q++){ 
      for (int l=0; l<_ermg._q; l++){
	if ( (*_it_Pi_ql < PRECISION)|| 
	     ((*_it_Pi_ql > 1-PRECISION)&&(*_it_Pi_ql < 1+PRECISION)) ){
	  *_it_val1_ql = 0;
	  *_it_val2_ql = 0;
	}
	else{  
	  *_it_val1_ql = log(*_it_Pi_ql);
	  *_it_val2_ql = log(1-(*_it_Pi_ql)); 
	}
	_it_Pi_ql++;
	_it_val1_ql++; 
	_it_val2_ql++;
      }
    }

  }

  void EmCore::initPiDenomPi()
  {
    _it_Pi_ql = _ermg._Pi.access(0);
    _it_denom_Pi_ql = _denom_Pi.access(0);
    for (int q=0; q<_ermg._q; q++){
      for (int l=0; l<_ermg._q; l++){
	*_it_Pi_ql = 0;
	*_it_denom_Pi_ql = 0;
	_it_Pi_ql++;
	_it_denom_Pi_ql++;
      }
    }
  }
}
