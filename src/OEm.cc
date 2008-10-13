/* OEm.cc
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


#include <OEm.h>
#include <vector>
#include <algorithm>
#include <cmath>


using namespace std;


namespace ermg {

  bool OEm::run(int nbclass)
  { 
    _q = nbclass;
    _curr_i=_n_init-1;

    init();    
    mFirstStep();
    _curr_i++; 

    bool stop=false;
    while((_curr_i<(_em_nitermax*_n))&&(!stop))
      {         
	eStep();
	mStep();
	_curr_i++;
	if (_curr_i%_n==0)
	  if (this->emConvergency())
	    stop=true;
      }
    mLikelihood();
    return true;
  }
  

  bool OEm::init()
  {
    // UGLY A MODIFIER!!!!!!!!!!!! 
    // A MODIFIER!!!!!!!!!!!!
    _ermg._cardinal_class.assign(_q, 0);
    for(int i=0; i<_n_init; i++)    
      _ermg._cardinal_class[ _ermg._class[i] ]++;
    for(int i=_n_init; i<_n; i++)   
      _ermg._class[i] = -1;
    // A MODIFIER!!!!!!!!!!!!
    // A MODIFIER!!!!!!!!!!!!

    _ermg.emInitTau();
    EmCore::init();
    _t_Tau2 = _t_Tau1;

    _dotTau.assign(_q, 0);
    _Theta.assign(_q, _q, 0); 
    _t_Tau1_back.assign(_q, 0); 

    return true;
  }

  
  void OEm::eStep()
  {
    int i = _curr_i%_n;
    bool real_online = (_curr_i<_n);

    // backup Tau
    if (real_online)
      for (int q=0; q<_ermg._q; q++) 
	_t_Tau1_back[q] = (*_t_Tau1).value(q,i);

    initBetaijql();
    eStepFixedPointVertex( i, real_online );
  }

  void OEm::mFirstStep()
  { 
    for (int i=0; i<_n_init; i++){
      mStepVertex( i, true);// stops at diagonal because symetry
    }
    mStepMaximization();
  }

  void OEm::mStep()
  {   
    mStepVertex( _curr_i%_n, (_curr_i<_n) ); 
    mStepMaximization();
  }

  void OEm::mStepVertex(int i, bool stop_at_diagonal)
  { 
    _it_val1_ql = _val1.access(0);
    _it_val2_ql = _val2.access(0);

    for (int q=0; q<_ermg._q; q++)
      _dotTau[q] += (*_t_Tau1).value(q,i);

    if (_ermg._Xadj[i].size()>0){
      for (int q=0; q<_ermg._q; q++){ 
	_it_Theta_ql = _Theta.access(q);
	_it_denom_Pi_ql = _denom_Pi.access(q);

	double t_Tau1_qi = (*_t_Tau1).value(q,i);
	double b_t_Tau1_qi = _t_Tau1_back[q]; 
	double diff_t_Tau1_qi = t_Tau1_qi-b_t_Tau1_qi; 
	if (!stop_at_diagonal)
	  diff_t_Tau1_qi = t_Tau1_qi;

	// update _dotTau for MStepMaximization
	if (!stop_at_diagonal)
	  _dotTau[q] -= b_t_Tau1_qi;
	_dotTau[q] += t_Tau1_qi;

	for (int l=0; l<_ermg._q; l++){  
	  mStepLikelihoodStepCore( _ermg._Xadj, diff_t_Tau1_qi, i, q, l, stop_at_diagonal );
	  // then considering the fact that we should consider a line and its symetric column
	  if (q>=l){
	    *_it_Theta_ql += _curr_pi;
	    *_it_denom_Pi_ql += _curr_denom;
	    _it_Theta_ql++;      
	    _it_denom_Pi_ql++;
	  }
 	  else{
 	    _Theta.value(l,q) += _curr_pi;
 	    _denom_Pi.value(l,q) += _curr_denom;	    
 	  }
	  _it_val1_ql++;
	  _it_val2_ql++;
	}
      }
    }
    else{
      for (int q=0; q<_ermg._q; q++){
	_it_denom_Pi_ql = _denom_Pi.access(q);
	
	double t_Tau1_qi = (*_t_Tau1).value(q,i);
	double b_t_Tau1_qi = _t_Tau1_back[q];	
	double diff_t_Tau1_qi = t_Tau1_qi-b_t_Tau1_qi; 
	if (!stop_at_diagonal)
	  diff_t_Tau1_qi = t_Tau1_qi;
      
	// update _dotTau for MStepMaximization
	if (!stop_at_diagonal)
	  _dotTau[q] -= b_t_Tau1_qi;
	_dotTau[q] += t_Tau1_qi;

	for (int l=0; l<_ermg._q; l++){
	  mStepLikelihoodStepDefault( diff_t_Tau1_qi, i, q, l, stop_at_diagonal ); 
	  if (q>=l){
	    //_denom_Pi.value(q,l) += _curr_denom;
	    *_it_denom_Pi_ql += _curr_denom;
	    _it_denom_Pi_ql++;
	  }
	  else
	    _denom_Pi.value(l,q) += _curr_denom;

	  _it_val2_ql++;
	}
      }
    }
  }


  void OEm::mStepMaximization()
  {
    bool real_online = (_curr_i<_n);
    double sumalpha = 0.;
    

    // Alpha from dotTau
    for (int l=0; l<_q; l++){ 
      if (real_online)
	_ermg._Alpha[l] = _dotTau[l]/(_curr_i+1);
      else	
	_ermg._Alpha[l] = _dotTau[l]/_n;
  
      if (_ermg._Alpha[l]<MSTEPPARAMMIN)
	_ermg._Alpha[l] = MSTEPPARAMMIN;    
      sumalpha += _ermg._Alpha[l];
    } 

    //! normalization
    for (int l=0; l<_q; l++)
      _ermg._Alpha[l] /= sumalpha;


    // Pi
    mStepMaximizationPi();

  
#ifdef VERBOSE 
    cerr<<"Alpha at end of Mstep"<<endl;
    for (int l=0; l<_q; l++)
      cerr<< _ermg._Alpha[l]<<" ";
    cerr<<endl;

    cerr<<"Pi at end of Mstep"<<endl;
    for (int q=0; q<_q; q++){
      for (int l=0; l<_q; l++) 
	cerr<< _ermg._Pi.value(q,l)<<" ";
      cerr<<endl;
    }
#endif
  }

  void OEm::mStepMaximizationPi()
  {
    // Pi from Theta / denom_Pi
    for (int q=0; q<_q; q++){
      _it_Pi_ql = _ermg._Pi.access(q);
      _it_denom_Pi_ql = _denom_Pi.access(q);
      _it_Theta_ql = _Theta.access(q);
      for (int l=0; l<=q; l++){ 

	if (*_it_denom_Pi_ql>PRECISION)
	  *_it_Pi_ql = *_it_Theta_ql / *_it_denom_Pi_ql;
	else
	  // class with a single vertex
	  *_it_Pi_ql = 0.5;

	if (*_it_Pi_ql<MSTEPPARAMMIN)
	  *_it_Pi_ql = MSTEPPARAMMIN;  
	if (*_it_Pi_ql>(1-MSTEPPARAMMIN))
	  *_it_Pi_ql = 1-MSTEPPARAMMIN;   
    
	_it_Pi_ql++;
	_it_denom_Pi_ql++;
	_it_Theta_ql++;
      }
    }
    
    // symmetrizing
    for (int q=0; q<_q; q++)
      for (int l=q+1; l<_q; l++) 
	_ermg._Pi.value(q,l) = _ermg._Pi.value(l,q);
  }
    
}

