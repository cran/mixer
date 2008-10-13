/* SOCEm.cc
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


#include <limits>
#include <SOCEm.h>

using namespace std;
namespace ermg {
  
  void SOCEm::eStepVertex(int i, bool stop_at_diagonal)
  { 
    double maxsumTau_qi = -numeric_limits<double>::max();
      
    if (_ermg._Xadj[i].size()>0){
      for (int q=0; q<_q; q++){
	_sumsTau_qi[q] = eStepCore( _ermg._Xadj, q, i, stop_at_diagonal) + log(_ermg._Alpha[q]);
	if (_sumsTau_qi[q] > maxsumTau_qi)
	  maxsumTau_qi = _sumsTau_qi[q];
      }
    }
    else{
      for (int q=0; q<_q; q++){	  
	_sumsTau_qi[q] = eStepDefault( q, i, stop_at_diagonal ) + log(_ermg._Alpha[q]);
	if (_sumsTau_qi[q] > maxsumTau_qi)
	  maxsumTau_qi = _sumsTau_qi[q];
      }
    }      


    // factorization
    double sumval=0; 
    double max=-10000000;
    for (int q=0; q<_q; q++){
      double val = exp( _sumsTau_qi[q]-maxsumTau_qi );
      sumval += val;
      _sumsTau_qi[q] = val;
      if (val>max){
	max=val;
      }
    }
      
    // drawing    
    int qstar = 0;
    double c=_sumsTau_qi[0], tmprand = double(rand())/double(RAND_MAX), bound=tmprand*sumval;
    while (c<bound){
      qstar++;
      c+=_sumsTau_qi[qstar];
    }
	  

    // classification
    _ermg._class[i] = qstar;
  }
   
  void SOCEm::mStepVertex(int i, bool stop_at_diagonal)
  { 
    OCEmT::mStepVertex(i, stop_at_diagonal);
    updatePonderate();
    mStepMqlUpdate(i, stop_at_diagonal);
  }


  void SOCEm::mStepMqlUpdate(int i, bool stop_at_diagonal)    
  {  
    int qstar = _ermg._class[i];
    _it_class=_ermg._class.begin();

    // real online
    if (stop_at_diagonal){
      for (int j=0; j<i; j++){
	if (qstar>=*_it_class){
	  _m_ql.value(qstar,*_it_class)++;
	}
	else{
	  _m_ql.value(*_it_class, qstar)++;
	}
 	_it_class++;
      }
      if (_ermg._enable_loop)
 	_m_ql.value(qstar,qstar)++;	  
    }
    else{
      // after real online
      for (int j=0; j<i; j++){
	if (_curr_old_classi>=*_it_class)
	  _m_ql.value(_curr_old_classi, *_it_class)-=_ponderate;
	else
	  _m_ql.value(*_it_class, _curr_old_classi)-=_ponderate;

	if (qstar>=*_it_class)
	  _m_ql.value(qstar, *_it_class)+=_ponderate;	  
	else
	  _m_ql.value(*_it_class, qstar)+=_ponderate;

	_it_class++;
      }
      if (_ermg._enable_loop){	    	    
	_m_ql.value(_curr_old_classi, _curr_old_classi)-=_ponderate;
	_m_ql.value(qstar,qstar)+=_ponderate;
      }
      _it_class++;
      for (int j=i+1; j<_n; j++){
	if (_curr_old_classi>=*_it_class){
	  //cerr<<i<<","<<j<<" -> "<<_curr_old_classi<<" "<<*_it_class<<endl;	    	    
	  _m_ql.value(_curr_old_classi, *_it_class)-=_ponderate;
	}
	else	    	    
	  _m_ql.value(*_it_class, _curr_old_classi)-=_ponderate;
	
	if (qstar>=*_it_class){
	  //cerr<<i<<","<<j<<" -> "<<_curr_old_classi<<" "<<*_it_class<<endl;	    	    
	  _m_ql.value(qstar, *_it_class)+=_ponderate;
	}
	else	    	    
	  _m_ql.value(*_it_class, qstar)+=_ponderate;
       
	_it_class++;
      }	 		  	
    }
  }

  void SOCEm::mStepMaximization()
  { 
    if (_curr_i<_n)
      for (int l=0; l<_q; l++){
	_ermg._Alpha[l] = _n_q[l]/(_curr_i+1);
	if (_ermg._Alpha[l]<MSTEPPARAMMIN)
	  _ermg._Alpha[l] = MSTEPPARAMMIN;
      }
    else
      for (int l=0; l<_q; l++){
	_ermg._Alpha[l] = _n_q[l]/(_n);
	if (_ermg._Alpha[l]<MSTEPPARAMMIN)
	  _ermg._Alpha[l] = MSTEPPARAMMIN;
      }
	
    // Pi
    mStepMaximizationPi();


#ifdef VERBOSE 
    cerr<<"Verif: ";
    double sum=0.;
    for (int q=0; q<_q; q++){
      for (int l=0; l<_q; l++){ 
        sum+=_n_ql.value(q,l);
      }
    }
    cerr<<" nbedges: "<<sum;
    sum=0;
    for (int q=0; q<_q; q++){
      sum+=_n_q[q];
    }
    cerr<<" nbvertices: "<<sum<<endl;

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

  void SOCEm::mStepMaximizationPi()
  {    
    for (int q=0; q<_q; q++){
      _it_Pi_ql = _ermg._Pi.access(q);
      _it_n_ql = _n_ql.access(q);
      _it_m_ql = _m_ql.access(q);
      for (int l=0; l<=q; l++){
	if(*_it_m_ql>0){
  	  *_it_Pi_ql = *_it_n_ql / *_it_m_ql;
	    
  	  if (*_it_Pi_ql<MSTEPPARAMMIN)
  	    *_it_Pi_ql = MSTEPPARAMMIN;  
  	  if (*_it_Pi_ql>(1-MSTEPPARAMMIN))
  	    *_it_Pi_ql = 1-MSTEPPARAMMIN;   
  	}
  	else
  	  *_it_Pi_ql = MSTEPPARAMMIN;
 	_it_Pi_ql++;
 	_it_n_ql++;
 	_it_m_ql++;
      }      
    }
    

    // symmetrizing
    for (int q=0; q<_q; q++)
      for (int l=q+1; l<_q; l++) 
	_ermg._Pi.value(q,l) = _ermg._Pi.value(l,q);
  }
 
}
