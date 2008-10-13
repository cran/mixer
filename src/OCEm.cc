/* OCEm.cc
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


#include <OCEm.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>


using namespace std;


namespace ermg {
  
  bool OCEmT::run(int nbclass)
  { 
    _q = nbclass;

    if (!this->init())
      return false;
    
    _curr_i=_n_init-1;
    this->mFirstStep();
    _curr_i++; 

    bool stop=false;
    while((_curr_i<(_em_nitermax*_n))&&(!stop))
      {         
	this->eStep();
	this->mStep();
	_curr_i++;
	if (_curr_i%_n==0)
	  if (this->emConvergency())
	    stop=true;	
      }
    this->emConcludeAndLikelihood();

    return true;
  }
  
  bool OCEmT::init()
  {
    EmCore::init();
    sufficientStatInit();
    return true;
  }

  void OCEmT::eStep()
  {
    int i = _curr_i%_n;
    _curr_old_classi = _ermg._class[i];
    bool real_online = (_curr_i<_n);
    initBetaijql();
    eStepVertex( i, real_online );
  }

  void OCEm::eStepVertex(int i, bool stop_at_diagonal)
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
    double maxval = 0.;
    int qstar = 0;
    for (int q=0; q<_q; q++){
      double val = exp( _sumsTau_qi[q]-maxsumTau_qi );
      if (val>maxval){
	qstar = q;
	maxval = val;
      }
    }
    
    // classification
    _ermg._class[i] = qstar;
  }
  
  double OCEmT::eStepCore( SparseMatrix< int, vector<int> >& Xadj, int q, int i, bool stop_at_diagonal )
  {
    int s = Xadj[i].size(), ind1, ind2;
    double sumTau_qi = 0.;
	
    int pos = Xadj.firstPositionAfterDiagonal(i);
    // avoiding Xadj[i][i]
    if (Xadj.fullDiagonal(i)){
      pos--;
    }


    //// LOWER than i
    _it_class = _ermg._class.begin(); 
    _it_Xadj_i = Xadj[i].begin();	     
    ind1=-1;  

    // all the non null elements
    for (int nbe=0; nbe<pos; nbe++){
      ind2 = *_it_Xadj_i;
      // all the null elements between
      for (int j=ind1+1; j<ind2; j++){  
	sumTau_qi += _val2.value( q, *_it_class );
	_it_class++;
      }
      
      // non null element j=ind2
      sumTau_qi += _val1.value( q, *_it_class );
      _it_class++;

      ind1 = ind2;
      _it_Xadj_i++;	  
    }
    // last null elements
    for (int j=ind1+1; j<i; j++){
      sumTau_qi += _val2.value( q, *_it_class );
      _it_class++;
    }

    //// EQUAL to i 
    if (Xadj.fullDiagonal(i)){
      // non null element
      if (_ermg._enable_loop){ 
	sumTau_qi += _val1.value( q, q );
      }
      // zapping Xadj[i][i]
      pos++;
      _it_Xadj_i++;
    }
    else{
      // null element
      if (_ermg._enable_loop){
	sumTau_qi += _val2.value( q, q );
      }
    }
    _it_class++;

    
    if (!stop_at_diagonal){
      //// GREATER than i
      ind1=i;
      // all the non null elements
      for (int nbe=pos; nbe<s; nbe++){
	ind2 = *_it_Xadj_i;
	// all the null elements between
	for (int j=ind1+1; j<ind2; j++){
	  sumTau_qi += _val2.value( q, *_it_class );
	  _it_class++;
	}
      
	// non null element j=ind2
	sumTau_qi += _val1.value( q, *_it_class );
	_it_class++;

	ind1 = ind2;
	_it_Xadj_i++;	  
      }
      // last null elements
      for (int j=ind1+1; j<_n; j++){
	sumTau_qi += _val2.value( q, *_it_class );
	_it_class++;
      }

    }    
    return sumTau_qi;
  }


  double OCEmT::eStepDefault(int q, int i, bool stop_at_diagonal )
  {	
    int ind1;
    double sumTau_qi = 0.;
  
    //// LOWER than i
    _it_class = _ermg._class.begin();  
    ind1=-1;  
    // all null elements
    for (int j=ind1+1; j<i; j++){
      sumTau_qi += _val2.value( q, *_it_class );
      _it_class++;
    }


    //// EQUAL to i
    // null element
    if (_ermg._enable_loop){
      sumTau_qi += _val2.value( q, q );
    }
    _it_class++;


    if (!stop_at_diagonal){
      //// GREATER than i
      ind1=i;  
      // all non null elements
      for (int j=ind1+1; j<_n; j++){
	sumTau_qi += _val2.value( q, *_it_class );
	_it_class++;
      }
    
    }
    return sumTau_qi;
  }


  void OCEmT::mFirstStep()
  { 
    for (int i=0; i<_n_init; i++)  
      this->mStepVertex( i, true );// stops at diagonal because symetry
    this->mStepMaximization();
  }

  void OCEmT::mStep()
  {   
    int i = _curr_i%_n;
    if (_curr_old_classi!=_ermg._class[i]){
      // unclassification
      if (_curr_i>=_n){
	nqDecremente( _curr_old_classi );
	this->mStepPrepareVertex( i );
	this->mStepVertex( i, false ); 
      }
      else
	this->mStepVertex( i, true );
      nqIncremente( _ermg._class[i] );       
      this->mStepMaximization();
    }
  }

  void OCEmT::mStepPrepareVertex(int i)
  {
    int s=_ermg._Xadj[i].size(), ind2;


    if (s>0){

      int pos =_ermg._Xadj.firstPositionAfterDiagonal(i);
      // avoiding _Xadj[i][i]
      if (_ermg._Xadj.fullDiagonal(i)){
	pos--;
      }

      //// LOWER than i
      _it_Xadj_i =  _ermg._Xadj[i].begin();
    
      // all the non null elements
      for (int nbe=0; nbe<pos; nbe++){
	ind2 = *_it_Xadj_i;      
	// all the null elements between
	 
	// non null element
	if (_curr_old_classi>=_ermg._class[ind2])
	  nqlDecremente(_curr_old_classi, _ermg._class[ind2]);
	else
	  nqlDecremente(_ermg._class[ind2], _curr_old_classi);

	_it_Xadj_i++;
      }
      // last null elements
      
      //// EQUAL to i
      if (_ermg._Xadj.fullDiagonal(i)){
	// non null element
	if (_ermg._enable_loop){
	  nqlDecremente(_curr_old_classi, _curr_old_classi);
	}      
	// zapping _Xadj[i][i]
	pos++;
	_it_Xadj_i++;	  
      } 

      //// GREATER than i

      // all the non null elements
      for (int nbe=pos; nbe<s; nbe++){
	ind2 = *_it_Xadj_i;
      
	// all the null elements between

	// non null element
	if (_curr_old_classi>=_ermg._class[ind2])
	  nqlDecremente(_curr_old_classi, _ermg._class[ind2]);
	else
	  nqlDecremente(_ermg._class[ind2], _curr_old_classi);

	_it_Xadj_i++;
      }
      // last null elements
    }
  }  

  void OCEmT::mStepVertex(int i, bool stop_at_diagonal)
  {
    int s=_ermg._Xadj[i].size(), ind2;
    int classi = _ermg._class[i];

    if (s>0){

      int pos =_ermg._Xadj.firstPositionAfterDiagonal(i);
      // avoiding _Xadj[i][i]
      if (_ermg._Xadj.fullDiagonal(i)){
	pos--;
      }

      //// LOWER than i
      _it_Xadj_i =  _ermg._Xadj[i].begin();
    
      // all the non null elements
      for (int nbe=0; nbe<pos; nbe++){
	ind2 = *_it_Xadj_i;      
	// all the null elements between
	 
	// non null element
	if (classi>=_ermg._class[ind2])
	  nqlIncremente(classi, _ermg._class[ind2]);
	else
	  nqlIncremente(_ermg._class[ind2], classi);
	
	_it_Xadj_i++;
      }
      // last null elements
      
      //// EQUAL to i
      if (_ermg._Xadj.fullDiagonal(i)){
	// non null element
	if (_ermg._enable_loop){
	  nqlIncremente(classi, classi);
	}      
	// zapping _Xadj[i][i]
	pos++;
	_it_Xadj_i++;	  
      } 

      if (!stop_at_diagonal){
	//// GREATER than i

	// all the non null elements
	for (int nbe=pos; nbe<s; nbe++){
	  ind2 = *_it_Xadj_i;
      
	  // all the null elements between

	  // non null element
	  if (classi>=_ermg._class[ind2])
	    nqlIncremente(classi, _ermg._class[ind2]);
	  else
	    nqlIncremente(_ermg._class[ind2], classi);
	  
	  _it_Xadj_i++;
	}
	// last null elements
      }
    }
  }  

  void OCEm::mStepMaximization()
  { 
    if (_curr_i<_n)
      for (int l=0; l<_q; l++){
	_ermg._Alpha[l] = _n_q[l]*1.0/(_curr_i+1);
	if (_ermg._Alpha[l]<MSTEPPARAMMIN)
	  _ermg._Alpha[l] = MSTEPPARAMMIN;
      }
    else
      for (int l=0; l<_q; l++){
	_ermg._Alpha[l] = _n_q[l]*1.0/(_n);
	if (_ermg._Alpha[l]<MSTEPPARAMMIN)
	  _ermg._Alpha[l] = MSTEPPARAMMIN;
      }
	

    // Pi
    mStepMaximizationPi();

#ifdef VERBOSE 
    cerr<<"Verif: ";
    int sum=0;
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

  void OCEm::mStepMaximizationPi()
  {
    for (int q=0; q<_q; q++){
      _it_Pi_ql = _ermg._Pi.access(q);
      _it_n_ql = _n_ql.access(q);
      if(_n_q[q]>0){
	for (int l=0; l<q; l++){ 
	  if(_n_q[l]>0){
	    *_it_Pi_ql = double(*_it_n_ql) / double(_n_q[q]*_n_q[l]);
	    if (*_it_Pi_ql<MSTEPPARAMMIN)
	      *_it_Pi_ql = MSTEPPARAMMIN;  
	    if (*_it_Pi_ql>(1-MSTEPPARAMMIN))
	      *_it_Pi_ql = 1-MSTEPPARAMMIN;   
	  }
	  else
	    *_it_Pi_ql = MSTEPPARAMMIN;
	  _it_Pi_ql++;
	  _it_n_ql++;
	}
	if (_ermg._enable_loop)	
	  *_it_Pi_ql =  double(*_it_n_ql) / double(0.5*_n_q[q]*(_n_q[q]-1) + _n_q[q]);
	else	
	  if(_n_q[q]>1)
	    *_it_Pi_ql =  double(*_it_n_ql) / double(0.5*_n_q[q]*(_n_q[q]-1));
	  else
	    *_it_Pi_ql = MSTEPPARAMMIN;	
	if (*_it_Pi_ql<MSTEPPARAMMIN)
	  *_it_Pi_ql = MSTEPPARAMMIN;  
	if (*_it_Pi_ql>(1-MSTEPPARAMMIN))
	  *_it_Pi_ql = 1-MSTEPPARAMMIN;
	
	_it_Pi_ql++;
	_it_n_ql++;
      }
      else{
	for (int l=0; l<_q; l++){ 
	  *_it_Pi_ql = MSTEPPARAMMIN;
	  _it_Pi_ql++;
	  _it_n_ql++;
	}
      }      
    }

    // symmetrizing
    for (int q=0; q<_q; q++)
      for (int l=q+1; l<_q; l++) 
	_ermg._Pi.value(q,l) = _ermg._Pi.value(l,q);
  }


  void OCEmT::mLikelihood()
  {
    this->initBetaijql();
    _ermg.resetLikelihood();
    
    for (int i=0; i<_n; i++){
      this->mLikelihoodVertex(i, true);
    }
  }

  void OCEmT::emConcludeAndLikelihood()
  {
    this->mLikelihood();
    for (int i=0; i<_n; i++){
      (*_t_Tau1).value(_ermg._class[i], i) = 1;
    }
  } 

 
 
  void OCEmT::mLikelihoodVertex(int i, bool stop_at_diagonal)
  {
    // almost same code than eStepCore but with a different treatment 
    // for the "equal to i" part

    int s=_ermg._Xadj[i].size(), ind1, ind2;

    _it_class = _ermg._class.begin();
    int classi = _ermg._class[i];

    if (s>0){
      double tmplikelihood = 0.;

      int pos =_ermg._Xadj.firstPositionAfterDiagonal(i);
      // avoiding _Xadj[i][i]
      if (_ermg._Xadj.fullDiagonal(i)){
	pos--;
      }

      //// LOWER than i
      _it_class = _ermg._class.begin(); 
      _it_Xadj_i =  _ermg._Xadj[i].begin();	     
      ind1=-1;  
    
      // all the non null elements
      for (int nbe=0; nbe<pos; nbe++){
	ind2 = *_it_Xadj_i;
      
	// all the null elements between
	for (int j=ind1+1; j<ind2; j++){
	  tmplikelihood += _val2.value(classi,*_it_class );
	  _it_class++;
	}
	// non null element
	tmplikelihood += _val1.value(classi, *_it_class );
	_it_class++;

	ind1 = ind2;	
	_it_Xadj_i++;
      }
      // last null elements
      for (int j=ind1+1; j<i; j++){
	tmplikelihood += _val2.value(classi, *_it_class );
	_it_class++;
      }

      //// EQUAL to i
      if (_ermg._Xadj.fullDiagonal(i)){
	// non null element
	if (_ermg._enable_loop){
	  // direct add, unvariant if or not directed
	  _ermg.addToCompleteLikelihood(_val1.value(classi, classi));
	}      
	// zapping _Xadj[i][i]
	pos++;
	_it_Xadj_i++;	  
      }
      else{
	// null element
	if (_ermg._enable_loop){
	  // direct add, unvariant if or not directed
	  _ermg.addToCompleteLikelihood(_val2.value(classi, classi));
	}
      } 
      _it_class++;     

      if (!stop_at_diagonal){
	//// GREATER than i 
	ind1=i; 

	// all the non null elements
	for (int nbe=pos; nbe<s; nbe++){
	  ind2 = *_it_Xadj_i;
      
	  // all the null elements between
	  for (int j=ind1+1; j<ind2; j++){
	    tmplikelihood += _val2.value(classi, *_it_class);
	    _it_class++;
	  }
	  // non null element
	  tmplikelihood += _val1.value(classi, *_it_class);
	  _it_class++;

	  ind1 = ind2;	
	  _it_Xadj_i++;
	}
	// last null elements
	for (int j=ind1+1; j<_ermg._n; j++){
	  tmplikelihood += _val2.value(classi, *_it_class);
	  _it_class++;
	}
      }
      _ermg.addToCompleteLikelihood( tmplikelihood );
      _ermg.addToCompleteLikelihood( log(_ermg._Alpha[classi]) );
    }
    else{
      double tmplikelihood = 0.;
   
      //// LOWER than i 
      ind1=-1;
      // all null elements
      for (int j=ind1+1; j<i; j++){
	tmplikelihood += _val2.value(classi, *_it_class);
	_it_class++;
      }

      //// EQUAL to i
      // null element
      if (_ermg._enable_loop){
	// direct add, unvariant if or not directed
	_ermg.addToCompleteLikelihood(_val2.value(classi, *_it_class));
      }    
      _it_class++;

      if (!stop_at_diagonal){
	//// GREATER than i 
	ind1=i;
	// all null elements
	for (int j=ind1+1; j<_ermg._n; j++){
	  tmplikelihood += _val2.value(classi, *_it_class);
	  _it_class++;
	}
	
	_ermg.addToCompleteLikelihood( tmplikelihood );	
	_ermg.addToCompleteLikelihood( log(_ermg._Alpha[classi]) );
      }
    }
  }

}

