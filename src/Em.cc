/* Em.cc
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


#include <Em.h>
#include <vector>
#include <algorithm>
#include <limits>


using namespace std;


namespace ermg {

  Em::Em( Ermg& ermg,
	  int em_nitermax, double em_eps,
	  int e_step_fixed_point_nitermax, double e_step_fixed_point_eps, bool classif  )     
    : EmCore(ermg, em_nitermax, em_eps)
  {
    _e_step_fixed_point_nitermax = e_step_fixed_point_nitermax;
    _e_step_fixed_point_eps = e_step_fixed_point_eps;
    _e_step_fixed_point_convergency = false;
    _aux_t_Tau = NULL;    
    _t_Tau2 = NULL;
    _classif = classif;
  }

  Em::~Em()
  {
    if (_aux_t_Tau){
      delete _aux_t_Tau;
      _aux_t_Tau=NULL;
    }
  }

  bool Em::run(int nbclass)
  {  
    _q = nbclass;

    int niter = 0;
    if (!this->init())
      return false;
  
    this->mStep();
    do{
    
      this->eStep();
      this->mStep();
      niter++;
    
    }
    while( (!this->emConvergency())&&(niter<_em_nitermax) );

    this->finish();
    return true;
  }


  bool Em::init( )
  {
    _ermg.emInitTau();
    
    this->EmCore::init();
    this->newAuxTau();

    return true;
  }
  
  void Em::finish( )
  {
    _ermg._t_Tau =_t_Tau1;
    _aux_t_Tau =_t_Tau2;  
    _t_Tau1 = NULL;  
    _t_Tau2 = NULL;
  }

  void Em::newAuxTau()
  {
    if (_aux_t_Tau==NULL)
      _aux_t_Tau = new BandMemMatrix<double>(_ermg._q, _ermg._n);
    else
      if (((_aux_t_Tau)->nLines()!=_q)||((_aux_t_Tau)->nCols()!=_n)){
	delete _aux_t_Tau;
	_aux_t_Tau = new BandMemMatrix<double>(_ermg._q, _ermg._n);
      }

    _t_Tau2 = _aux_t_Tau;
  }


  void Em::eStep()
  {
    int niter = 0;   
    _e_step_fixed_point_convergency = false;

    while( (!_e_step_fixed_point_convergency)&&(niter<_e_step_fixed_point_nitermax) ){
      this->eStepFixedPoint();
      this->eStepFixedPointSwitchTau2toTau1();
      niter++;
    }  
#ifdef VERBOSE
    cerr<<"Estep end after "<<niter<<" fixed point iterations"<<endl;
#endif
  }

  void Em::eStepFixedPoint()
  {
    this->initBetaijql();

    _e_step_fixed_point_convergency = true;
    for (int i=0; i<_ermg._n; i++){
      this->eStepFixedPointVertex( i );
    }
  }

  void Em::eStepFixedPointSwitchTau2toTau1()
  {  
    BandMemMatrix<double>* ptrbmm = _t_Tau1;
    _t_Tau1 = _t_Tau2;
    _t_Tau2 = ptrbmm;
  }
 

  void Em::eStepFixedPointVertex(int i, bool stop_at_diagonal)
  { 
    double sum_t_Tau2_qi = 0;
    double maxsumTau_qi = -numeric_limits<double>::max();

    _it_val1_ql = _val1.access(0);
    _it_val2_ql = _val2.access(0);

    int qmax=0;

    if (_ermg._Xadj[i].size()>0){
      for (int q=0; q<_ermg._q; q++){
	double sumTau_qi = 0.;
	for (int l=0; l<_ermg._q; l++){ 
	  sumTau_qi += eStepFixedPointStepCore( _ermg._Xadj, *_it_val1_ql, *_it_val2_ql, 
						i, l, l==q, stop_at_diagonal );
	  _it_val1_ql++;
	  _it_val2_ql++;
      
	}
	_sumsTau_qi[q] = sumTau_qi + log(_ermg._Alpha[q]);
	if (_sumsTau_qi[q] > maxsumTau_qi){
	  maxsumTau_qi = _sumsTau_qi[q];
	  qmax = q;
	}
      }
    }
    else{
      for (int q=0; q<_ermg._q; q++){	
	double sumTau_qi = 0.;      
	for (int l=0; l<_ermg._q; l++){ 	  
	  sumTau_qi += eStepFixedPointStepDefault( *_it_val2_ql,
						   i, l, l==q, stop_at_diagonal );
	  _it_val2_ql++;      
	}
	_sumsTau_qi[q] = sumTau_qi + log(_ermg._Alpha[q]);
	if (_sumsTau_qi[q] > maxsumTau_qi){
	  maxsumTau_qi = _sumsTau_qi[q];
	  qmax = q;
	}
      }
    }

    if (!_classif){
      // factorization
      for (int q=0; q<_ermg._q; q++){
	double val = exp( _sumsTau_qi[q]-maxsumTau_qi );
	(*_t_Tau2).value(q,i) = val;
	sum_t_Tau2_qi += val;
      }
    }
    else{
      // MAP
      for (int q=0; q<_ermg._q; q++)
       if (q==qmax)
         (*_t_Tau2).value(q,i) = 1;
       else
         (*_t_Tau2).value(q,i) = 0;
      sum_t_Tau2_qi = 1;

    }
    this->eStepFixedPointStepNormalizeTauAndConvergency(i, sum_t_Tau2_qi);  
  } 


  double Em::eStepFixedPointStepCore( SparseMatrix< int, vector<int> >& Xadj, 
				      double val1, double val2,
				      int i, int l, bool l_is_q,
				      bool stop_at_diagonal )
  {
    int s = Xadj[i].size(), ind1, ind2;
    double sumTau_qi = 0.;
	
    int pos = Xadj.firstPositionAfterDiagonal(i);
    // avoiding Xadj[i][i]
    if (Xadj.fullDiagonal(i)){
      pos--;
    }

    //// LOWER than i (Tau2 concerned) 
    _it_Xadj_i = Xadj[i].begin();	     
    ind1=-1;  

#ifdef ACCELERATE
    _it_t_Tau2_lj = (*_t_Tau2).access(l, ind1+1);
#else
    _it_t_Tau2_lj = (*_t_Tau1).access(l, ind1+1);
#endif
    // all the non null elements
    for (int nbe=0; nbe<pos; nbe++){
      ind2 = *_it_Xadj_i;
      // all the null elements between
      for (int j=ind1+1; j<ind2; j++){
	sumTau_qi += val2*(*_it_t_Tau2_lj);

	_it_t_Tau2_lj++;
      }
      
      // non null element j=ind2
      sumTau_qi += val1*(*_it_t_Tau2_lj);

      _it_t_Tau2_lj++;
      ind1 = ind2;
      _it_Xadj_i++;	  
    }
    // last null elements
    for (int j=ind1+1; j<i; j++){
      sumTau_qi += val2*(*_it_t_Tau2_lj);

      _it_t_Tau2_lj++;
    }


    //// EQUAL to i	
    if (Xadj.fullDiagonal(i)){
      if (l_is_q){
	// non null element
	if (_ermg._enable_loop){ 
	  sumTau_qi += val1;
	}
      }
      // zapping Xadj[i][i]
      pos++;
      _it_Xadj_i++;
    }
    else{
      if (l_is_q){
	// null element
	if (_ermg._enable_loop){
	  sumTau_qi += val2;
	}
      }
    }

    if (!stop_at_diagonal){
      
      //// GREATER than i (Tau1 concerned)
      ind1=i;
      _it_t_Tau1_lj = (*_t_Tau1).access(l, ind1+1);
      // all the non null elements
      for (int nbe=pos; nbe<s; nbe++){
	ind2 = *_it_Xadj_i;
	// all the null elements between
	for (int j=ind1+1; j<ind2; j++){
	  sumTau_qi += val2*(*_it_t_Tau1_lj);
	  
	  _it_t_Tau1_lj++;
	}
	
	// non null element j=ind2
	sumTau_qi += val1*(*_it_t_Tau1_lj);
	_it_t_Tau1_lj++;
	ind1 = ind2;
	_it_Xadj_i++;	  
      }
      // last null elements
      for (int j=ind1+1; j<_ermg._n; j++){
	sumTau_qi += val2*(*_it_t_Tau1_lj);
	
	_it_t_Tau1_lj++;
      }
  
    }
    return sumTau_qi;
  }


  double  Em::eStepFixedPointStepDefault( double val2, 
					  int i, int l, bool l_is_q,
					  bool stop_at_diagonal )
  {	
    int ind1;
    double sumTau_qi = 0.;
  
    //// LOWER than i (Tau2 concerned)  
    ind1=-1;
#ifdef ACCELERATE
    _it_t_Tau2_lj = (*_t_Tau2).access(l, ind1+1);
#else
    _it_t_Tau2_lj = (*_t_Tau1).access(l, ind1+1);
#endif
    // all non null elements
    for (int j=ind1+1; j<i; j++){
      sumTau_qi += val2*(*_it_t_Tau2_lj);
      _it_t_Tau2_lj++;
    }


    //// EQUAL to i 
    if (l_is_q){
      // null element
      if (_ermg._enable_loop){
	sumTau_qi += val2;
      }
    }


    if (!stop_at_diagonal){
      
      //// GREATER than i (Tau1 concerned)
      ind1=i;  
      _it_t_Tau1_lj = (*_t_Tau1).access(l, ind1+1);
      // all non null elements
      for (int j=ind1+1; j<_ermg._n; j++){
	sumTau_qi += val2*(*_it_t_Tau1_lj); 
	_it_t_Tau1_lj++;
      }
      
    }
    return sumTau_qi;
  }


  void Em::eStepFixedPointStepNormalizeTauAndConvergency(int i, double sum_t_Tau2_qi)
  {
    double tmp;
    double tmp_sum_t_Tau2_qi = 0.;
    bool tonormalize = false;
    // pre-normalization, checking low Tau
    for (int q=0; q<_ermg._q; q++){
      double& tau2_qi_ref = (*_t_Tau2).value(q,i);
      tmp = tau2_qi_ref / sum_t_Tau2_qi;
      if (tmp<TAUMIN){
	tau2_qi_ref = TAUMIN;
	tmp_sum_t_Tau2_qi += TAUMIN;
	tonormalize = true; 
      }
      else{
	tau2_qi_ref = tmp;
	tmp_sum_t_Tau2_qi += tmp;
      }
    }
    // normalization 
    if (tonormalize)
      for (int q=0; q<_ermg._q; q++){
	double& tau2_qi_ref = (*_t_Tau2).value(q,i);
	tau2_qi_ref = tau2_qi_ref / tmp_sum_t_Tau2_qi;
      }

    // Convergency
    for (int q=0; q<_ermg._q; q++){
      if ( fabs((*_t_Tau2).value(q,i)-(*_t_Tau1).value(q,i))> _e_step_fixed_point_eps ){
	_e_step_fixed_point_convergency = false;
      }
    }
  }


  void Em::mFirstStep()
  {
    int back_n = _n;
    _n = _n_init;
    for (int i=0; i<_n_init; i++){
      mStepLikelihoodVertex( i, true);// stops at diagonal because symetry
    }
    mStepMaximization();
    _n = back_n;
  }
  
  void Em::mStep()
  { 
    this->initPiDenomPi();
    _ermg.resetLikelihood();
    for (int i=0; i<_ermg._n; i++){
      this->mStepLikelihoodVertex(i);
    }
    this->mStepMaximization();
  }

  void Em::mLikelihood()
  { 
    _ermg.resetLikelihood();
    for (int i=0; i<_ermg._n; i++){
      this->mStepLikelihoodVertex(i, false);
    }
  }

  void Em::mStepLikelihoodVertex(int i, bool toestim)
  {	
    _it_val1_ql = _val1.access(0);
    _it_val2_ql = _val2.access(0);
    
    for (int q=0; q<_ermg._q; q++){
      if (toestim){
	_it_Pi_ql = _ermg._Pi.access(q);
	_it_denom_Pi_ql = _denom_Pi.access(q);
      }    
      double t_Tau1_qi = (*_t_Tau1).value(q,i);
      if (_ermg._Xadj[i].size()>0)
	_ermg.addToCompleteLikelihood( t_Tau1_qi*log(_ermg._Alpha[q])
				      + mStepLikelihoodVertexRun(t_Tau1_qi, i, q, toestim) );
      else
	_ermg.addToCompleteLikelihood( t_Tau1_qi*log(_ermg._Alpha[q])
				       + mStepLikelihoodVertexRunDefault(t_Tau1_qi, i, q, toestim) );      
      if (t_Tau1_qi>0.)
	_ermg._entropy += t_Tau1_qi*log(t_Tau1_qi);
    } 
  }
 

  double Em::mStepLikelihoodVertexRun(double t_Tau1_qi, int i, int q, bool toestim)
  {
    double tmplikelihood = 0.;
    for (int l=0; l<_ermg._q; l++){      
      tmplikelihood += mStepLikelihoodStepCore( _ermg._Xadj, t_Tau1_qi, i, q, l, true );
      if (toestim){
	if (q>=l){
	  *_it_Pi_ql += _curr_pi;
	  *_it_denom_Pi_ql += _curr_denom;
	  _it_Pi_ql++;      
	  _it_denom_Pi_ql++;
	}
	else{
	  _ermg._Pi.value(l,q) += _curr_pi;
	  _denom_Pi.value(l,q) += _curr_denom;
	}
      }
      _it_val1_ql++;
      _it_val2_ql++;
    }
    return tmplikelihood;
  }

  double Em::mStepLikelihoodVertexRunDefault(double t_Tau1_qi, int i, int q, bool toestim)
  { 
    double tmplikelihood = 0.;
    for (int l=0; l<_ermg._q; l++){
      tmplikelihood += mStepLikelihoodStepDefault( t_Tau1_qi, i, q, l, true );
      if (toestim){
	if (q>=l){	  
	  *_it_denom_Pi_ql += _curr_denom;
	  _it_denom_Pi_ql++;
	}
	else
	  _denom_Pi.value(l,q) += _curr_denom;
      }
      _it_val2_ql++;	    
    }
    return tmplikelihood;
  }

  
          
  double Em::mStepLikelihoodStepCore( SparseMatrix< int, std::vector<int> >& Xadj, 
				      double t_Tau1_qi, int i, int q, int l, bool stop_at_diagonal )    
  {
    double c;
    int s=_ermg._Xadj[i].size(), ind1, ind2;
    double tmplikelihood=0;

    // initialize the beta_ijql as 2 possible values      
    double val1 = *_it_val1_ql;
    double val2 = *_it_val2_ql;

    _curr_pi=0;
    _curr_denom=0;
    
    int pos =_ermg._Xadj.firstPositionAfterDiagonal(i);
    // avoiding _Xadj[i][i]
    if (_ermg._Xadj.fullDiagonal(i)){
      pos--;
    }
    
    
    //// LOWER than i
    _it_Xadj_i =  _ermg._Xadj[i].begin();	     
    ind1=-1;  
    
    _it_t_Tau1_lj = (*_t_Tau1).access(l, ind1+1);
    // all the non null elements
    for (int nbe=0; nbe<pos; nbe++){
      ind2 = *_it_Xadj_i;
      
      // all the null elements between
      for (int j=ind1+1; j<ind2; j++){
	c = t_Tau1_qi * (*_it_t_Tau1_lj);
	tmplikelihood += c*val2;
	_curr_denom += c;
	_it_t_Tau1_lj++;
      }
      // non null element
      c = t_Tau1_qi*(*_it_t_Tau1_lj);
      tmplikelihood += c*val1;
      _curr_pi += c;
      _curr_denom += c;
      ind1 = ind2;
      _it_t_Tau1_lj++;	
      _it_Xadj_i++;
    }
    // last null elements
    for (int j=ind1+1; j<i; j++){
      c = t_Tau1_qi*(*_it_t_Tau1_lj);
      tmplikelihood += c*val2;
      _curr_denom += c;
      _it_t_Tau1_lj++;
    }


    //// EQUAL to i
    bool l_is_q = (l==q);	
    if (_ermg._Xadj.fullDiagonal(i)){
      if (l_is_q){
	// non null element
	if (_ermg._enable_loop){
	  c = t_Tau1_qi;
	  // direct add, unvariant if or not directed
	  _ermg.addToCompleteLikelihood(c*val1); 
	  _curr_pi += c;
	  _curr_denom += c;
	}
      }
      // zapping _Xadj[i][i]
      pos++;
      _it_Xadj_i++;	  
    }
    else{
      if (l_is_q){
	// null element
	if (_ermg._enable_loop){
	  c = t_Tau1_qi;
	  // direct add, unvariant if or not directed
	  _ermg.addToCompleteLikelihood(c*val2);
	  _curr_denom += c;
	}
      }
    }
    _it_t_Tau1_lj++;

    if (!stop_at_diagonal){
      //// GREATER than i 
      ind1=i; 
      // all the non null elements
      for (int nbe=pos; nbe<s; nbe++){
	ind2 = *_it_Xadj_i;
      
	// all the null elements between
	for (int j=ind1+1; j<ind2; j++){
	  c = t_Tau1_qi * (*_it_t_Tau1_lj);
	  tmplikelihood += c*val2;
	  _curr_denom += c;
	  _it_t_Tau1_lj++;
	}
	// non null element
	c = t_Tau1_qi*(*_it_t_Tau1_lj);
	tmplikelihood += c*val1;
	_curr_pi += c;
	_curr_denom += c;
	ind1 = ind2;
	_it_t_Tau1_lj++;	
	_it_Xadj_i++;
      }
      // last null elements
      for (int j=ind1+1; j<_ermg._n; j++){
	c = t_Tau1_qi*(*_it_t_Tau1_lj);
	tmplikelihood += c*val2;
	_curr_denom += c;
	_it_t_Tau1_lj++;
      }
    }
    return tmplikelihood;
  }


  double Em::mStepLikelihoodStepDefault( double t_Tau1_qi, int i, int q, int l, bool stop_at_diagonal )
  {
    double c;
    int ind1;
    double tmplikelihood=0;

    // initialize the only possible beta_ijql value
    double val2 = *_it_val2_ql;

    _curr_denom=0;


    //// LOWER than i 
    ind1=-1;
    _it_t_Tau1_lj = (*_t_Tau1).access(l, ind1+1);
    // all null elements
    for (int j=ind1+1; j<i; j++){
      c = t_Tau1_qi * (*_it_t_Tau1_lj);
      tmplikelihood += c*val2;
      _curr_denom += c;
      _it_t_Tau1_lj++;
    }


    //// EQUAL to i
    bool l_is_q = (l==q); 
    if (l_is_q){
      // null element
      if (_ermg._enable_loop){
	c = t_Tau1_qi;
	// direct add, unvariant if or not directed
	_ermg.addToCompleteLikelihood(c*val2);
	_curr_denom += c;
      }
    }    
    _it_t_Tau1_lj++;


    if (!stop_at_diagonal){
      //// GREATER than i 
      ind1=i;
      // all null elements
      for (int j=ind1+1; j<_ermg._n; j++){
	c = t_Tau1_qi * (*_it_t_Tau1_lj);
	tmplikelihood += c*val2;
	_curr_denom += c;
	_it_t_Tau1_lj++;
      }
    }
    return tmplikelihood;
  }

  void Em::mStepMaximization()
  {
    double sumalpha = 0.;
    _it_t_Tau1_lj = (*_t_Tau1).access(0);
    for (int l=0; l<_q; l++){ 
      double sumT_lj = 0.;
      for (int j=0; j<_n; j++){    
	//sumT_lj += (*_t_Tau1).value(l,j);
	sumT_lj += *_it_t_Tau1_lj;
	_it_t_Tau1_lj++;  
      }
      _ermg._Alpha[l] = sumT_lj/_n;
  
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

  void Em::mStepMaximizationPi()
  {
    for (int q=0; q<_q; q++){
      _it_Pi_ql = _ermg._Pi.access(q);
      _it_denom_Pi_ql = _denom_Pi.access(q);
      for (int l=0; l<=q; l++){ 
	if (*_it_denom_Pi_ql>PRECISION){
	  *_it_Pi_ql /= *_it_denom_Pi_ql;
	}
	else
	  // class with a single vertex
	  *_it_Pi_ql = 0.5;

	if (*_it_Pi_ql<MSTEPPARAMMIN)
	  *_it_Pi_ql = MSTEPPARAMMIN;  
	if (*_it_Pi_ql>(1-MSTEPPARAMMIN))
	  *_it_Pi_ql = 1-MSTEPPARAMMIN;   
    
	_it_Pi_ql++;
	_it_denom_Pi_ql++;
      }
    }

    // symmetrizing
    for (int q=0; q<_q; q++)
      for (int l=q+1; l<_q; l++) 
	_ermg._Pi.value(q,l) = _ermg._Pi.value(l,q);
  }

}
