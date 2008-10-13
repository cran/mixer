/* Emdg.cc
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
#include <Emd.h>
#include <limits>

using namespace std;


namespace ermg {

  void Emd::eStepFixedPointVertex(int i, bool stop_at_diagonal)
  { 
    Ermdg& ermdg_ref = dynamic_cast< Ermdg& >(_ermg);
    double sum_t_Tau2_qi = 0;
    double maxsumTau_qi = -numeric_limits<double>::max();

    _it_val1_ql = _val1.access(0);
    _it_val2_ql = _val2.access(0);

    int qmax=0;

    if (ermdg_ref._Xadj[i].size()>0){
      if (ermdg_ref._tXadj[i].size()>0){
	for (int q=0; q<_q; q++){
	  double sumTau_qi = 0.;	
	  for (int l=0; l<_q; l++){ 
	    sumTau_qi += eStepFixedPointStepCore( ermdg_ref._Xadj, *_it_val1_ql, *_it_val2_ql, 
						  i, l, l==q );
	    sumTau_qi += eStepFixedPointStepCore( ermdg_ref._tXadj, _val1.value(l,q), _val2.value(l,q), 
						  i, l, l==q );
	    _it_val1_ql++;
	    _it_val2_ql++;      
	  }	
	  _sumsTau_qi[q] = sumTau_qi + log(ermdg_ref._Alpha[q]);
	  if (_sumsTau_qi[q] > maxsumTau_qi){
	    maxsumTau_qi = _sumsTau_qi[q];
	    qmax = q;
	  }
	}
      }
      else{
	for (int q=0; q<_q; q++){
	  double sumTau_qi = 0.;	
	  for (int l=0; l<_q; l++){ 
	    sumTau_qi += eStepFixedPointStepCore( ermdg_ref._Xadj, *_it_val1_ql, *_it_val2_ql, 
						  i, l, l==q );
	    sumTau_qi += eStepFixedPointStepDefault( _val2.value(l,q), 
						     i, l, l==q );
	    _it_val1_ql++;
	    _it_val2_ql++;	  
	  }	
	  _sumsTau_qi[q] = sumTau_qi + log(ermdg_ref._Alpha[q]);
	  if (_sumsTau_qi[q] > maxsumTau_qi){
	    maxsumTau_qi = _sumsTau_qi[q];
	    qmax = q;
	  }
	}      
      }
    }
    else{
      if (ermdg_ref._tXadj[i].size()>0){
	for (int q=0; q<_q; q++){
	  double sumTau_qi = 0.;	
	  for (int l=0; l<_q; l++){ 
	    sumTau_qi += eStepFixedPointStepDefault( *_it_val2_ql, 
						     i, l, l==q );
	    sumTau_qi += eStepFixedPointStepCore( ermdg_ref._tXadj, _val1.value(l,q), _val2.value(l,q), 
						  i, l, l==q );
	    _it_val2_ql++;	  
	  }	
	  _sumsTau_qi[q] = sumTau_qi+ log(ermdg_ref._Alpha[q]);
	  if (_sumsTau_qi[q] > maxsumTau_qi){
	    maxsumTau_qi = _sumsTau_qi[q];
	    qmax = q;
	  }
	}      
      }
      else{
	for (int q=0; q<_q; q++){	
	  double sumTau_qi = 0.;     
	  for (int l=0; l<_q; l++){ 	  
	    sumTau_qi += eStepFixedPointStepDefault( *_it_val2_ql, 
						     i, l, l==q );	  
	    sumTau_qi += eStepFixedPointStepDefault( _val2.value(l,q), 
						     i, l, l==q );	  
	    _it_val2_ql++;      
	  }	
	  _sumsTau_qi[q] = sumTau_qi + log(ermdg_ref._Alpha[q]);
	  if (_sumsTau_qi[q] > maxsumTau_qi){
	    maxsumTau_qi = _sumsTau_qi[q];
	    qmax = q;
	  }
	}
      }
    }

    if (!_classif){
      // factorization
      for (int q=0; q<_q; q++){
	// with factorization
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

  void Emd::mStep()
  { 
    this->initPiDenomPi();
    _ermg.resetLikelihood();
    for (int i=0; i<_n; i++){
      this->mStepLikelihoodVertex(i);
      this->mStepLikelihoodComplementaryVertex(i); 
    }
    this->mStepMaximization();
  }

  void Emd::mLikelihood()
  { 
    _ermg.resetLikelihood();
    for (int i=0; i<_n; i++){
      this->mStepLikelihoodVertex(i, false);
      this->mStepLikelihoodComplementaryVertex(i); 
    }
  }  

  double Emd::mStepLikelihoodVertexRun(double t_Tau1_qi, int i, int q, bool toestim)
  {
    double tmplikelihood = 0.;
    for (int l=0; l<_ermg._q; l++){      
      tmplikelihood += mStepLikelihoodStepCore( _ermg._Xadj, t_Tau1_qi, i, q, l, false );
      if (toestim){
	*_it_Pi_ql += _curr_pi;
	*_it_denom_Pi_ql += _curr_denom;
	_it_Pi_ql++;      
	_it_denom_Pi_ql++;
      }
      _it_val1_ql++;
      _it_val2_ql++;
    }
    return tmplikelihood;
  }

  double Emd::mStepLikelihoodVertexRunDefault(double t_Tau1_qi, int i, int q, bool toestim)
  { 
    double tmplikelihood = 0.;
    for (int l=0; l<_ermg._q; l++){
      tmplikelihood += mStepLikelihoodStepDefault( t_Tau1_qi, i, q, l, false );
      if (toestim){
	*_it_denom_Pi_ql += _curr_denom;
	_it_denom_Pi_ql++;
      }
      _it_val2_ql++;
    }
    return tmplikelihood;
  }

  void Emd::mStepLikelihoodComplementaryVertex(int i)
  {
    Ermdg& ermdg_ref = dynamic_cast< Ermdg& >(_ermg);

    if (ermdg_ref._tXadj[i].size()>0){
      double tmplikelihood = 0.;
      for (int q=0; q<_q; q++){ 
	double t_Tau1_qi = (*_t_Tau1).value(q,i);	
	for (int l=0; l<_q; l++){
	  tmplikelihood += mStepLikelihoodComplementaryStepCore( ermdg_ref._tXadj, t_Tau1_qi, i, q, l );
	} 
      }  
      _ermg.addToCompleteLikelihood( tmplikelihood );
    }
    else{
      double tmplikelihood = 0.;
      for (int q=0; q<_q; q++){
	double t_Tau1_qi = (*_t_Tau1).value(q,i);      
	for (int l=0; l<_q; l++){
	  tmplikelihood += mStepLikelihoodComplementaryStepDefault( t_Tau1_qi, i, q, l );
	}
      }
      _ermg.addToCompleteLikelihood( tmplikelihood );
    }
  }

  double Emd::mStepLikelihoodComplementaryStepCore( SparseMatrix< int, std::vector<int> >& tXadj, 
						    double t_Tau1_qi, int i, int q, int l )
  {
    double c;
    int s=tXadj[i].size(), ind1, ind2;
    double tmplikelihood = 0.;

    // initialize the beta_ijql as 2 possible values 
    double val1 = _val1.value(l,q);
    double val2 = _val2.value(l,q);

    int pos =tXadj.firstPositionAfterDiagonal(i);
    // avoiding tXadj[i][i]
    if (tXadj.fullDiagonal(i)){
      pos--;
    }


    //// LOWER than i
    _it_Xadj_i = tXadj[i].begin();	     
    ind1=-1;  
    
    _it_t_Tau1_lj = (*_t_Tau1).access(l, ind1+1);
    // all the non null elements
    for (int nbe=0; nbe<pos; nbe++){
      ind2 = *_it_Xadj_i;
      
      // all the null elements between
      for (int j=ind1+1; j<ind2; j++){
	c = t_Tau1_qi * (*_it_t_Tau1_lj);
	tmplikelihood += c*val2;
	_it_t_Tau1_lj++;
      }
      // non null element
      c = t_Tau1_qi*(*_it_t_Tau1_lj);
      tmplikelihood += c*val1;
      ind1 = ind2;
      _it_t_Tau1_lj++;	
      _it_Xadj_i++;
    }
    // last null elements
    for (int j=ind1+1; j<i; j++){
      c = t_Tau1_qi*(*_it_t_Tau1_lj);
      tmplikelihood += c*val2;
      _it_t_Tau1_lj++;
    }


    //// EQUAL to i (no process)
    if (tXadj.fullDiagonal(i)){
      // zapping _Xadj[i][i]
      pos++;
      _it_Xadj_i++;
    }	
    _it_t_Tau1_lj++;

 
    //// GREATER than i  (Tau1 concerned)
    ind1=i; 
    //_it_t_Tau1_lj = (*_t_Tau1).access(l, ind1+1);

    // all the non null elements
    for (int nbe=pos; nbe<s; nbe++){
      ind2 = *_it_Xadj_i;
      
      // all the null elements between
      for (int j=ind1+1; j<ind2; j++){
	c = t_Tau1_qi * (*_it_t_Tau1_lj);
	tmplikelihood += c*val2;
	_it_t_Tau1_lj++;
      }
      // non null element
      c = t_Tau1_qi*(*_it_t_Tau1_lj);
      tmplikelihood += c*val1;
      ind1 = ind2;
      _it_t_Tau1_lj++;	
      _it_Xadj_i++;
    }
    // last null elements
    for (int j=ind1+1; j<_n; j++){
      c = t_Tau1_qi*(*_it_t_Tau1_lj);
      tmplikelihood += c*val2;
      _it_t_Tau1_lj++;
    }

    return tmplikelihood;
  }

  double Emd::mStepLikelihoodComplementaryStepDefault( double t_Tau1_qi, int i, int q, int l )
  {
    double c;
    int ind1;
    double tmplikelihood = 0.;

    // initialize the only possible beta_ijql value
    double val2 = _val2.value(l,q);
    
    //// LOWER than i 
    ind1=-1;
    _it_t_Tau1_lj = (*_t_Tau1).access(l, ind1+1);
    // all null elements
    for (int j=ind1+1; j<i; j++){
      c = t_Tau1_qi * (*_it_t_Tau1_lj);
      tmplikelihood += c*val2;
      _it_t_Tau1_lj++;
    }    
    
    //// EQUAL to i (no process)
    _it_t_Tau1_lj++;
        
    //// LOWER than i 
    // all null elements
    ind1=i;  
    for (int j=ind1+1; j<_n; j++){
      c = t_Tau1_qi * (*_it_t_Tau1_lj);
      tmplikelihood += c*val2;
      _it_t_Tau1_lj++;
    }
    
    return tmplikelihood;
  }

  void Emd::mStepMaximizationPi()
  {
    _it_Pi_ql = _ermg._Pi.access(0);
    _it_denom_Pi_ql = _denom_Pi.access(0);
    for (int q=0; q<_ermg._q; q++){
      for (int l=0; l<_ermg._q; l++){ 
	
	if (*_it_denom_Pi_ql>PRECISION)
	  *_it_Pi_ql /= *_it_denom_Pi_ql;
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
  }
}
