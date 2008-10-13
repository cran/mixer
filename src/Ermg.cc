/* Ermg.cc
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

#include <Ermg.h>
#include <EmCore.h>

#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iterator>
#include <limits>
#include <set>

using namespace std;


namespace ermg {

  Ermg::Ermg( const vector< vector<int> >& spm,  
	      bool enable_loop, 
	      int kmeans_nitermax,
	      int kmeans_nbclass, 
	      int cah_minnbclass, int cah_maxnbclass )
    : _Xadj( spm ), _G(_Xadj)
  {
#ifdef VERBOSE
    _Xadj.stderrPrint();
#endif
    
    _n = _Xadj.nLines();
    _enable_loop=enable_loop;
    _kmeans_nitermax = kmeans_nitermax;
    _kmeans_nbclass = kmeans_nbclass;
    _cah_minnbclass = cah_minnbclass;
    _cah_maxnbclass = cah_maxnbclass;
    _q = -1;

    _t_Tau = NULL;
    
    _gr = NULL;
  }


  Ermg::~Ermg()
  {
    if (_t_Tau){
      delete _t_Tau;
      _t_Tau=NULL;
    }
  }

  bool Ermg::em( EmCore& em, int nbclass )
  {
    if (!this->emInit( nbclass ))
      return false; 
    return em.run( nbclass );
  }

  bool Ermg::emInit( int nbclass )
  {
    _q = nbclass;

    _Alpha.assign(_q, 0.);
    _Pi.assign(_q, _q, 0.);

    this->newTau();
    this->emInitClassAlpha();

    return true;
  }
  
  bool Ermg::emInitClassAlpha()
  {
    if (_q>1)
      _class = this->savedClass( _q );
    else
      _class.assign(_n,0);
    _cardinal_class.assign(_q, 0);
    
    std::vector<int>::iterator it_class = _class.begin();
    
    int curr_num=0;
    vector<int> corresp(_kmeans_nbclass,-1);
    
    int defined_class=0;
    for (int i=0; i<_n; i++){
      if (*it_class>=0){
	defined_class++;
	if (corresp[*it_class]==-1){
	  corresp[*it_class]=curr_num;
	  _cardinal_class[*it_class]=0;
	  curr_num++;
	}
	*it_class = corresp[*it_class];
	_cardinal_class[*it_class]++;
      }
      it_class++;
    }

    int q=0;
    for(q=0; q<_q; q++)
      _Alpha[q] = _cardinal_class[q]*1.0/ defined_class;
    q=0;
    while( (_cardinal_class[q] > 1) && (q<int(_cardinal_class.size()-1))){
	q++;
    }
    if( (_cardinal_class[q] > 1) && (q==int(_cardinal_class.size()-1)) )
	q++;

    if (q!=int(_cardinal_class.size()))
      return false;
    return true;
  }

  void Ermg::newTau()
  {
    if (_t_Tau==NULL)
      _t_Tau = new BandMemMatrix<double>(_q, _n, 0.);
    else
      if (((_t_Tau)->nLines()!=_q)||((_t_Tau)->nCols()!=_n)){
	delete _t_Tau;
	_t_Tau = new BandMemMatrix<double>(_q, _n, 0.);
      }
  }

  void Ermg::emInitTau()
  {
    std::vector<int>::iterator it_class;
    _it_t_Tau_lj = (*_t_Tau).access(0);
    for (int l=0; l<_q; l++){
      it_class = _class.begin();
      for (int j=0; j<_n; j++){
	if ((*it_class)==l){
	  *_it_t_Tau_lj = 1;
	  _it_t_Tau_lj++;
	  it_class++;
	}
	//if (_class[j]==l){
	//(*_t_Tau).value(l,j) = 1;
	//}
	else{
	  *_it_t_Tau_lj = 0;
	  _it_t_Tau_lj++;
	  it_class++;
	  //(*_t_Tau).value(l,j) = 0;
	}
      }
    }
  }



  void Ermg::resetLikelihood()
  {
    _complete_likelihood = 0.;
    _entropy = 0.;
  }


  double Ermg::entropy(int q)
  {
    double entropy = 0.;
    _it_t_Tau_lj = (*_t_Tau).access(q);
    for (int i=0; i<_n; i++){
      //entropy += (*_it_t_Tau_lj)*log(*_it_t_Tau_lj);
      //_it_t_Tau_lj++;
      entropy += ((*_t_Tau).value(q,i))*log((*_t_Tau).value(q,i));    
    }
    return -entropy;
  }


  void Ermg::merge()
  {
    if (_q<3){
      cerr<<"Merge impossible with only "<<_q<<" classes"<<endl;
      exit(1);
    }   
    _kmeans_nbclass = _q; 
  
    this->constructClassFromTau(true);
    this->constructCardinalFromClass();
    if ( this->constructBarycentersFromClass() ){    
      _cah_maxnbclass = _q; 
      _cah_minnbclass = _q-1;
      _save_class.clear();
      this->cah();
    }
    else{
      cerr<<"Merge warning: a single class assigned for all the vertices"<<endl;
      _save_class.clear();
      _class.assign(_n,0);
      _save_class.push_back(_class);
    }
    _q--;
  }


  void Ermg::split()
  {
    int cselect = 0;
    double emax = -numeric_limits<double>::max();
    for(int c=0; c<_q; c++){    
      double e = this->entropy(c);
      if (e>emax){
	emax=e;
	cselect=c;
      }
    }

    constructClassFromTau(true);
#ifdef VERBOSE
    constructCardinalFromClass();
    cerr<<cselect<<" of cardinal "<<_cardinal_class[cselect]<<" ";
#endif
    // filling concerned
    _kmeans_concerned_classes.clear();
    _kmeans_concerned_classes.push_back(cselect);
    _kmeans_concerned_classes.push_back(_q);
    _kmeans_concerned_vertices.clear();
    for (int i=0; i<_n; i++)
      if (_class[i]==cselect)
	_kmeans_concerned_vertices.push_back(i);
    _q++;

    if (_kmeans_concerned_vertices.size()>4)
      splitSpectral();
    //splitKmeans();
  
    _save_class.clear();
    _save_class.push_back(_class);
  }

  void Ermg::splitSpectral()
  {
    _kmeans_nbclass++;
    // Focus on the second vp
    SparseMatrix< int, std::vector<int> > subXadj(_Xadj, _kmeans_concerned_vertices);
    //_Xadj.stderrPrint();
    //cerr<<endl;
    //subXadj.stderrPrint();
    vector<double>& vp = subXadj.secondEigenVector();
    if (!vp.empty()){
      for (int j=0; j<subXadj.nLines(); j++){
	if (vp[j]>0)
	  _class[_kmeans_concerned_vertices[j]] = _kmeans_concerned_classes[0];
	else
	  _class[_kmeans_concerned_vertices[j]] = _q-1;
      }      
    }
  }

  void Ermg::splitKmeans()
  {
    int save_kmeans_nbclass = _kmeans_nbclass+1;
    if (_kmeans_concerned_vertices.size()>0){
      _kmeans_nbclass = _kmeans_concerned_classes.size();
      this->initBarycenters(_kmeans_nbclass, _n); 
      this->kmeansCore();
    }
    _kmeans_nbclass = save_kmeans_nbclass;
  }


  void Ermg::noSplitMerge()
  {
    this->constructClassFromTau();
  }

 
  void Ermg::kmeans()
  {  
    _class.assign(_n,0);

    _kmeans_concerned_classes.clear();
    for(int c=0; c<_kmeans_nbclass; c++)
      _kmeans_concerned_classes.push_back(c);
    _kmeans_concerned_vertices.clear();
    _kmeans_concerned_vertices.reserve(_n);
    for(int i=0; i<_n; i++)
      _kmeans_concerned_vertices.push_back(i);

    this->initBarycenters(_kmeans_nbclass, _n); 
    this->kmeansCore();
  }


  void Ermg::kmeansCore()
  {  
    int niter = 0;
    double entropy;
  
    this->kmeansInitializeBarycenters();

    // first iterate
    _cardinal_class.assign(_kmeans_nbclass, 0);
    entropy = 0;
    _it_class = _class.begin();
    for (int i=0; i<int(_kmeans_concerned_vertices.size()); i++){
      entropy += this->kmeansStep(_kmeans_concerned_vertices[i], true);
      _it_class++;
    }
    this->kmeansNormalizeBarycentersAndConvergency(false);
#ifdef VERBOSE
    cerr<<"Entropy :"<<entropy<<endl;
    cerr<<"=====> KMeans cardinal class"<<endl;
    copy(_cardinal_class.begin(),_cardinal_class.end(),ostream_iterator<int>(cerr, " "));
    cerr<<endl;
#endif

    // next  
    do{
      this->kmeansSwitchBarycenters();
      _cardinal_class.assign(_kmeans_nbclass, 0);
      entropy = 0; 
      _it_class = _class.begin();
      for (int i=0; i<int(_kmeans_concerned_vertices.size()); i++){
	entropy += this->kmeansStep(_kmeans_concerned_vertices[i]);
	_it_class++;
      }
      niter++;
#ifdef VERBOSE
      cerr<<"Entropy :"<<entropy<<endl;
      cerr<<"=====> KMeans cardinal class"<<endl;
      copy(_cardinal_class.begin(),_cardinal_class.end(),ostream_iterator<int>(cerr, " "));
      cerr<<endl; 
#endif 
    } 
    while((this->kmeansNormalizeBarycentersAndConvergency(true)==false) && (niter<_kmeans_nitermax));
#ifdef VERBOSE
    cerr<<"Niter kmeans: "<<niter<<endl;
#endif 
    this->kmeansSwitchBarycenters();


#ifdef VERBOSE
    cerr<<"KMeans class"<<endl;
    for (int i=0; i<int(_kmeans_concerned_vertices.size()); i++)
      cerr<<i<<" => "<<_class[i]<<endl;; 
    cerr<<"Cardinal"<<endl;
    for (int c=0; c<_kmeans_nbclass; c++)
      cerr<<"class "<<_kmeans_concerned_classes[c]<<" : "<<_cardinal_class[c]<<endl;  
#endif  
  }


  double Ermg::kmeansStep(int i, bool first)
  {
    double min_score_i = INF, entropy=0.;
    int min_c = -1;

 
    // calculating the distance between element _Xadj[i] 
    // and barycentre G1[c,.]
    vector<int> exaequo_min;
    for (int c=0; c<_kmeans_nbclass; c++){

      double score_ic = kmeansDistance(i, c);	  
    
      if (score_ic<min_score_i){
	exaequo_min.clear();
	min_score_i = score_ic;
	min_c = c;
      }
      else if ( (score_ic<(min_score_i+PRECISION)) && (score_ic>(min_score_i-PRECISION)) )
	exaequo_min.push_back(c);    
    }
  
    // in case of ex-aequo
    if (exaequo_min.size()>0){
      exaequo_min.push_back(min_c);    
      int nbexaequo = exaequo_min.size();
      double tmprand = double(rand())/double(RAND_MAX)*nbexaequo;
      int c=1;
      while (c<tmprand){
	c++;
      }
      min_c = exaequo_min[c-1];
    }
    entropy = min_score_i;   
  
 
    // Update _class
    *_it_class = _kmeans_concerned_classes[min_c];
    //_class[i] = _kmeans_concerned_classes[min_c];
    _cardinal_class[ min_c ]++;
  
    // updating _G
    this->kmeansUpdateBarycenters(i, min_c);

    return entropy;
  }

  void Ermg::kmeansInitializeBarycenters()
  { 
      if (_kmeans_nbclass<_n)
	kmeansInitializeBarycentersCore();
      else
	kmeansInitializeBarycentersDefault();	
  }

  void Ermg::kmeansInitializeBarycentersCore()
  {    
    vector<int> selected(_kmeans_nbclass);
    this->kmeansRandomSelect( selected );
    this->initBarycentersFromSelected(selected);  
  }

  void Ermg::kmeansInitializeBarycentersDefault()
  {
    vector<int> selected(_kmeans_nbclass);
    for (int i=0; i<_n; i++)
      selected[i] = i; // cas CAHNOKMEANS
    this->initBarycentersFromSelected(selected);    
  }


  void Ermg::kmeansRandomSelect(vector<int>& selected)
  {
    double tmpcoeff;
    double step = double(_kmeans_concerned_vertices.size())/double(_kmeans_nbclass);
    for (int c=0; c<_kmeans_nbclass; c++){
      tmpcoeff = double(rand())/double(RAND_MAX);
      selected[c] = _kmeans_concerned_vertices[ int( (c+tmpcoeff)*step ) ];
    }
#ifdef VERBOSE      
    cerr<<"Vertices selected: ";
    for (int c=0; c<_kmeans_nbclass; c++)
      cerr<<selected[c]<<" ";
    cerr<<endl;
#endif
  }


  void Ermg::cahInitNoKmeans()
  {
    _kmeans_nbclass = _n;
    this->initBarycenters(_kmeans_nbclass, _n);
    this->kmeansInitializeBarycenters();
  
    _cardinal_class.assign(_kmeans_nbclass, 0); 
    for (int i=0; i<_kmeans_nbclass; i++){
      _cardinal_class[i] = 1;
    }
    _class.assign(_n,0);
    for (int i=0; i<_n; i++){
      _class[i] = i;
    }
  } 
 
  void Ermg::cahInitNoKmeansFromSubXadj(int n_init)
  { 
    int back_n = _n;
    
    _n=n_init;
    _kmeans_nbclass = _n;
    
    this->initBarycenters(_kmeans_nbclass, _n);
    this->kmeansInitializeBarycenters();
    
    _cardinal_class.assign(_kmeans_nbclass, 0); 
    for (int i=0; i<_kmeans_nbclass; i++){
      _cardinal_class[i] = 1;
    }
    _class.assign(_n,0);
    for (int i=0; i<_n; i++){
      _class[i] = i;
    }
    
    _n = back_n;
  }
  

  int Ermg::currentNbClass()
  {
    if (_kmeans_nbclass == _n)
      return _n;
    else{
      int cnb = 0;
      for (_it_card_c1=_cardinal_class.begin(); _it_card_c1!=_cardinal_class.end();  _it_card_c1++)
	if (*_it_card_c1>0)
	  cnb ++;
      return cnb;
    }
  }

  void Ermg::cah()
  {
    _curr_nbclass = this->currentNbClass();

    while (_curr_nbclass>_cah_maxnbclass+1){
      this->cahStep();
    }
  
    if (_curr_nbclass <= _cah_maxnbclass){
      if (_curr_nbclass <= _cah_minnbclass)
	_save_class.push_back(_class);
      else{
	for (int i=_cah_maxnbclass; i>=_curr_nbclass; i--)
	  _save_class.push_back(_class);	
      }
    }

    while (_curr_nbclass>_cah_minnbclass){
      this->cahStep();
      _save_class.push_back(_class);
    }
  }

  void Ermg::cahForSubXadj(int n_init)
  {
    this->cahInitNoKmeansFromSubXadj(n_init);
    
    _curr_nbclass = this->currentNbClass();

    while (_curr_nbclass>_cah_maxnbclass+1){
      this->cahStep();
    }
  
    vector<int> aux(_n, -1);
    if (_curr_nbclass <= _cah_maxnbclass){
      if (_curr_nbclass <= _cah_minnbclass){
	copy(_class.begin(), _class.end(), aux.begin());
	_save_class.push_back(aux);
      }
      else{
	for (int i=_cah_maxnbclass; i>=_curr_nbclass; i--){
	  copy(_class.begin(), _class.end(), aux.begin());
	  _save_class.push_back(aux);	
	}
      }
    }

    while (_curr_nbclass>_cah_minnbclass){
      this->cahStep();
      copy(_class.begin(), _class.end(), aux.begin());
      _save_class.push_back(aux);   
    }
  }

  void Ermg::replaceCahForSubXadjByLoad(const vector<int>& inputclass )
  {
    _curr_nbclass = _cah_minnbclass;
    _kmeans_nbclass = _cah_minnbclass;

    vector<int> aux(_n, -1);
    int m = *max_element(inputclass.begin(), inputclass.end());
    if ( m>(_cah_minnbclass-1) ){
      cerr<<"Error: in input class file, the classes goes from 0 to "<<m<<"  >= "<<_cah_minnbclass<<" classes"<<endl;
      exit(1);
    }
    copy(inputclass.begin(), inputclass.end(), aux.begin());
    _save_class.push_back(aux);
  }

   
  void Ermg::cahStep()
  {
    if (_curr_nbclass==2)
      this->cahSingleClass();
    else{
      double min_score = INF;
      int min_c1=0, min_c2=0;
    
      int nbexaequo=0;
      vector<int> exaequo_min;

      _it_card_c1 = _cardinal_class.begin();
      for (int c1=0; c1<_kmeans_nbclass; c1++){
	_it_card_c2 = _cardinal_class.begin();
	//if(_cardinal_class[c1]!=0){
	if(*_it_card_c1!=0){

	  for (int c2=0; c2<c1; c2++){
	    //if(_cardinal_class[c2]!=0){
	    if(*_it_card_c2!=0){
	    
	      //double score_c1c2 = (_cardinal_class[c1]*_cardinal_class[c2])/double(_cardinal_class[c1]+_cardinal_class[c2])
	      double score_c1c2 = ((*_it_card_c1)*(*_it_card_c2))/double((*_it_card_c1)+(*_it_card_c2))
		* this->cahDistance( c1, c2 );

	      if (score_c1c2<min_score){
		min_c1 = c1;
		min_c2 = c2;
		min_score = score_c1c2; 
	      
		nbexaequo=0;
		exaequo_min.clear();
		exaequo_min.push_back(c1);
		exaequo_min.push_back(c2);
	      
	      }	  
	      else if ( (score_c1c2<(min_score+PRECISION))
			&& (score_c1c2>(min_score-PRECISION)) ) {
		nbexaequo++;
		exaequo_min.push_back(c1);
		exaequo_min.push_back(c2);
	      }	  	
	    }

	    _it_card_c2++;
	  }
	}

	_it_card_c1++;
      }

      // if a single class 
      if ((min_c1==0) && (min_c2==0)){
	cerr<<"Hierarchical Clustering warning: a single class assigned for all the vertices"<<endl;
	this->cahSingleClass();
      }
      else{
	// in case of ex-aequo
	if (exaequo_min.size()>1){
	  double tmprand = double(rand())/double(RAND_MAX)*(exaequo_min.size()/2);
	  int c=1;
	  while (c<tmprand){
	    c++;
	  }
	  min_c1 = exaequo_min[(c-1)*2];
	  min_c2 = exaequo_min[(c-1)*2+1];
	}
      
	int c;
	if (min_c1>min_c2){
	  c=min_c2;
	  min_c2=min_c1;
	  min_c1=c;
	}
	int card = _cardinal_class[min_c1] + _cardinal_class[min_c2];
	int card1 = _cardinal_class[min_c1];
	int card2 = _cardinal_class[min_c2];
#ifdef VERBOSE	
	cerr<<min_c1<<" and "<<min_c2<<" of cardinal "<<card1<<" and "<<card2<<" ";
#endif  
	this->cahUpdateBarycenters(min_c1, card1, min_c2, card2);
      
	_cardinal_class[min_c1] = card;
	_cardinal_class[min_c2] = 0;
      
	this->cahUpdateClass(min_c1, min_c2);
      }
    }
  }

  void Ermg::cahSingleClass()
  {
    fill(_cardinal_class.begin(), _cardinal_class.end(), 0);
    _cardinal_class[0] = _n;
    _curr_nbclass--;
  }


  void Ermg::cahUpdateClass(int c1, int c2)
  {
    replace( _class.begin(), _class.end(), c2, c1 );
    _curr_nbclass--; 
  }

  

  vector<int>& Ermg::savedClass(int q){    
    if (_save_class.size()>1)
      return _save_class[_cah_maxnbclass-q];
    else
      return _save_class[0];
  }


  void Ermg::constructClassFromTau(bool tosave)
  {
    _class.assign(_n,0);
    for (int i=0; i<_n; i++){
      double max=0;
      int imax=-1;
      for (int q=0; q<_q; q++){
	if( (*(_t_Tau)).value(q,i)>max ){
	  imax = q;
	  max = (*(_t_Tau)).value(q,i);
	}
      }
      _class[i] = imax;
    } 
    if (tosave)
      _save_class.push_back( _class );

#ifdef VERBOSE 
    cerr<<"Constructed Class from Tau:"<<endl;
    for (int i=0; i<_n; i++)
      cerr<<_class[i]<<" ";
    cerr<<endl; 
#endif  
  }


  int Ermg::constructCardinalFromClass()
  {
    _cardinal_class.assign(_q, 0);
    for (int i=0; i<_n; i++)
      _cardinal_class[ _class[i] ]++;
    int nbnnc = 0;
    for (int q=0; q<_q; q++)
      if (_cardinal_class[q]>0)
	nbnnc++;
    return nbnnc;
  }

  bool Ermg::constructBarycentersFromClass()
  { 
    if (this->constructCardinalFromClass()>1){
      this->initBarycenters(_kmeans_nbclass, _n);
      this->initBarycentersFromClass();
      return true;
    }
    else
      return false;
  }


  ostream&  operator<<(ostream& out, Ermg& ermg)
  {
    out<<"# Q\n"<<ermg._q<<endl;

    out<<"# Alpha\n";
    for (int q=0; q<ermg._q; q++)
      out<<ermg._Alpha[q]<<" ";
    out<<endl;
 
    out<<"# Pi\n";
    for (int q=0; q<ermg._q; q++){
      for (int l=0; l<ermg._q; l++)
	out<<ermg._Pi.value(q,l)<<" ";
      out<<endl;
    }
    out.precision(10);
    out<<"# Incomplete Likelihood Approximation\n"<<ermg.completeLikelihood()<<endl;
    out<<"# Entropy\n"<<ermg.entropy()<<endl;
    out<<"# ICL\n"<<ermg.ICL()<<endl;

    out.precision(5);
    out<<"# Tau\n";
    if (ermg._gr){
      for (int i=0; i<ermg._n; i++){
	out<<ermg._gr->getLabel(i)<<" ";
	for (int q=0; q<ermg._q; q++)
	  out<<(*(ermg._t_Tau)).value(q,i)<<" ";
	out<<endl;
      }
    }
    else{
      for (int i=0; i<ermg._n; i++){
	out<<i<<" ";
	for (int q=0; q<ermg._q; q++)
	  out<<(*(ermg._t_Tau)).value(q,i)<<" ";
	out<<endl;
      }
    }
    

#ifdef VERBOSE
    vector<int> cardinal_class(ermg._q, 0);
    for (int i=0; i<ermg._n; i++){
      //    out<<i<<" => ";
      double max=0;
      int imax=-1;
      for (int q=0; q<ermg._q; q++){
	if( (*(ermg._t_Tau)).value(q,i)>max ){
	  imax = q;
	  max = (*(ermg._t_Tau)).value(q,i);
	}
      }
      cardinal_class[imax]++;
    } 

    cerr<<"Cardinal"<<endl;
    for (int c=0; c<ermg._q; c++)
      cerr<<"class "<<c<<" : "<<cardinal_class[c]<<endl; 
    cerr<<endl;
#endif  
    return out; 
  }

  istream&  operator>>(istream& in, Ermg& ermg)
  {
    string line;
    string currstr;
    istringstream linestream;

    // # Q
    getline(in,line);
    linestream.clear(); linestream.str(line);
    linestream >> currstr;
    linestream >> currstr;
    if (currstr!="Q"){
      ermg.errFormatFile();
    }
    // q
    getline(in,line);
    linestream.clear(); linestream.str(line);
    linestream >> currstr;
    if (currstr.length()==0){
      ermg.errFormatFile();
    }
    ermg._q = atoi( currstr.c_str() );


    if (int(ermg._Alpha.size())!=ermg._q)
      ermg._Alpha.assign(ermg._q, 0);
    if (((ermg._Pi).nLines()!=ermg._q)||((ermg._Pi).nCols()!=ermg._q))
      ermg._Pi.assign(ermg._q, ermg._q, 0.);
  
    ermg.newTau();
  
    // # Alpha
    getline(in,line);
    linestream.clear(); linestream.str(line);
    linestream >> currstr;
    linestream >> currstr;
    if (currstr!="Alpha"){
      ermg.errFormatFile();
    }
    // Alpha    
    getline(in,line);
    linestream.clear(); linestream.str(line);
    for (int q=0; q<ermg._q; q++){
      linestream >> currstr;
      if (currstr.length()==0){
	ermg.errFormatFile();
      }
      ermg._Alpha[q] = atof( currstr.c_str() );
    }

    // # Pi
    getline(in,line);
    linestream.clear(); linestream.str(line);
    linestream >> currstr;
    linestream >> currstr;
    if (currstr!="Pi"){
      ermg.errFormatFile();
    }    
    // Pi
    for (int q=0; q<ermg._q; q++){  
      getline(in,line);
      linestream.clear(); linestream.str(line);    
      for (int l=0; l<ermg._q; l++){
	linestream >> currstr;
	if (currstr.length()==0){
	  ermg.errFormatFile();
	}
	ermg._Pi.value(q,l) = atof( currstr.c_str() );
      }
    }
    // # Likelihood
    getline(in,line);
    linestream.clear(); linestream.str(line);
    linestream >> currstr;
    linestream >> currstr;
    linestream >> currstr;
    if (currstr!="Likelihood"){
      ermg.errFormatFile();
    }
    // Complete Likelihood
    getline(in,line);
    linestream.clear(); linestream.str(line);
    linestream >> currstr;
    if (currstr.length()==0){
      ermg.errFormatFile();
    }
    double incomplikeapp = atof( currstr.c_str() );

    // # Entropy
    getline(in,line);
    linestream.clear(); linestream.str(line);
    linestream >> currstr;
    linestream >> currstr;
    linestream >> currstr;
    if (currstr!="Entropy"){
      ermg.errFormatFile();
    }
    // Entropy
    getline(in,line);
    linestream.clear(); linestream.str(line);
    linestream >> currstr;
    if (currstr.length()==0){
      ermg.errFormatFile();
    }
    ermg._entropy = atof( currstr.c_str() );
    ermg._complete_likelihood = incomplikeapp + ermg._entropy;

    // # ICL
    getline(in,line);
    getline(in,line);

    // # Tau
    getline(in,line);
    linestream.clear(); linestream.str(line);
    linestream >> currstr;
    linestream >> currstr;
    if (currstr!="Tau"){
      ermg.errFormatFile();
    } 
    // Tau   
    for (int i=0; i<ermg._n; i++){  
      getline(in,line);
      linestream.clear(); linestream.str(line);
      linestream >> currstr; // label
      for (int q=0; q<ermg._q; q++){
	linestream >> currstr;
	if (currstr.length()==0){
	  ermg.errFormatFile();
	}
	(*(ermg._t_Tau)).value(q,i) = atof( currstr.c_str() );
      }
    }

    ermg._cah_maxnbclass = ermg._q;
    ermg._cah_minnbclass = ermg._q;
    ermg._kmeans_nbclass = ermg._q;

    return in;
  }

void Ermg::outFile(const std::string& ofile, const std::string& odir)
{      
      std::stringstream strout;
      if (odir.length()>0){
	std::string::size_type islash = ofile.find_last_of("/");
	strout<<odir<<"/"<<ofile.substr(islash+1, ofile.size()-4-islash-1)<<"_Q"<<_q<<".model";
      }
      else{
	strout<<ofile.substr(0, ofile.size()-4)<<"_Q"<<_q<<".model";
      }
      std::ofstream fout;
      fout.open(strout.str().c_str(), std::ios::out);
      fout<<*this;
      fout.close();
}

void Ermg::outFile(const std::string& ofile, const Graph *g, const std::string& odir)
{
      std::stringstream strout;
      if (odir.length()>0){
	std::string::size_type islash = ofile.find_last_of("/");
	strout<<odir<<"/"<<ofile.substr(islash+1, ofile.size()-4-islash-1)<<"_Q"<<_q<<".model";
      }
      else{
	strout<<ofile.substr(0, ofile.size()-4)<<"_Q"<<_q<<".model";
      }
      std::ofstream fout;
      fout.open(strout.str().c_str(), std::ios::out);
      this->printWithLabels(fout, g);
      fout.close();
}

void Ermg::inFile(const std::string& ifile)
{
  std::ifstream fin;
  fin.open(ifile.c_str(), std::ios::in);
  if (!fin.is_open()){
    std::cerr<<"Unable to open "<<ifile.c_str()<<std::endl;
    exit(1);
  }
  fin>>*this;
  fin.close();
}



void Ermg::saveTo(ErmgDescriptor& ed_out)
{
  ed_out._q = _q;
  ed_out._complete_likelihood = _complete_likelihood; 
  ed_out._entropy = _entropy;
  ed_out._icl = ICL();
  ed_out._Alpha = _Alpha;
  ed_out._Pi = _Pi;
  ed_out._t_Tau = *(_t_Tau);
}


void Ermg::loadFrom(const ErmgDescriptor& ed_in)
{
  _q = ed_in._q;
  _complete_likelihood = ed_in._complete_likelihood; 
  _entropy = ed_in._entropy; 
  _Alpha = ed_in._Alpha;
  _Pi = ed_in._Pi; 

  newTau();
  *(_t_Tau) = ed_in._t_Tau;
  
  _cah_maxnbclass = _q;
  _cah_minnbclass = _q;
  _kmeans_nbclass = _q;
}

}
