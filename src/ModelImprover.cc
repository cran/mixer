
/* ModelImprover.cc
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

#include <ModelImprover.h>
#include <limits>
#include <unistd.h>
#include <time.h>
using namespace std;


namespace ermg {

  ModelImprover::ModelImprover( Ermg& ermg, 
				int qmin, int qmax,
				int nitermax, bool silent )
    : _curr_ermg(ermg)
  {
    _qmin = qmin;
    _qmax = qmax; 
    _nitermax = nitermax;
    _silent = silent;
    _ed_rep = vector<ErmgDescriptor>(_qmax-_qmin+1);
    int cm = _qmax-_qmin;
    _imove_usedtosplit.assign(cm+1, -1);
    _imove_usedtomerge.assign(cm+1, -1);
    _imove.assign(cm+1, 0);
  }


  void ModelImprover::initialLikelihoods(EmCore& em, bool nokmeans)
  {
    if (em.requiredNumberOfInitializedVertices()==_curr_ermg.nbVertices()){
      if (nokmeans){
	_curr_ermg.cahInitNoKmeans();
      }
      else{
	if (!_silent)
	  cout<<"K-means...\t"<<flush; 
	_curr_ermg.kmeans();
	if (!_silent)
	  cout<<"done"<<endl;
      } 
      if (!_silent)
	cout<<"Hierarchical Clustering...\t\t"<<flush; 
      _curr_ermg.cah();
      if (!_silent)
	cout<<"done"<<endl;
    }
    else{
      if (!_silent)
	cout<<"Hierarchical Clustering ...\t\t"<<flush; 
      _curr_ermg.cahForSubXadj( em.requiredNumberOfInitializedVertices() );      
      if (!_silent)
	cout<<"done"<<endl;
    }
  
    for (int c=_qmin; c<=_qmax; c++){ 
      if (!_silent)	
	cout<<"\rEM for Q="<<c<< flush; 
      _curr_ermg.em( em, c );
      _curr_ermg.saveTo( _ed_rep[c-_qmin] );      
    }
    if (!_silent)
      cout<<"\rEM for Q="
	<<_qmin
	<<" to "
	<<_qmax
	<<"...\tdone"<<endl;
  }


  void ModelImprover::loadLikelihoods(const string& spmfile, const string& idir)
  {
    string::size_type islash = spmfile.find_last_of("/");
    for (int c=_qmin; c<=_qmax; c++){ 
      stringstream strin;
      strin<<idir<<"/"<<spmfile.substr(islash+1, spmfile.size()-4-islash-1)<<"_Q"<<c<<".model";
      if (!_silent)	
	cout<<"Loading "<<strin.str()<<"...\t"<< flush; 
      _curr_ermg.inFile( strin.str() );
      _curr_ermg.saveTo( _ed_rep[c-_qmin] );      
      if (!_silent)
	cout<<"done"<<endl;
    }
  }


//   const vector<int> ModelImprover::notIncreasedLikelihood()
//   {  
//     vector<int> tmpnum_ed;
//     for (int c=1; c<=_qmax-_qmin; c++){
//       if (_ed_rep[c].incompleteLikelihoodApproximation()<_ed_rep[c-1].incompleteLikelihoodApproximation()){
// 	tmpnum_ed.push_back( c );
// #ifdef VERBOSE
// 	cerr<<_qmin+c<<" not increase"<<endl;
// #endif
//       }
//     }
// #ifdef VERBOSE
//     cerr<<tmpnum_ed.size()<<" notIncreased Likelihoods"<<endl;
// #endif
//     return tmpnum_ed;
//   }
  //
//   const vector<int> ModelImprover::notConvexLikelihood()
//   { 
//     vector<int> tmpnum_ed;
//     for (int c=1; c<_qmax-_qmin; c++){
//       if ( (_ed_rep[c-1].incompleteLikelihoodApproximation()+_ed_rep[c+1].incompleteLikelihoodApproximation()
// 	    - 2*_ed_rep[c].incompleteLikelihoodApproximation())>0 ){
// 	tmpnum_ed.push_back( c );
// #ifdef VERBOSE
// 	cerr<<_qmin+c<<" not convex"<<endl;
// #endif
//       }
//     }
// #ifdef VERBOSE
//     cerr<<tmpnum_ed.size()<<" notConvex Likelihoods"<<endl;
// #endif
//     return tmpnum_ed;
//   }

  const vector<int> ModelImprover::notConvexLikelihood(ModelImproverLog& implog)
  { 
    vector<int> tmpnum_ed;
    int c=0;
    int cm = _qmax-_qmin;
    while(c<cm){
      double max = -std::numeric_limits<double>::max();
      int kmax = c+1;
      for (int k=c+1; k<=cm; k++){
	double pk = (_ed_rep[k].incompleteLikelihoodApproximation()-_ed_rep[c].incompleteLikelihoodApproximation()) / (k-c);
	// 	cerr<<"between "<<c+1<<" and "<<k+1<<": "<<
	// 	  _ed_rep[k].incompleteLikelihoodApproximation()<<"-"<<_ed_rep[c].incompleteLikelihoodApproximation()<<"/"
	// 	    <<(k-c)<<"="<<pk<<endl;
	if (pk>max){
	  kmax = k;
	  max = pk;
	}
      }
      implog.inenvelopeLog(c);
      implog.inenvelopeLog(kmax);
      //      cerr<<c+1<<" and "<<kmax+1<<" in envelope"<<endl;
      for (int k=c+1; k<kmax; k++){
#ifdef VERBOSE
	cerr<<_qmin+k<<" not convex"<<endl;
#endif
	tmpnum_ed.push_back(k);
      }
      c=kmax;
    }
    tmpnum_ed.push_back(cm);
    return tmpnum_ed;
  }


  void ModelImprover::improveLikelihoodsSpeed(EmCore& em)
  { 
    if (!_silent)
      cout<<"Improvement...\t"<< flush;
    int cm = _qmax-_qmin;
    vector<int> candidate(1,0);
    int niter = 0;
    bool statuquo;

    do{ 
      //clock_t start, end;
      //start = clock();
      ModelImproverLog implog(_qmin, _qmax);
      if (!_silent)
	cout<<"\rImprovement...\t"<<niter+1<<flush;
      statuquo = true;	
      bool mstatuquo, sstatuquo; 
      int c=0;
      while(c<cm){
	//0// envelope calculus
	double max = -std::numeric_limits<double>::max();
	int kmax = c+1;
	for (int k=c+1; k<=cm; k++){
	  double pk = (_ed_rep[k].incompleteLikelihoodApproximation()-_ed_rep[c].incompleteLikelihoodApproximation()) / (k-c);
	  if (pk>max){
	    kmax = k;
	    max = pk;
	  }
	}

	implog.inenvelopeLog(c);
	implog.inenvelopeLog(kmax);
	if (kmax-c>1){
	  //1// envelope extrema
	  candidate[0] = c+1;	  
	  sstatuquo = splitStrategy( em,  candidate, implog );
	  candidate[0] = kmax-1;
	  mstatuquo = mergeStrategy( em,  candidate, implog );
	  if ( (!mstatuquo) || (!sstatuquo) )
	    statuquo = false;
	}
	if (kmax-c>2){
	  int kselect;
	  //2.1// most decreasing point
	  double diffmax = 0.;
	  for (int k=c+1; k<kmax; k++){
	    double diff = _ed_rep[k-1].incompleteLikelihoodApproximation()-_ed_rep[k].incompleteLikelihoodApproximation();
	    if (diff>diffmax){
	      diffmax = diff;
	      kselect = k;
#ifdef VERBOSE
	      cerr<<_qmin+k<<" not increase"<<endl;
#endif	    
	    }
	  }
	  if (false==true){
	  //if ( (diffmax>0) && (_imove_usedtosplit[kselect-1] != _imove[kselect-1]) ){
	    candidate[0] = kselect;
	    sstatuquo = splitStrategy( em,  candidate, implog );
	    implog.selectstatusLog(kselect);	    
	  }
	  else{
	    //2.2// or most away from the envelope
	    double y1 = _ed_rep[c].incompleteLikelihoodApproximation();
	    double y2 = _ed_rep[kmax].incompleteLikelihoodApproximation();
	    int x1 = c, x2=kmax;
	    double a = -(y2-y1)/(x2-x1);
	    double b = 1.;
	    double cste = -(y1-a*x1);
	    // distance to envelope
	    double dmax = 0;
	    for (int k=c+1; k<kmax; k++){
	      #ifdef VERBOSE
	      cerr<<_qmin+k<<" not convex"<<endl;
	      #endif
	      double y3 = _ed_rep[k].incompleteLikelihoodApproximation();
	      int x3 = k;
	      double d = abs( (a*x3+b*y3+cste)/sqrt(a*a+b*b) );
	      //cerr<<"distance: "<<d<<endl;
	      if (d>dmax){
		dmax = d;
		kselect = k;
	      }
	    }
	    #ifdef VERBOSE
	    cerr<<_qmin+kselect<<" is the most away from the envelope"<<endl;
	    #endif
	    if ( (_imove_usedtomerge[kselect+1] == _imove[kselect+1])
		 &&(_imove_usedtosplit[kselect-1] == _imove[kselect-1]) ){
	      //2.3// or random drawing
	      double tmprand = double(rand())/double(RAND_MAX)*(kmax-c-2);
	      int krand=1;
	      while (krand<tmprand){
		krand++;
	      }
	      #ifdef VERBOSE
	      cerr<<"Random replacement of "<<kselect+1;
	      #endif
	      if (c+krand>=kselect)
		kselect = c+krand+1;
	      else
		kselect = c+krand;
	      implog.randomstatusLog(kselect);
	      #ifdef VERBOSE
	      cerr<<" by "<<kselect+1<<endl;
	      #endif
	    }
	    else{
	      #ifdef VERBOSE
	      cerr<<"No random replacement of "<<kselect+1<<endl;
	      #endif
	      implog.selectstatusLog(kselect);
	    }

	    if (kselect-1!=c){
	      candidate[0] = kselect;
	      sstatuquo = splitStrategy( em,  candidate, implog );
	    }
	    if (kselect+1!=kmax){
	      candidate[0] =  kselect;	  
	      mstatuquo = mergeStrategy( em,  candidate, implog );
	    }
	    if ( (!mstatuquo) || (!sstatuquo) )
	      statuquo = false;
	  }
	}
	
	c=kmax;
      }
      // 3 extrema
      candidate[0] = cm;
      sstatuquo = splitStrategy( em,  candidate, implog );

      niter++;
      implog.likeLog(_ed_rep);
      //cerr<<implog;
      //end = clock();
      //double elapsed = ((double)end - start) / CLOCKS_PER_SEC; /* Conversion en seconde */
      //cerr<<"#Time: "<<elapsed<<endl<<endl<<endl;
    } while( (!statuquo)&&(niter<_nitermax) );    
  
    if (!_silent)
      cout<<"\rImprovement...\tdone                          "<<endl;   
  }
   
  void ModelImprover::improveLikelihoods(EmCore& em)
  {  
    if (!_silent)
      cout<<"Improvement...\t"<< flush; 
    int niter;
    bool statuquo;
    vector<int> candidate(1,0);
    niter=0;
    do{
      //clock_t start, end;
      //start = clock();
      ModelImproverLog implog(_qmin, _qmax);
      if (!_silent)
 	cout<<"\rImprovement...\t"<<niter+1<<flush;
      statuquo = true; 
      vector<int> num_ed =  this->notConvexLikelihood(implog);
      for (vector<int>::iterator it=num_ed.begin(); it!=num_ed.end(); it++){
 	candidate[0] = *it;
 	bool mstatuquo = mergeStrategy( em,  candidate, implog );
 	bool sstatuquo = splitStrategy( em, candidate, implog );
 	if ( (!mstatuquo) || (!sstatuquo) )
 	  statuquo = false;
      }
      niter++;
      implog.likeLog(_ed_rep);
      //cerr<<implog;
      //end = clock();
      //double elapsed = ((double)end - start) / CLOCKS_PER_SEC; /* Conversion en seconde */
      //cerr<<"#Time: "<<elapsed<<endl<<endl<<endl;;
    } while( (!statuquo)&&(niter<_nitermax) );    
    
    if (!_silent)
      cout<<"\rImprovement...\tdone                          "<<endl;   
  }    

  bool ModelImprover::mergeStrategy(EmCore& em, const vector<int>& num_ed, ModelImproverLog& implog)
  {
    bool statuquo = true;

    vector<int>::const_reverse_iterator tmpit = num_ed.rbegin();
    for (;tmpit!=num_ed.rend();tmpit++){

      if (*tmpit+1<=(_qmax-_qmin)){
#ifdef VERBOSE
	cerr<<"merging Q="<<_qmin+ *tmpit+1<<" ";
#endif
	if (_imove_usedtomerge[*tmpit+1] != _imove[*tmpit+1]){
	  _imove_usedtomerge[*tmpit+1] = _imove[*tmpit+1];
	  _curr_ermg.loadFrom( _ed_rep[*tmpit+1] );
	  _curr_ermg.merge();
	  _curr_ermg.em( em, _curr_ermg.nbClasses() );
	  if ( (_curr_ermg.incompleteLikelihoodApproximation()-_ed_rep[*tmpit].incompleteLikelihoodApproximation())>PRECISION ){	    
#ifdef VERBOSE
	    cerr<<"accepted"<<endl;
#endif
	    statuquo = false;
	    _curr_ermg.saveTo( _ed_rep[*tmpit] );	  
	    _imove[*tmpit]++;
	  }
#ifdef VERBOSE
	  else
	    cerr<<"rejected"<<endl;
#endif 
	  implog.mergeLog(statuquo, *tmpit);
	}
#ifdef VERBOSE
	else 
	  cerr<<" previously done"<<endl;
#endif 
      }  
    }
    return statuquo;
  }

  bool ModelImprover::splitStrategy(EmCore& em, const vector<int>& num_ed, ModelImproverLog& implog)
  {
    ErmgDescriptor ed;

    int _snitermax;
    // if dynamic_cast<SpecErmg*>(_ermg)
    _snitermax = 1;
    //else
    //_snitermax = SNITERMAX;

    bool statuquo = true;
    vector<int>::const_iterator tmpit = num_ed.begin();
    for (;tmpit!=num_ed.end();tmpit++){
      if (*tmpit-1>=0){
#ifdef VERBOSE
	cerr<<"spliting Q="<<_qmin+ *tmpit-1<<" ";
#endif
	if (_imove_usedtosplit[*tmpit-1] != _imove[*tmpit-1]){
	  _imove_usedtosplit[*tmpit-1] = _imove[*tmpit-1];
	  double ilikea =_ed_rep[*tmpit].incompleteLikelihoodApproximation();

	  //! selecting the best split candidate
	  bool sstatuquo = true;
	  int sniter = 0;
	  double max = -INF;	
	  em.setMinimalNiter() ; 
	  while( (sstatuquo)&&(sniter<_snitermax) ){

	    _curr_ermg.loadFrom( _ed_rep[*tmpit-1] );
	    _curr_ermg.split();
	    _curr_ermg.em( em, _curr_ermg.nbClasses() );
	    double lik = _curr_ermg.incompleteLikelihoodApproximation();
	    if ( lik>max ){
	      _curr_ermg.saveTo( ed );
	      max = lik;
	    }
	    if ( (lik-ilikea)>PRECISION )
	      sstatuquo = false;
	
	    sniter++;
	  }

	  //! performing EM on the candidate
	  em.unsetMinimalNiter() ;
	  _curr_ermg.loadFrom( ed );
	  _curr_ermg.em( em, _curr_ermg.nbClasses() );	  
	  if ( (_curr_ermg.incompleteLikelihoodApproximation()-ilikea)>PRECISION ){
#ifdef VERBOSE
	    cerr<<"accepted"<<endl;
#endif
	    statuquo = false;
	    _curr_ermg.saveTo( _ed_rep[*tmpit] );
	    _imove[*tmpit]++;
	  }
#ifdef VERBOSE
	  else
	    cerr<<"rejected"<<endl;
#endif 
	  implog.splitLog(statuquo, *tmpit);	  
	}
#ifdef VERBOSE
	else 
	  cerr<<" previously done"<<endl;
#endif 
      }
    }
    return statuquo;
  }


  void ModelImprover::printLikelihoods()
  {
    if (!_silent)
      for (int c=0; c<=_qmax-_qmin; c++)
	cout<<"# incomplete Likelihood Approximation"<<endl
	    <<_qmin+c<<"\t"<<_ed_rep[c].incompleteLikelihoodApproximation()<<endl;
    cout<<endl<<endl;
  }

  void ModelImprover::outFile(const string& ofile, const Graph *g, const string& odir )
  {
    if (!_silent)
      cout<<"Saving...\t"<< flush; 
    for (int c=0; c<=_qmax-_qmin; c++){
      _curr_ermg.loadFrom( _ed_rep[c] );
      _curr_ermg.outFile( ofile, g, odir );
    }

    stringstream strout;
    if (odir.length()>0){
      string::size_type islash = ofile.find_last_of("/");
      strout<<odir<<"/"<<ofile.substr(islash+1, ofile.size()-4-islash-1)<<".likelihoods";
    }
    else{
      strout<<ofile.substr(0, ofile.size()-4)<<".likelihoods";
    }
    ofstream fout;
    fout.open(strout.str().c_str(), ios::out);
    for (int c=0; c<=_qmax-_qmin; c++){
      fout<<"# Q="<<_qmin+c
	  <<" Incomplete Likelihood Approximation:\n"<<_ed_rep[c].incompleteLikelihoodApproximation()<<endl;
    }
    fout<<endl<<endl;
  
    for (int c=0; c<=_qmax-_qmin; c++){
      fout<<"# Q="<<_qmin+c
	  <<" Entropy:\n"<<_ed_rep[c]._entropy<<endl;
    } 
    fout<<endl<<endl;

    for (int c=0; c<=_qmax-_qmin; c++){
      fout<<"# Q="<<_qmin+c
	  <<" ICL:\n"<<_ed_rep[c]._icl<<endl;
    }  
    fout.close(); 
    if (!_silent)
      cout<<"done"<<endl; 
  }

}
