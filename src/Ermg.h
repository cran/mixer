/* Ermg.h
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
 *       \brief Erdös Reni Mixture for Graphs
 *
 *          13/03/2006
 *          
*/

#ifndef ERMG_ERMG_H
#define ERMG_ERMG_H
/*!
  \class Ermg libermg/Ermg.h
  \brief Erdös Reni Mixture for Graphs class
*/


#include <cmath>
#include <Barycenters.h>
#include <SparseMatrix.h>
#include <BandMemMatrix.h>
#include <Graph.h>


namespace ermg {

  const double INF=1000000000.0;
  const double TAUMIN=1e-10;
  const double MSTEPPARAMMIN=1e-6;

  class EmCore;

  class ErmgDescriptor{    
  public:
    //! classes number
    int _q;
    //! complete likelihood
    double _complete_likelihood;
    //! entropy
    double _entropy;
    //! ICL
    double _icl;
    //! alpha
    std::vector<double> _Alpha;
    //! Pi
    BandMemMatrix<double> _Pi;
    //! transposed Tau
    BandMemMatrix<double> _t_Tau;
    //! return the imcomplete Likelihood Approximation
    double incompleteLikelihoodApproximation(){
      return (_complete_likelihood - _entropy);
    }
    //! returns the ICL
    double ICL(){ 
      return _icl;
    }
  };


  class Ermg
  {   
  protected:
    //! complete likelihood
    double _complete_likelihood;
    //! entropy
    double _entropy;

    //! graph as a sparse matrix
    SparseMatrix< int, std::vector<int> > _Xadj;
    //! vertices number
    int _n;
    //! classes number
    int _q;

    //! allows edge from a vertex to itself
    bool _enable_loop;

    //! initialization before em
    bool emInit( int nbclass );   
    //! initialization of _class and _alpha from _save_class before em
    bool emInitClassAlpha();
    //! initialization of Tau from the CAH _class
    void emInitTau();
    //! allocate Tau
    void newTau();


    //! transposed Tau at iteration k
    BandMemMatrix<double>* _t_Tau;
    //! iterator points to _t_Tau[l,j]
    BandMemMatrix<double>::BandCursor _it_t_Tau_lj;


    //! alpha
    std::vector<double> _Alpha;
 
    //! Pi
    BandMemMatrix<double> _Pi;



    //! Barycenters
    Barycenters _G;
    //! initialize the Barycenters object
    virtual void initBarycenters(int nbclass, int n){
      _G.init(nbclass, n);
    }
     //! computes the barycenters from class
    virtual void initBarycentersFromClass(){
      _G.initFromClass(_class, _cardinal_class);
    }
    //! initiallize barycenters from selected vertices
    virtual void initBarycentersFromSelected(const std::vector<int>& selected){
      _G.initFromSelected(selected);
    }

    //! class of each vertex
    std::vector<int> _class;
    //! save of the class at the end of CAH_step
    std::vector< std::vector<int> > _save_class;
    //! iterator points to _class
    std::vector<int>::iterator _it_class;
    //! cardinal of each class
    std::vector<int> _cardinal_class;
    //! iterator points to _cardinal_class
    std::vector<int>::iterator _it_card_c1;
    //! iterator points to _cardinal_class
    std::vector<int>::iterator _it_card_c2;
    //! returns the number of unempty classes
    int currentNbClass();
  
  
    //! kmeans: number of class
    int _kmeans_nbclass;
    //! kmeans: concerned classes
    std::vector<int> _kmeans_concerned_classes;
    //! kmeans: concerned vertices
    std::vector<int> _kmeans_concerned_vertices;
    //! CAH: current number of class
    int _curr_nbclass;
    //! CAH: min number of class to save
    int _cah_minnbclass;
    //! CAH: max number of class to save
    int _cah_maxnbclass;
    //! maximum number of iteration of kmeans
    int _kmeans_nitermax;


    //! kmeans core code
    void kmeansCore();
    //! kmeans step 
    double kmeansStep(int i, bool first=false);
    //! computes the kmeans distance vertex i and barycenter c
    virtual  double kmeansDistance(int i, int c){ 
      return _G.kmeansDistance(i, c);
    }
    //! update the barycenter min_c
    virtual void kmeansUpdateBarycenters(int i, int min_c){  
      _G.kmeansUpdate(i, min_c);
    }
    //! normalizes and check convergency
    virtual bool kmeansNormalizeBarycentersAndConvergency( bool check_convergency ){  
      return _G.normalizeAndConvergency(_cardinal_class, check_convergency);
    }
    //! switch last barycenters to previous
    virtual void kmeansSwitchBarycenters(){
      _G.switchG2toG1();
    }
    //! initiallize barycenters
    void kmeansInitializeBarycenters();
    //! initiallize barycenters from random vertices of _Xadj
    void kmeansInitializeBarycentersCore();
    //! initiallize barycenters as all the vertices of _Xadj
    void kmeansInitializeBarycentersDefault();
    //! select randomly the _nbclass vertices for the first random initialization
    void kmeansRandomSelect(std::vector<int>& selected);


    //! CAH step
    void cahStep();
    //! updates class by replacing c2 by c1
    void cahUpdateClass(int c1, int c2);
    //! operates when a single represented class 
    void cahSingleClass();
    //! cah distance between classes c1 and c2
    virtual double cahDistance(int c1,  int c2){
      return _G.cahDistance(c1, c2);
    }    
    //! update the barycenter min_c1 from itself and min_c2
    virtual void cahUpdateBarycenters(int min_c1, int cardinal_min_c1, int min_c2, int cardinal_min_c2){
      _G.cahUpdate(min_c1, cardinal_min_c1, min_c2, cardinal_min_c2);
    }
    //! prepares the cah for a sub Xadj of size n_init
    void cahInitNoKmeansFromSubXadj(int n_init);


    //! construct _class from Tau
    void constructClassFromTau(bool tosave=true);
    //! construct _cardinal from Class; returns the number of cardinal>0
    int constructCardinalFromClass();
    //! construct _G from Class
    bool constructBarycentersFromClass();
 

    //! performs the split step with kmeans
    void splitKmeans();
    //! performs the split step with spectral clustering
    void splitSpectral();

    // reset the likelihood
    void resetLikelihood();
    // add l to complete likelihood
    void addToCompleteLikelihood(double l){
      _complete_likelihood += l;
    }

    //! associated Graph if required
    const Graph* _gr;

    //! print with the Graph labels
    void printWithLabels(std::ostream& out, const Graph *gr){
      _gr = gr;
      out<<*this;
      _gr = NULL;
    }

    //! error file message
    void errFormatFile(){
      std::cerr<<"Error format file"<<std::endl;
      exit(1);
    }

    //! forbidden default constructor
    Ermg();



  public:
    //! constructor
    Ermg( const std::vector< std::vector<int> >& spm, bool enable_loop,
	  int kmeans_nitermax,
	  int kmeans_nbclass, 
	  int cah_minnbclass, int cah_maxnbclass );
 
    //! destructor
    virtual ~Ermg();

    //! copy constructor
    Ermg(const Ermg& ermg);
    //! operator=
    Ermg& operator=(const Ermg& ermg);

    //! returns the number of classes
    int nbClasses() const{
      return _q;
    } 
    //! returns the size of _Xadj, i.e. the number of vertices
    int nbVertices() const{
      return _n;
    }

    //! performs the kmeans step
    void kmeans();
    //! prepares the cah when no kmeans before
    void cahInitNoKmeans();
    //! performs the cah
    void cah();
    //! performs the cah on a sub Xadj of size n_init
    void cahForSubXadj(int n_init);
    //! replace the cah by loading classes
    void replaceCahForSubXadjByLoad(const std::vector<int>& inputclass );
    //! saves class
    std::vector<int>& savedClass(int q=0);
    
    
    //! performs the em step 
    bool em(EmCore& em, int nbclass);
    
    //! performs the split step
    void split();
    //! performs the merge step
    void merge();
    //! performs a simple model relaod step
    void noSplitMerge();
    //! returns the entropy for the class q 
    double entropy(int q);


    //! returns the complete likelihood
    double completeLikelihood(){
      return _complete_likelihood;
    }
    //! returns the entropy
    double entropy(){
      return _entropy;
    }
    //! returns the incomplete likelihood approximation
    double incompleteLikelihoodApproximation(){
      return (_complete_likelihood - _entropy);
    }
    //! returns the ICL
    virtual double ICL(){ 
      return this->completeLikelihood() - (_q-1)/2.0*log(double(_n)) - _q*(_q+1)/4.0*log(double(_n*(_n-1)/2.0));
    }
    //! returns the BIC
    virtual double BIC(){ 
      return this->incompleteLikelihoodApproximation() - (_q-1)/2.0*log(double(_n)) - _q*(_q+1)/4.0*log(double(_n*(_n-1)/2.0));
    }

    //! saves the model in the file ofile in directory odir
    void outFile(const std::string& ofile, const std::string& odir = std::string());  
    //! saves the model in the file ofile in directory odir with Graph labels
    void outFile(const std::string& ofile, const Graph *g, const std::string& odir = std::string());
    //! loads the model from the file ifile
    void inFile(const std::string& ifile);

    //! save the model to a ErmgDescriptor
    void saveTo(ErmgDescriptor& ed_out);
    //! load the model from a ErmgDescriptor
    void loadFrom(const ErmgDescriptor& ed_in);

  
    friend 
      std::ostream& operator<<(std::ostream& out, Ermg& ermg);
    friend 
       std::istream&  operator>>( std::istream& in, Ermg& ermg);

  
    friend class EmCore;
    friend class Em;
    friend class Emd;
    friend class OCEmT;
    friend class OCEm;
    friend class SOCEm;
    friend class OEm;
  };
}
#endif

