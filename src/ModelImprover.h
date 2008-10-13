/* ModelImprover.h
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
 *       \brief Erdös Reni Mixture for Graphs selection
 *
 *          13/10/2006
 *          
*/

#ifndef ERMG_MODELIMPROVER_H
#define ERMG_MODELIMPROVER_H
/*!
  \class ModelImprover.h libermg/ModelImprover.h
  \brief Erdös Reni Mixture for Graphs Selection  class
*/


#include <cmath>
#include <Emd.h>
#include <Ermdg.h>

namespace ermg {

  const int SNITERMAX=5;
  class ModelImproverLog;
  class ModelImprover
  {   
  protected: 
    //! minimal number of classes
    int _qmin;
    //! maximal number of classes
    int _qmax;
    //! maximal number of iteration 
    int _nitermax;
    //! silent mode
    bool _silent;

    //! Ermg Desriptor repository
    std::vector<ErmgDescriptor> _ed_rep;
    //! curent Ermg
    Ermg& _curr_ermg;

    //! index in terms of number of moves 
    //! for each number of classes
    std::vector<int> _imove;
    //! index in terms of number of moves used for the previous split
    std::vector<int> _imove_usedtosplit;
    //! index in terms of number of moves used for the previous merge
    std::vector<int> _imove_usedtomerge;

    //! checks the models with a decreasing likelihood
    //const std::vector<int> notIncreasedLikelihood();
    //! checks the models with a unconvex likelihood
    const std::vector<int> notConvexLikelihood(ModelImproverLog& implog);
    //! merge
    bool mergeStrategy(EmCore& em, const std::vector<int>& num_ed, ModelImproverLog& implog);
    //! split
    bool splitStrategy(EmCore& em, const std::vector<int>& num_ed, ModelImproverLog& implog);

  public:
    //! constructor
    ModelImprover( Ermg& ermg,
		   int qmin, int qmax, 
		   int nitermax, bool silent=false );
    //! destructor
    ~ModelImprover(){
    }
    //! performs the initial EMs with/out kmeans
    void initialLikelihoods(EmCore& em, bool nokmeans=true);

    //! loads the initial models
    void loadLikelihoods(const std::string& spmfile, const std::string& idir);

    void improveLikelihoodsSpeed(EmCore& em);
    void improveLikelihoods(EmCore& em);

    void printLikelihoods();

    void outFile(const std::string& ofile, const Graph *g, const std::string& odir = std::string());

    ErmgDescriptor& getModelForNbclass(int q){
      return _ed_rep[q-_qmin];
    }
  };

  class ModelImproverLog
  {
  private:
    int _qmin; 
    int _qmax;    
    std::vector<int> _splitlog;
    std::vector<int> _mergelog;
    std::deque<bool> _inenvelopelog;
    std::vector<char> _statuslog;
    std::vector<double> _likelog;
    
    ModelImproverLog(){
    }
  public:
    ModelImproverLog(int qmin, int qmax){
      _qmin = qmin;
      _qmax = qmax;
      init();
    }
    ~ModelImproverLog(){}

    void init(){
      _splitlog.assign(_qmax-_qmin+1, 0);
      _mergelog.assign(_qmax-_qmin+1, 0);
      _inenvelopelog.assign(_qmax-_qmin+1, 0);
      _statuslog.assign(_qmax-_qmin+1, '-');
      _likelog.assign(_qmax-_qmin+1, 0.);
    }
    void splitLog(bool sstatuquo, int k){
      if (sstatuquo)
	_splitlog[k] = -1;
      else
	_splitlog[k] = 1;
    }
    void mergeLog(bool mstatuquo, int k){
      if (mstatuquo)
	_mergelog[k] = -1;
      else
	_mergelog[k] = 1;
    }
    void inenvelopeLog(int k){
      _inenvelopelog[k] = true;
    }     
    void likeLog(std::vector<ErmgDescriptor>& ed_rep){
      for (int c=0; c<=_qmax-_qmin; c++){
	_likelog[c] = ed_rep[c].incompleteLikelihoodApproximation();
      }
    }
    void randomstatusLog(int k){
      _statuslog[k] = 'R';
    }
    void selectstatusLog(int k){
      _statuslog[k] = 'S';
    }    
    friend std::ostream&  operator<<(std::ostream& out, ModelImproverLog& mil){
      out<<"#Q\tsplit\tmerge\tenv\tstatus\tlogl"<<std::endl;
      for (int c=0; c<=mil._qmax-mil._qmin; c++){	
	out<<mil._qmin+c<<"\t"<<mil._splitlog[c]<<"\t"<<mil._mergelog[c]<<"\t"
	   <<mil._inenvelopelog[c]<<"\t"<<mil._statuslog[c]<<"\t"
	   <<mil._likelog[c]<<std::endl;
      }
      return out; 
    }	  
  };

}
#endif
