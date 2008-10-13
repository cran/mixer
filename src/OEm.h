/* OEm.h
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
 *       \brief Online EM algorithm
 *
 *          10/01/2007
 *          
*/

#ifndef ERMG_OEM_H
#define ERMG_OEM_H
/*!
  \class OEm libermg/OEm.h
  \brief Online EM algorithm
*/

#include <Em.h>


namespace ermg {

  class OEm : public Em
  {
  protected:
    //! current index
    int _curr_i;

    //! dot_Tau[q]
    std::vector<double> _dotTau;
    //! theta[q,l]
    BandMemMatrix<double> _Theta; 
    //! iterator points to _Theta[l,j]
    BandMemMatrix<double>::BandCursor _it_Theta_ql; 
    //! backup of t_Tau1 for vertex i
    std::vector<double> _t_Tau1_back;

    //! before em    
    virtual bool init();

    //! e step
    virtual void eStep();

    //! first m step
    virtual void mFirstStep();
    //! m step
    virtual void mStep();     
    //! prepare m step for vertex i by substraction for old class of i 
    void mStepPrepareVertex(int i, int old_classi);   
    //! m step for vertex i
    virtual void mStepVertex(int i, bool stop_at_diagonal); 
    //! m step and likelihood step
    virtual void mStepMaximization();
    //! m step maximization for Pi
    virtual void mStepMaximizationPi();



  public:   
    OEm( Ermg& ermg, int n_init, int em_nitermax=3 )
      :  Em(ermg){
      _curr_i=0;
      _n_init=n_init;
      _em_nitermax=em_nitermax;
    }    
      
      bool run(int nbclass);      
      
      // returns the number of vertices to be initialized
      virtual int requiredNumberOfInitializedVertices() const{
	return _n_init;
      }
  };

}
#endif
