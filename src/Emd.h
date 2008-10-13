/* Emd.h
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
 *       \brief EM algorithm for directed graphs
 *
 *          10/01/2007
 *          
*/

#ifndef ERMG_EMD_H
#define ERMG_EMD_H
/*!
  \class Ermg libermg/Em.h
  \brief EM algorithm
*/


#include <Ermdg.h>
#include <Em.h>

namespace ermg {

  class Emd : public Em
  {
  protected:
    //! iterator points to _tXadj[i]
    std::vector<int>::iterator _it_tXadj_i;
    
    //! fixed point iteration for vertex i
    virtual void eStepFixedPointVertex(int i, bool stop_at_diagonal=false);
    
    // m Step
    virtual void mStep();
    //! run the m step core calculus
    virtual double mStepLikelihoodVertexRun(double t_Tau1_qi, int i, int q, bool toestim);
    //! run the m step core calculus for empty  Xadj[i]
    virtual double mStepLikelihoodVertexRunDefault(double t_Tau1_qi, int i, int q, bool toestim);
    //! m step and likelihood complementary step using _tXadj only for the likelihood for vertex i
    void mStepLikelihoodComplementaryVertex(int i);  
    //! core calculus of m step and likelihood step
    double mStepLikelihoodComplementaryStepCore( SparseMatrix< int, std::vector<int> >& Xadj, double t_Tau1_qi, int i, int q, int l );
    //! core calculus of m step and likelihood step for empty  Xadj[i]
    double mStepLikelihoodComplementaryStepDefault( double t_Tau1_qi, int i, int q, int l );
    //! m step maximization for Pi
    virtual void mStepMaximizationPi();

    //! computes the complete likelihood
    virtual void mLikelihood();

  public:    
    Emd( Ermdg& ermdg,
	 int em_nitermax, double em_eps,
	 int e_step_fixed_point_nitermax, double e_step_fixed_point_eps, bool classif=false )
      :  Em(ermdg, em_nitermax, em_eps, e_step_fixed_point_nitermax, e_step_fixed_point_eps, classif)
      {}
      
      virtual ~Emd(){};
  };

}
#endif
