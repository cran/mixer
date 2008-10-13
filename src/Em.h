/* Em.h
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
 *       \brief EM algorithm
 *
 *          10/01/2007
 *          
*/

#ifndef ERMG_EM_H
#define ERMG_EM_H
/*!
  \class Em libermg/Em.h
  \brief EM algorithm
*/

#include <EmCore.h>


namespace ermg {

  class Em : public EmCore
  {
  protected:
    // classification method
    bool _classif;

    //! epsilon for convergency of the e_step fixed point
    double _e_step_fixed_point_eps;
    //! convergency of the e_step fixed point
    bool _e_step_fixed_point_convergency;
    //! maximum number of iteration of the e_step fixed point
    int _e_step_fixed_point_nitermax;
    //! copy of maximum number of iteration of the e_step fixed point
    int _back_e_step_fixed_point_nitermax;


    //! transposed Tau at iteration k+1
    BandMemMatrix<double>* _aux_t_Tau;

    //! pointer to transposed Tau at iteration k+1
    BandMemMatrix<double>* _t_Tau2;
    //! iterator points to _t_Tau2[l,j]
    BandMemMatrix<double>::BandCursor _it_t_Tau2_lj;

  
    //! e step
    virtual void eStep();
    //! m step with likelihood computation
    virtual void mStep();
    //! checks if in m first step
    bool _m_firststep;
    //! first m step
    virtual void mFirstStep();

    //! initialization
    virtual bool init( );  
    //! allocate aux_Tau 
    void newAuxTau();
    //! complete the algorithm
    void finish(); 

    //! fixed point iteration
    void eStepFixedPoint();
    //! fixed point iteration for vertex i
    virtual void eStepFixedPointVertex(int i, bool stop_at_diagonal=false);
    //! core calculus of fixed point iteration
    double eStepFixedPointStepCore(SparseMatrix< int, std::vector<int> >& Xadj, double val1, double val2, int i, int l, bool l_is_q, bool stop_at_diagonal=false);
    //! core calculus of fixed point iteration for empty  Xadj[i]
    double eStepFixedPointStepDefault(double val2, int i, int l, bool l_is_q, bool stop_at_diagonal=false);
    //! fixed point iteration normalization and convergency
    void eStepFixedPointStepNormalizeTauAndConvergency(int i, double sum_t_Tau2_qi); 
    //! move _t_Tau2 to _t_Tau1
    void eStepFixedPointSwitchTau2toTau1();



    //! m step and likelihood step for vertex i
    void mStepLikelihoodVertex(int i, bool toestim=true);
    //! useful for mStepLikelihoodStepCore/Default
    double _curr_pi;
    double _curr_denom;
    //! run the m step core calculus
    virtual double mStepLikelihoodVertexRun(double t_Tau1_qi, int i, int q, bool toestim);
    //! run the m step core calculus for empty  Xadj[i]
    virtual double mStepLikelihoodVertexRunDefault(double t_Tau1_qi, int i, int q, bool toestim);
    //! core calculus of m step and likelihood step
    double mStepLikelihoodStepCore( SparseMatrix< int, std::vector<int> >& Xadj, double t_Tau1_qi, int i, int q, int l, bool stop_at_diagonal=false );
    //! core calculus of m step and likelihood step for empty  Xadj[i]
    double mStepLikelihoodStepDefault( double t_Tau1_qi, int i, int q, int l, bool stop_at_diagonal=false );
    //! m step maximization
    virtual void mStepMaximization();
    //! m step maximization for Pi
    virtual void mStepMaximizationPi();

    //! computes the complete likelihood
    virtual void mLikelihood();
    //! likelihood step for vertex i
    void mLikelihoodVertex(int i);
    
    //! forbidden default constructor
    Em();

  public:
    Em( Ermg& ermg,
	int em_nitermax=0, double em_eps=PRECISION,
	int e_step_fixed_point_nitermax=0, double e_step_fixed_point_eps=PRECISION, bool classif=false ); 
    
    virtual ~Em();
    
    //! performs the em step
    virtual bool run(int nbclass);


    //! set minimal nb iteration
    virtual void setMinimalNiter(){
      this->EmCore::setMinimalNiter();
      _back_e_step_fixed_point_nitermax = _e_step_fixed_point_nitermax;
      _e_step_fixed_point_nitermax = 1;
    }
    
    //! unset minimal nb iteration
    virtual void unsetMinimalNiter(){
      this->EmCore::unsetMinimalNiter();
      _e_step_fixed_point_nitermax = _back_e_step_fixed_point_nitermax;
    }

    //! returns the number of e_step fixed point iterations
    int eStepFixedPointNitermax() const{
      return _e_step_fixed_point_nitermax;
    }
    //! sets the number of e_step fixed point iterations
    void eStepFixedPointSetNitermax(int e_step_fixed_point_nitermax){
      _e_step_fixed_point_nitermax = e_step_fixed_point_nitermax;
    }
  };

}
#endif
