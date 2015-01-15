#include <math.h>
#include <cmath>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include "matrix.h"

extern "C" {
void VEstep_( 
	      int &N, int &Q, 
	      double *X_, 
	      double *D_, double *E_, double *CstMat_, 
	      double &mincut, double &maxcut,
	      double &dmaxInner,
	      int    &fpnbiter,
	      int    &directed,
	      double *tau_,
	      int    &nit
	     );
}
 
void VEstep_( 
	     int &N, int &Q, 
	     double *X_, 
	     double *D_, double *E_, double *CstMat_, 
	     double &mincut, double &maxcut,
	     double &dmaxInner,
	     int    &fpnbiter,	      
	     int    &directed,
	     double *tau_,
	     int    &nit
	    )
{

   GetRNGstate(); // Initialize random number generator
  
  // Construct matrices
  Matrix X  ( N, N, X_     , true ); // Do not deallocate X_, D_, ..., CstMat_ 
  Matrix D  ( Q, Q, D_     , true ); 
  Matrix E  ( Q, Q, E_     , true ); 
  Matrix Cst( N, Q, CstMat_, true ); 
  Matrix Tau( N, Q, tau_,    true );

  // Local variables
  Matrix TauD(N,Q), TauE(N,Q), Tau_tE(N,Q);
  double dGamma    = (double) (N*Q);
  double dGammaOld    = 0.9;
  double eps       = DBL_MIN;
  double MaxDouble = DBL_MAX; 
  Matrix TauOld(N, Q, 0.0);
  int i;
  bool old_decrease, decrease = true, go_on=true;

  for( i=0;  
       (dGamma < MaxDouble) && (dGamma > dmaxInner) && (go_on) && (i<fpnbiter); 
       i++
     ) {
    
    dGammaOld = dGamma;
    TauOld    = Tau;

    TauD = Tau * D;
    TauE = Tau * E;

    if ( directed ) { 
      Tau_tE = Tau * E.t();
    }
    
    // RestrictOp are matriciel operations which don't involve diagonal terms
    Tau = Cst + RestrictSum(TauD,"col")  + RestrictProd ( X.t(), TauE );
    if ( directed ) { 
      Tau = Tau + RestrictProd ( X, Tau_tE );
    }

    Tau = min( Tau, maxcut);
    Tau = max( Tau, mincut);
    Tau = exp( Tau );
    Tau = Normalize(Tau, "row");
    Tau = max( Tau, eps ); 
    
    dGamma =  sum( abs(Tau - TauOld) );

    // std::cout << "[C, dir=" << directed <<  "] iter = " << i 
    //	    << ", dGamma = " << dGamma << std::endl;
    if ( i > 1 ) {
      old_decrease = decrease;
      decrease  = ( (dGamma - dGammaOld) < 0 );
      if ( old_decrease & !(decrease) ) {
	go_on = false;
	// std::cout << "Warning: Dgamma increase, stop E-step" <<  std::endl;
      }
    } else if ( i == 1 ) {
      decrease = ( (dGamma - dGammaOld) < 0 );
    }       

  }

  // Must be stored in the output value
  if (go_on)
    Tau.copyData   ( tau_ );
  else
    TauOld.copyData( tau_ );
 
  nit=i-1;

  void PutRNGstate(); // send RNG status to R
}
