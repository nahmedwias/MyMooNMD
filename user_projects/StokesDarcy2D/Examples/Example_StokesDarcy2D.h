#ifndef _EXAMPLE_STOKESDARCY2D_
#define _EXAMPLE_STOKESDARCY2D_

#include<Example_NSE2D.h>
#include<Example_Darcy2D.h>
#include<Example_CD2D.h>

#include <memory>

/**
 * currently supported examples are
 * 
 * 0: sin_cos_polynomial_BJS
 * 1: polynomial_ut0
 * 2: riverbed
 * 3: ana_sol_Stokes_Darcy.3
 */

class Example_StokesDarcy2D : public Example2D
{

  // coefficients for Stokes in 'problem_coefficients' in Example2D part of 
  // this class
  CoeffFct2D *mixed_darcy_coeffs;
  CoeffFct2D *primal_darcy_coeffs;
  
  public:
  
  Example_StokesDarcy2D();

  ~Example_StokesDarcy2D(){};

  /** getters */
  std::shared_ptr<Example_NSE2D>   get_stokes_example() const;
   // for mixed darcy
  std::shared_ptr<Example_Darcy2D> get_mixed_darcy_example() const;
  // for primal darcy
  std::shared_ptr<Example_CD2D>    get_primal_darcy_example() const;
};


#endif