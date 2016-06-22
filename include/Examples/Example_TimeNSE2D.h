/** ************************************************************************ 
*
* @class Example_TimeNSE2D
* @brief store all functions needed to describe a (Navier-)Stokes example 
* 
* The standard
* constructor of this class will fill the vectors (in Example2D_TimeProblem2D) with pointers
* to the functions needed to fully describe a particular example.
* 
* @author 
* @date 
* 
* @ruleof0
* 
 ************************************************************************  */
#ifndef _Example_TimeNSE2D_
#define _Example_TimeNSE2D_

#include <Example_NonStationary2D.h>

class Example_TimeNSE2D : public Example_NonStationary2D
{
public:
  /** @brief default constructor
   * This intializes a Time dependent (Navier-)Stokes example in 2D. 
   * It is chosen according to example_code.
   */
  Example_TimeNSE2D(int example_code,
                    const ParameterDatabase& user_input_parameter_db);
  /** @brief initialize your own example
   * 
   * Create an example with all vectors already defined.
   */
  Example_TimeNSE2D(std::vector <DoubleFunct2D*> exact,
                   std::vector <BoundCondFunct2D*> bc,
                   std::vector <BoundValueFunct2D*> bd, 
                   CoeffFct2D *coeffs,
                   bool timedependentrhs, bool timedependentcoeffs,
                   std::vector <DoubleFunct2D*> init_cond)
  : Example_NonStationary2D(exact, bc, bd, coeffs, timedependentrhs, 
                          timedependentcoeffs, init_cond)
  {

  };
};
#endif // _Example_TimeNSE2D_
