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
#include <functional>

class Time_NSE2D; //forward declaration

class Example_TimeNSE2D : public Example_NonStationary2D
{
public:
  /** @brief default constructor
   * This intializes a Time dependent (Navier-)Stokes example in 2D. 
   * It is chosen according to example_code.
   */
  Example_TimeNSE2D(const ParameterDatabase& user_input_parameter_db);
  /** @brief initialize your own example
   * 
   * Create an example with all vectors already defined.
   */
  Example_TimeNSE2D(std::vector <DoubleFunct2D*> exact,
                   std::vector <BoundCondFunct2D*> bc,
                   std::vector <BoundValueFunct2D*> bd, 
                   CoeffFct2D coeffs,
                   bool timedependentrhs, bool timedependentcoeffs,
                   std::vector <DoubleFunct2D*> init_cond);

  /// Apply the function stored as post processing routine.
  void do_post_processing(Time_NSE2D& tnse2d, double& val) const;

  /// Return kinematic viscosity, if set.
  double get_nu() const;

  private:
  /// Function doing the post processing for a stationary example.
  /// TODO put Time_NSE2D argument const as soon as FEFunctions can be copied properly!
  std::function<void(Time_NSE2D&, double& val)> post_processing_stat;
  
  std::function<void(Time_NSE2D &)> post_processing_stat_old;
  /// TODO Function doing the post processing for a time dependent example.

};
#endif // _Example_TimeNSE2D_
