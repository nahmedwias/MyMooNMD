/** ************************************************************************ 
*
* @class Example_TimeLinElastic2D
* @brief store all functions needed to describe a Linear elasticity example
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
#ifndef _Example_TimeLinElastic2D_
#define _Example_TimeLinElastic2D_

#include <Example_NonStationary2D.h>
#include <functional>

class Time_LinElastic2D; //forward declaration

class Example_TimeLinElastic2D : public Example_NonStationary2D
{
public:
  /** @brief default constructor
   * This intializes a Time dependent (Navier-)Stokes example in 2D. 
   * It is chosen according to example_code.
   */
  Example_TimeLinElastic2D(const ParameterDatabase& user_input_parameter_db);

  /// Apply the function stored as post processing routine.
  void do_post_processing(Time_LinElastic2D& tlinelastic2d) const;

  private:
  /// Function doing the post processing for a stationary example.
  /// TODO put Time_LinearElastic2D argument const as soon as FEFunctions can be copied properly!
  std::function<void(Time_LinElastic2D &)> post_processing_stat;
  /// TODO Function doing the post processing for a time dependent example.

};
#endif // _Example_TimeLinElastic2D_
