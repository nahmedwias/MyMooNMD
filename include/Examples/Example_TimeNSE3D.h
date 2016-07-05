/** ************************************************************************ 
*
* @class Example_TimeNSE3D
* @brief store all functions needed to describe a (Navier-)Stokes example 
* 
* The standard constructor of this class will fill the vectors (in Example2D) 
* with pointers to the functions needed to fully describe a particular example.
* 
* @author 
* @date   
* 
* @ruleof0
* 
 ************************************************************************  */

#ifndef _EXAMPLE_TimeNSE3D_
#define _EXAMPLE_TimeNSE3D_

#include <Example_NonStationary3D.h>
#include <functional>

class Time_NSE3D;

class Example_TimeNSE3D : public Example_NonStationary3D
{
public:
  /** @brief default constructor
   * 
   * This intializes a (Navier-)Stokes example in 2D. It is chosen according
   * to example_code.
   */
  Example_TimeNSE3D(int example_code,
                    const ParameterDatabase& user_input_parameter_db);
  
  /** @brief initialize your own example
   * 
   * Create an example with all vectors already defined.
   */
  Example_TimeNSE3D(std::vector <DoubleFunct3D*> exact,
                   std::vector <BoundCondFunct3D*> bc,
                   std::vector <BoundValueFunct3D*> bd, 
                   CoeffFct3D *coeffs,                    
                   bool timedependentrhs, bool timedependentcoeffs, 
                   std::vector <DoubleFunct3D*> init_cond)
  : Example_NonStationary3D(exact, bc, bd, coeffs,  timedependentrhs, 
                          timedependentcoeffs, init_cond) 
  {
  };
  
  /// Apply the function stored as post processing routine.
  void do_post_processing(Time_NSE3D& tnse3d) const;
  
  /// Return kinematic viscosity, if set.
  double get_nu() const;

  //Declaration of special member functions - rule of zero
    
  //! Default copy constructor. Performs deep copy.
  Example_TimeNSE3D(const Example_TimeNSE3D&) = default;
  
  //! Default move constructor.
  Example_TimeNSE3D(Example_TimeNSE3D&&) = default;
  
  //! Default copy assignment operator. Performs deep copy.
  Example_TimeNSE3D& operator=(const Example_TimeNSE3D&) = default;
  
  //! Default move assignment operator
  Example_TimeNSE3D& operator=(Example_TimeNSE3D&&) = default;
  
  //! Default destructor.
  ~Example_TimeNSE3D() = default;

private:
  /// Function doing the post processing for a stationary example.
  /// TODO @ULRICH: 
  // put NSE3D argument const as soon as FEFunctions can be copied properly!
  std::function<void(Time_NSE3D &)> post_processing_stat;
  /// TODO @ULRICH Function doing the post processing for a time dependent example.
  
  // Diffusion coefficient = kinematic viscosity. Should replace
  // former global parameter 1/RE_NR.
  double nu;
};
#endif // _EXAMPLE_TimeNSE2D_
