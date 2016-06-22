/** ************************************************************************ 
*
* @class Example_TimeCD3D
* @brief store all functions needed to describe a convection--diffusion example 
* 
* Depending on the value of TDatabase::ParamDB->EXAMPLE, the standard 
* constructor of this class will fill the vectors (in Example2D_TimeProblem2D) with pointers
* to the functions needed to fully describe a particular example.
* 
* So far the following examples are implemented and enabled:
* 	0 - 
* 	1 - 
* 	2 - 
*
*
* @author  
* @date    
* 
* @ruleof0
* 
 ************************************************************************  */
#ifndef _Example_TimeCD3D_
#define _Example_TimeCD3D_

#include <Example_NonStationary3D.h>

class Example_TimeCD3D : public Example_NonStationary3D
{
public:
    /**
   * @brief default constructor
   * This intializes a convection-diffusion example in 2D. It is chosen 
   * according to example_code.
   */
  Example_TimeCD3D(int example_code,
                   const ParameterDatabase& user_input_parameter_db);
  
  /** @brief initialize your own example
   * 
   * Create an example with all vectors already defined.
   */
  Example_TimeCD3D(std::vector <DoubleFunct3D*> exact,
                   std::vector <BoundCondFunct3D*> bc,
                   std::vector <BoundValueFunct3D*> bd, 
                   CoeffFct3D *coeffs,                    
                   bool timedependentrhs, bool timedependentcoeffs, 
                   std::vector <DoubleFunct3D*> init_cond)
  : Example_NonStationary3D(exact, bc, bd, coeffs,  timedependentrhs, 
                          timedependentcoeffs, init_cond) {};
  
  //Declaration of special member functions - rule of zero
  //! Default copy constructor. Performs deep copy.
  Example_TimeCD3D(const Example_TimeCD3D&) = default;
  
  //! Default move constructor.
  Example_TimeCD3D(Example_TimeCD3D&&) = default;
  
  //! Default copy assignment operator. Performs deep copy.
  Example_TimeCD3D& operator=(const Example_TimeCD3D&) = default;
  
  //! Default move assignment operator
  Example_TimeCD3D& operator=(Example_TimeCD3D&&) = default;
  
  //! Default destructor.
  ~Example_TimeCD3D() = default;  
};
#endif // _Example_TimeCD3D_
