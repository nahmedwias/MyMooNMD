/** ************************************************************************ 
*
* @class Example_NSE2D
* @brief store all functions needed to describe a (Navier-)Stokes example 
* 
* The standard
* constructor of this class will fill the vectors (in Example2D) with pointers
* to the functions needed to fully describe a particular example.
* 
* @author    Ulrich Wilbrandt, 
* @date      14.03.15
* 
* @ruleof0
* 
 ************************************************************************  */

#ifndef __EXAMPLE_NSE2D__
#define __EXAMPLE_NSE2D__

#include<Example2D.h>


class Example_NSE2D : public Example2D
{
  public:
    /** @brief default constructor
     * 
     * This intializes a (Navier-)Stokes example in 2D. It is chosen according
     * to example_code.
     */
    Example_NSE2D(int example_code,
                  const ParameterDatabase& user_input_parameter_db);

    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example_NSE2D(std::vector <DoubleFunct2D*> exact,
                  std::vector <BoundCondFunct2D*> bc,
                  std::vector <BoundValueFunct2D*> bd, CoeffFct2D *coeffs)
      : Example2D(exact, bc, bd, coeffs) {};
  
    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    Example_NSE2D(const Example_NSE2D&) = default;

    //! Default move constructor.
    Example_NSE2D(Example_NSE2D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    Example_NSE2D& operator=(const Example_NSE2D&) = default;

    //! Default move assignment operator
    Example_NSE2D& operator=(Example_NSE2D&&) = default;

    //! Default destructor.
    ~Example_NSE2D() = default;
};


#endif // __EXAMPLE_NSE2D__
