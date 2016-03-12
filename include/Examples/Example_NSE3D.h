/** ************************************************************************ 
*
* @class Example_NSE3D
* @brief store all functions needed to describe a (Navier-)Stokes example 
* 
* Depending on the value of TDatabase::ParamDB->EXAMPLE, the standard 
* constructor of this class will fill the vectors (in Example3D) with pointers
* to the functions needed to fully describe a particular example.
* 
* @date      02.06.2015
* 
* @ruleof0
* 
 ************************************************************************  */

#ifndef __EXAMPLE_NSE3D__
#define __EXAMPLE_NSE3D__

#include<Example3D.h>


class Example_NSE3D : public Example3D
{
  public:
    /** @brief default constructor
     * 
     * This intializes a (Navier-)Stokes example in 3D. It is chosen according
     * to TDatabase::ParamDB->EXAMPLE.
     */
    Example_NSE3D();
    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example_NSE3D(std::vector <DoubleFunct3D*> exact,
                  std::vector <BoundCondFunct3D*> bc,
                  std::vector <BoundValueFunct3D*> bd, CoeffFct3D *coeffs)
      : Example3D(exact, bc, bd, coeffs) {};
  
    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    Example_NSE3D(const Example_NSE3D&) = default;

    //! Default move constructor.
    Example_NSE3D(Example_NSE3D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    Example_NSE3D& operator=(const Example_NSE3D&) = default;

    //! Default move assignment operator
    Example_NSE3D& operator=(Example_NSE3D&&) = default;

    //! Default destructor.
    ~Example_NSE3D() = default;
};


#endif // __EXAMPLE_NSE3D__
