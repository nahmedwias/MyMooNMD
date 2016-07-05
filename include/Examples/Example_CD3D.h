/** ************************************************************************ 
*
* @class Example_CD3D
* @brief store all functions needed to describe a convection--diffusion example 
* 
* Depending on an integer example code, the standard
* constructor of this class will fill the vectors (in Example3D) with pointers
* to the functions needed to fully describe a particular example.
* 
* @author    Ulrich Wilbrandt, 
* @date      09.06.15
* 
* @ruleof0
* 
 ************************************************************************  */


#ifndef __EXAMPLE_CD3D__
#define __EXAMPLE_CD3D__

#include<Example3D.h>


class Example_CD3D : public Example3D 
{
  public:
    /** @brief default constructor
     * 
     * This intializes a convection-diffusion example in 3D.
     * It is chosen according to the given example_code.
     */
    Example_CD3D(int example_code,
                 const ParameterDatabase& user_input_parameter_db);
    
    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example_CD3D(std::vector <DoubleFunct3D*> exact,
                 std::vector <BoundCondFunct3D*> bc,
                 std::vector <BoundValueFunct3D*> bd, CoeffFct3D *coeffs)
    : Example3D(exact, bc, bd, coeffs) {};


    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    Example_CD3D(const Example_CD3D&) = default;

    //! Default move constructor.
    Example_CD3D(Example_CD3D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    Example_CD3D& operator=(const Example_CD3D&) = default;

    //! Default move assignment operator
    Example_CD3D& operator=(Example_CD3D&&) = default;

    //! Default destructor.
    ~Example_CD3D() = default;
};


#endif // __EXAMPLE_CD3D__
