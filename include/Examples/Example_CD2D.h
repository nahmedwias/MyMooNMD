/** ************************************************************************ 
*
* @class Example_CD2D
* @brief store all functions needed to describe a convection--diffusion example 
* 
* Depending on the value of TDatabase::ParamDB->EXAMPLE, the standard 
* constructor of this class will fill the vectors (in Example2D) with pointers
* to the functions needed to fully describe a particular example.
* 
* So far the following examples are implemented and enabled:
* Stationary examples
* 	0 - Simple Sine Laplace (analytical example)
* 	1 - Two interior Layers  (Knopp, Lube, Rapin, CMAME 2002)
* 	2 - Hemker Rotating Body (Hemker 1996)
* 	3 - Sharp Boundary Layer (Kuzmin & Moeller 2005)
*
* Time dependent examples
* 	101 - Exponential Function (analytical example)
* 	102 - Sine Sine Sine (analytical example)
* 	103 - Sine Cosine (analytical example)
*
*
* @author    Ulrich Wilbrandt, 
* @date      13.03.15
* 
* @ruleof0
* 
 ************************************************************************  */


#ifndef __EXAMPLE_CD2D__
#define __EXAMPLE_CD2D__

#include<Example2D.h>


class Example_CD2D : public Example2D 
{
  public:
    /** @brief default constructor
     * 
     * This intializes a convection-diffusion example in 2D. It is chosen 
     * according to TDatabase::ParamDB->EXAMPLE.
     */
    Example_CD2D();
    
    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example_CD2D(std::vector <DoubleFunct2D*> exact,
                  std::vector <BoundCondFunct2D*> bc,
                  std::vector <BoundValueFunct2D*> bd,
                  CoeffFct2D *coeffs,
                  std::vector <DoubleFunct2D*> init_cond =
                      std::vector<DoubleFunct2D*>())
    : Example2D(exact, bc, bd, coeffs, init_cond) {};

    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    Example_CD2D(const Example_CD2D&) = default;

    //! Default move constructor.
    Example_CD2D(Example_CD2D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    Example_CD2D& operator=(const Example_CD2D&) = default;

    //! Default move assignment operator
    Example_CD2D& operator=(Example_CD2D&&) = default;

    //! Default destructor.
    ~Example_CD2D() = default;
};


#endif // __EXAMPLE_CD2D__
