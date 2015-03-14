/** ************************************************************************ 
*
* @class Example_CD2D
* @brief store all functions needed to describe a convection--diffusion example 
* 
* Depending on the value of TDatabase::ParamDB->EXAMPLE, the standard 
* constructor of this class will fill the vectors (in Example2D) with pointers
* to the functions needed to fully describe a particular example.
* 
* @author    Ulrich Wilbrandt, 
* @date      13.03.15
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
                  std::vector <BoundValueFunct2D*> bd, CoeffFct2D *coeffs)
    : Example2D(exact, bc, bd, coeffs) {};

    ~Example_CD2D(){};
};


#endif // __EXAMPLE_CD2D__