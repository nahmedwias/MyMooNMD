/** ************************************************************************ 
*
* @class Example_CD3D
* @brief store all functions needed to describe a convection--diffusion example 
* 
* Depending on the value of TDatabase::ParamDB->EXAMPLE, the standard 
* constructor of this class will fill the vectors (in Example3D) with pointers
* to the functions needed to fully describe a particular example.
* 
* @author    Ulrich Wilbrandt, 
* @date      09.06.15
 ************************************************************************  */


#ifndef __EXAMPLE_CD3D__
#define __EXAMPLE_CD3D__

#include<Example3D.h>


class Example_CD3D : public Example3D 
{
  public:
    /** @brief default constructor
     * 
     * This intializes a convection-diffusion example in 3D. It is chosen 
     * according to TDatabase::ParamDB->EXAMPLE.
     */
    Example_CD3D();
    
    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example_CD3D(std::vector <DoubleFunct3D*> exact,
                 std::vector <BoundCondFunct3D*> bc,
                 std::vector <BoundValueFunct3D*> bd, CoeffFct3D *coeffs)
    : Example3D(exact, bc, bd, coeffs) {};

    ~Example_CD3D(){};
};


#endif // __EXAMPLE_CD3D__