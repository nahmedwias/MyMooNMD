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
  
    ~Example_NSE3D(){};
};


#endif // __EXAMPLE_NSE3D__