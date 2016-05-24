/** ************************************************************************ 
*
* @class Example_NSE3D
* @brief store all functions needed to describe a (Navier-)Stokes example 
* 
* The standard
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
#include <functional>

class NSE3D; //forward declaration

class Example_NSE3D : public Example3D
{
  public:
    /** @brief default constructor
     * 
     * This intializes a (Navier-)Stokes example in 3D. It is chosen according
     * to example_code.
     */
    Example_NSE3D(int example_code);
    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example_NSE3D(std::vector <DoubleFunct3D*> exact,
                  std::vector <BoundCondFunct3D*> bc,
                  std::vector <BoundValueFunct3D*> bd, CoeffFct3D *coeffs)
      : Example3D(exact, bc, bd, coeffs) {};

    /// Apply the function stored as post processing routine.
    void do_post_processing(NSE3D& nse3d) const;

    /// Return kinematic viscosity, if set.
    double get_nu() const;
  
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

  private:
    /// Function doing the post processing for a stationary example.
    /// TODO put NSE3D argument const as soon as FEFunctions can be copied properly!
    std::function<void(NSE3D &)> post_processing_stat;
    /// TODO Function doing the post processing for a time dependent example.

    // Diffusion coefficient = kinematic viscosity. Should replace
    // former global parameter 1/RE_NR.
    double nu;
};


#endif // __EXAMPLE_NSE3D__
