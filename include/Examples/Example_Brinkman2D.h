/** ************************************************************************ 
*
* @class Example_Brinkman2D
* @brief store all functions needed to describe a Brinkman example
* 
* The standard constructor of this class will fill the vectors
* (in Example2D) with pointers to the functions needed to fully
* describe a particular example.
* 
* @author    Ulrich Wilbrandt, 
* @date      14.03.15
* 
* @ruleof0
* 
 ************************************************************************  */

#ifndef __EXAMPLE_Brinkman2D__
#define __EXAMPLE_Brinkman2D__

#include<Example2D.h>
#include <functional>

class Brinkman2D; //forward declaration

class Example_Brinkman2D : public Example2D
{
  public:
    /** @brief default constructor
     * 
     * This intializes a Brinkman example in 2D. It is chosen according
     * to example_code.
     */
    Example_Brinkman2D(int example_code,
                       const ParameterDatabase& user_input_parameter_db);

    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example_Brinkman2D(std::vector <DoubleFunct2D*> exact,
                  std::vector <BoundCondFunct2D*> bc,
                  std::vector <BoundValueFunct2D*> bd, CoeffFct2D *coeffs)
      : Example2D(exact, bc, bd, coeffs) {};
  
    /// Apply the function stored as post processing routine.
    void do_post_processing(Brinkman2D& brinkman2d) const;

    /// Return kinematic viscosity, if set.
    double get_nu() const;

    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    Example_Brinkman2D(const Example_Brinkman2D&) = default;

    //! Default move constructor.
    Example_Brinkman2D(Example_Brinkman2D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    Example_Brinkman2D& operator=(const Example_Brinkman2D&) = default;

    //! Default move assignment operator
    Example_Brinkman2D& operator=(Example_Brinkman2D&&) = default;

    //! Default destructor.
    ~Example_Brinkman2D() = default;

  private:
    /// Function doing the post processing for a stationary example.
    /// TODO put Brinkman2D argument const as soon as FEFunctions can be copied properly!
    std::function<void(Brinkman2D &)> post_processing_stat;
    /// TODO Function doing the post processing for a time dependent example.

    // Diffusion coefficient = kinematic viscosity. Should replace
    // former global parameter 1/RE_NR.
    double nu;
};


#endif // __EXAMPLE_Brinkman2D__
