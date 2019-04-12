/** ************************************************************************
 *
 * @class Example_Brinkman2D
 * @brief store all functions needed to describe a Brinkman2D example
 *
 * The standard constructor of this class will fill the vectors
 * (in Example2D) with pointers to the functions needed to fully
 * describe a particular example.
 *
 * @author    Alfonso Caiazzo & Laura Blank,
 * @date      18.07.2016
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
    explicit Example_Brinkman2D(
      const ParameterDatabase& user_input_parameter_db);
    
    /** @brief initialize your own example
     *
     * Create an example with all vectors already defined.
     */
    Example_Brinkman2D(const std::vector<DoubleFunct2D*>& exact,
                       const std::vector<BoundCondFunct2D*>& bc,
                       const std::vector<BoundValueFunct2D*>& bd,
                       const CoeffFct2D coeffs)
    : Example2D(exact, bc, bd, coeffs) {};
    
    /// Apply the function stored as post processing routine.
    void do_post_processing(Brinkman2D& brinkman2d) const;

    /// Get physical parameter from the database
    double get_viscosity() const;
    double get_effective_viscosity() const;
    double get_permeablity() const;
    
    /// Get residual-based equal-order stabilization weight, if set in the .dat-file.
    double get_stab() const;

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
    
    double stab_weight;

  private:
    /// Function doing the post processing for a stationary example.
    /// TODO put Brinkman2D argument const as soon as FEFunctions can be copied properly!
    std::function<void(Brinkman2D &)> post_processing_stat;
    /// TODO Function doing the post processing for a time dependent example.
};


#endif // __EXAMPLE_Brinkman2D__
