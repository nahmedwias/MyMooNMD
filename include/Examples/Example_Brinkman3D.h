/** ************************************************************************
 *
 * @class Example_Brinkman3D
 * @brief store all functions needed to describe a Brinkman example
 *
 * The standard constructor of this class will fill the vectors
 * (in Example3D) with pointers to the functions needed to fully
 * describe a particular example.
 *
 * @author    Alfonso Caiazzo & Laura Blank,
 * @date      06.10.2016
 *
 * @ruleof0
 *
 ************************************************************************  */

#ifndef __EXAMPLE_Brinkman3D__
#define __EXAMPLE_Brinkman3D__

#include<Example3D.h>
#include <functional>

class Brinkman3D; //forward declaration

class Example_Brinkman3D : public Example3D
{
public:
    /** @brief default constructor
     *
     * This intializes a Brinkman example in 3D. It is chosen according
     * to example_code.
     */
    explicit Example_Brinkman3D(
      const ParameterDatabase& user_input_parameter_db);
    
    /** @brief initialize your own example
     *
     * Create an example with all vectors already defined.
     */
    Example_Brinkman3D(const std::vector<DoubleFunct3D*>& exact,
                       const std::vector<BoundCondFunct3D*>& bc,
                       const std::vector<BoundValueFunct3D*>& bd,
                       const CoeffFct3D& coeffs)
    : Example3D(exact, bc, bd, coeffs) {};
    
    /// Apply the function stored as post processing routine.
    void do_post_processing(Brinkman3D& brinkman3d) const;

    /// Get physical parameter from the database
    double get_viscosity() const;
    double get_effective_viscosity() const;
    double get_permeablity() const;

    /// Get residual-based equal-order stabilization weight, if set in the .dat-file.
    double get_stab() const;

    //Declaration of special member functions - rule of zero
    
    //! Default copy constructor. Performs deep copy.
    Example_Brinkman3D(const Example_Brinkman3D&) = default;
    
    //! Default move constructor.
    Example_Brinkman3D(Example_Brinkman3D&&) = default;
    
    //! Default copy assignment operator. Performs deep copy.
    Example_Brinkman3D& operator=(const Example_Brinkman3D&) = default;
    
    //! Default move assignment operator
    Example_Brinkman3D& operator=(Example_Brinkman3D&&) = default;
    
    //! Default destructor.
    ~Example_Brinkman3D() = default;
    
    double stab_weight;

  private:
    /// Function doing the post processing for a stationary example.
    /// TODO put Brinkman3D argument const as soon as FEFunctions can be copied properly!
    std::function<void(Brinkman3D &)> post_processing_stat;
    /// TODO Function doing the post processing for a time dependent example.
};


#endif // __EXAMPLE_Brinkman3D__
