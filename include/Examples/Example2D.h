/** ************************************************************************ 
*
* @class     Example2D
* @brief     store all functions needed to describe an example
* 
* Mainly this stores the exact solution (if available), boundary condition,
* boundary data, and the coefficients of the pde. Note that you almost always
* want to create an object of a derived class, rather than of type Example2D 
* itself.
* 
* Below we use std::vector in case there are multiple solution components (e.g.
* velocity components and pressure).
* 
* @author    Ulrich Wilbrandt, 
* @date      13.03.15
* 
* @ruleof0
* 
 ************************************************************************  */

#ifndef __EXAMPLE2D__
#define __EXAMPLE2D__

#include <Constants.h>
#include <vector>


class Example2D 
{
  protected:
    /**
    * @brief default constructor 
    * 
    * This is used only by the classes derived from this class.
    */
    Example2D();
  
  public:
    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example2D(std::vector <DoubleFunct2D*> exact,
              std::vector <BoundCondFunct2D*> bc,
              std::vector <BoundValueFunct2D*> bd, 
              CoeffFct2D *coeffs);

    /* functions representing the exact solution */
    std::vector <DoubleFunct2D*> exact_solution;
    /* functions representing the boundary conditions */
    std::vector <BoundCondFunct2D*> boundary_conditions;
    /* functions representing the boundary data */
    std::vector <BoundValueFunct2D*> boundary_data;
    /* functions representing the coefficients of the pde */
    CoeffFct2D *problem_coefficients;
    //void *example_info();

    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    Example2D(const Example2D&) = default;

    //! Default move constructor.
    Example2D(Example2D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    Example2D& operator=(const Example2D&) = default;

    //! Default move assignment operator
    Example2D& operator=(Example2D&&) = default;

    //! Default destructor.
    virtual ~Example2D() = default;


    // Getter functions

    const std::vector <DoubleFunct2D*> & get_exact() const 
    { return exact_solution; }

    DoubleFunct2D* get_exact(unsigned int i) const
    { return exact_solution.at(i); }

    BoundCondFunct2D * const * get_bc() const
    { return &boundary_conditions[0]; }

    BoundCondFunct2D* get_bc(unsigned int i) const
    { return boundary_conditions.at(i); }

    BoundValueFunct2D * const * get_bd() const
    { return &boundary_data[0]; }

    BoundValueFunct2D* get_bd(unsigned int i) const
    { return boundary_data.at(i); }

    virtual CoeffFct2D* get_coeffs() const;
};

#endif // __EXAMPLE2D__

