/** ************************************************************************ 
*
* @class     Example3D
* @brief     store all functions needed to describe an example
* 
* Mainly this stores the exact solution (if available), boundary condition,
* boundary data, and the coefficients of the pde. Note that you almost always
* want to create an object of a derived class, rather than of type Example3D 
* itself.
* 
* Below we use std::vector in case there are multiple solution components (e.g.
* velocity components and pressure).
* 
* @date      01.06.2015
 ************************************************************************  */

#ifndef __EXAMPLE3D__
#define __EXAMPLE3D__

#include <MooNMD_Io.h>
#include <Constants.h>
#include <vector>


class Example3D 
{
  protected:
    /**
    * @brief default constructor 
    * 
    * This is used only by the classes derived from this class.
    */
    Example3D();
  
  public:
    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example3D(std::vector <DoubleFunct3D*> exact,
              std::vector <BoundCondFunct3D*> bc,
              std::vector <BoundValueFunct3D*> bd, CoeffFct3D *coeffs);

    ~Example3D();

    /* functions representing the exact solution */
    std::vector <DoubleFunct3D*> exact_solution;
    /* functions representing the boundary conditions */
    std::vector <BoundCondFunct3D*> boundary_conditions;
    /* functions representing the boundary data */
    std::vector <BoundValueFunct3D*> boundary_data;
    /* functions representing the coefficients of the pde */
    CoeffFct3D *problem_coefficients;
    
    //void *example_info();

    /** getters */
    const std::vector <DoubleFunct3D*> & get_exact() const 
    { return exact_solution; }

    DoubleFunct3D* get_exact(unsigned int i) const
    { return exact_solution.at(i); }

    BoundCondFunct3D** get_bc()
    { return &boundary_conditions[0]; }

    BoundCondFunct3D* get_bc(unsigned int i) const
    { return boundary_conditions.at(i); }

    BoundValueFunct3D** get_bd()
    { return &boundary_data[0]; }

    BoundValueFunct3D* get_bd(unsigned int i) const
    { return boundary_data.at(i); }

    CoeffFct3D* get_coeffs() const
    { return problem_coefficients; }
};

#endif // __EXAMPLE3D__

