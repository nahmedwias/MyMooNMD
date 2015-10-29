#include <Example3D.h>

Example3D::Example3D() 
 : exact_solution(), boundary_conditions(), boundary_data(), 
   problem_coefficients(NULL)
{
}

Example3D::Example3D(std::vector <DoubleFunct3D*> exact,
                     std::vector <BoundCondFunct3D*> bc,
                     std::vector <BoundValueFunct3D*> bd, CoeffFct3D *coeffs)
 : exact_solution(exact), boundary_conditions(bc), boundary_data(bd),
   problem_coefficients(coeffs)
{ 
}
