#include <Example2D.h>

Example2D::Example2D() 
 : exact_solution(), boundary_conditions(), boundary_data(), 
   problem_coefficients(NULL)
{
}

Example2D::Example2D(std::vector <DoubleFunct2D*> exact,
                     std::vector <BoundCondFunct2D*> bc,
                     std::vector <BoundValueFunct2D*> bd, CoeffFct2D *coeffs)
 : exact_solution(exact), boundary_conditions(bc), boundary_data(bd),
   problem_coefficients(coeffs)
{ 
}
