#include <Example_NonStationary2D.h>

Example_NonStationary2D::Example_NonStationary2D() 
 : Example2D(), timeDependentRhs(true), 
   timeDependentCoeffs(true), initialCOndtion()
{
}

Example_NonStationary2D::Example_NonStationary2D(std::vector <DoubleFunct2D*> exact,
              std::vector <BoundCondFunct2D*> bc, std::vector <BoundValueFunct2D*> bd, 
              CoeffFct2D *coeffs, bool timedependentrhs, bool timedependentcoeffs, 
              std::vector <DoubleFunct2D*> init_cond)
 :Example2D(exact, bc, bd, coeffs)
 , timeDependentRhs(timedependentrhs) 
 , timeDependentCoeffs(timedependentcoeffs)
 , initialCOndtion(init_cond)
{
  
}
