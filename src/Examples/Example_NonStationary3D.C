#include <Example_NonStationary3D.h>

Example_NonStationary3D::Example_NonStationary3D(const ParameterDatabase & db)
 : Example3D(db)
 , timeDependentRhs()
 , timeDependentCoeffs()
 , initialCondtion()
{

}

Example_NonStationary3D::Example_NonStationary3D(std::vector <DoubleFunct3D*> exact,
                          std::vector <BoundCondFunct3D*> bc,
                          std::vector <BoundValueFunct3D*> bd, 
                          CoeffFct3D coeffs,
                          bool timedependentrhs, bool timedependentcoeffs,
                          std::vector <DoubleFunct3D*> init_cond)
: Example3D(exact, bc, bd, coeffs)
  , timeDependentRhs(timedependentrhs)
  , timeDependentCoeffs(timedependentcoeffs)
  , initialCondtion(init_cond)
{

}
