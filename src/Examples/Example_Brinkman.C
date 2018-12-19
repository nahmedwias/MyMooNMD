#include <cmath> // e.g. sin
#include <Example.h>
#include <MooNMD_Io.h>

Example Example::Brinkman(const ParameterDatabase& param_db)
{
  int example = param_db["example"];
  double reynolds_number = param_db["reynolds_number"];
  /// @todo include Brinkman parameters into database and use them here
  double viscosity = 1.; // param_db["viscosity"];
  double effective_viscosity = 1.; // param_db["effective_viscosity"];
  double permeability = 1.; // param_db["permeability"];
  // flag, false means the laplace term is discretized with the gradient, true
  // means it is discretized with the deformation tensor
  bool deformation_tensor = false;

  // indicate if two or three space dimensions
  bool two_d = true;
#ifdef __3D__
  two_d = false;
#endif

  // function needed to create the PDECoefficients
  std::function<void(const Point&, double, std::vector<double>&)> f;
  // vectors which will contain exactly one entry each.
  std::vector<BoundaryCondition> bc;     // boundary condition
  std::vector<BoundaryData> bd;          // boundary data
  std::vector<AnalyticalFunction> exact; // exact solution (possibly unknown)
  std::string name;                      // to print some information later on

  switch(example)
  {
    case 0: // Poisseuille flow
    {
      if(!two_d)
        ErrThrow("Poiseuille flow example in 3D not yet implemented");
      // boundary condition:
      bc.push_back(BoundaryCondition(DIRICHLET));
      auto u1 = [](const Point& point, FunctionEvaluation& v)
      {
        const double y = point.y();
        v.set(4 * y * (1 - y), 0, MultiIndex::D00);
        v.set(0., 0, MultiIndex::D10);
        v.set(4 - 8 * y, 0, MultiIndex::D01);
        v.set(0., 0, MultiIndex::D20);
        v.set(0., 0, MultiIndex::D11);
        v.set(-8., 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(u1));
      // exact solution, second velocity component (all zero)
      exact.push_back(AnalyticalFunction(0.0));
      // exact solution, pressure component
      auto pr = [reynolds_number](const Point& point, FunctionEvaluation& v)
      {
        const double x = point.x();
        v.set(8 * reynolds_number * (0.5 - x), 0, MultiIndex::D00);
        v.set(-8 * reynolds_number, 0, MultiIndex::D10);
        v.set(0., 0, MultiIndex::D01);
        v.set(0., 0, MultiIndex::D20);
        v.set(0., 0, MultiIndex::D11);
        v.set(0., 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(pr));
      // boundary data, automatically computed from known solution
      // bd.push_back(BoundaryData(bc, exact));
      // Boundary data for first velocity component
      auto g1 =
          [deformation_tensor](unsigned int component, double t, double)
      {
        return (component == 1 || component == 3) ? 4 * t * (1 - t) : 0.;
      };
      bd.push_back(BoundaryData(g1));
      // Boundary data for first velocity component (all zero)
      bd.push_back(BoundaryData(0.0));
      bd.push_back(BoundaryData(0.0)); // all zero for pressure
      // the coefficient function
      f = [reynolds_number, viscosity, effective_viscosity, permeability](
          const Point&, double, std::vector<double>& coeffs)
      {
        coeffs[0] = 1. / reynolds_number; // reynolds number
        coeffs[1] = 0.;            // right hand side, first component
        coeffs[2] = 0.;            // right hand side, second component
        coeffs[3] = 0.;            // divergence
        coeffs[4] = viscosity;
        coeffs[5] = effective_viscosity;
        coeffs[6] = permeability;
      };
      break;
    }
    case 1:
    {
      ErrThrow("Brinkman example 1: to be implemented");
    }
    case 2: // driven cavity
    {
      ErrThrow("Brinkman example 2: to be implemented");
    }
    case 3:
    {
      ErrThrow("Brinkman example 3: to be implemented");
      break;
    }
    default:
      ErrThrow("unknown example ", example);
      break;
  }

  PDECoefficients coeffs(f, Problem_type::NavierStokes, false, false);

  Output::info<2>("Example", "using the Navier--Stokes example: ", name);
  Output::dash<2>("The Reynolds number is ", reynolds_number);
  return Example(std::move(bc), std::move(bd), std::move(coeffs),
                 std::move(exact));
}
