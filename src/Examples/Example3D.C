#include <Example3D.h>


ParameterDatabase Example3D::default_example_database()
{
  Output::print<5>("creating a default Example3D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NSE3D database as well.
  ParameterDatabase db("Example3D parameter database");

  db.add("example", 0,
      "Choose which example to run. \nNote that depending on the type of "
      "problem you want to solve, different values are meaningful here. See "
      "the derived classes of 'Example3D'.", -5, 200);

  /** TDatabase::ParamDB->RE_NR */
  db.add("reynolds_number", 1.,
      "Reynolds number: dimensionless number which describes how viscous the  "
      "flow is (laminar, turbulent). Re = U.L/nu, where nu is the kinematic "
      "viscosity (=mu/rho)."
      "The higher it is, the more turbulent "
      "the flow is. Reynolds number can often be in the order of "
      "magnitude of millions. Then, the classical NSE is not relevant anymore"
      "to describe the flow. One should in these cases use turbulence Models,"
      "like Smagorinsky or k-eps. Maximum value is 1000,"
      "which already corresponds to a turbulent flow. Default value is 1."
      "Note that this is also equal to 1/eps, where eps is the DIMENSIONLESS "
      "viscosity (sometimes misleadingly named as nu) "
      "which sometimes appears in the dimensionless Navier-Stokes equations"
      "instead of 1/Re.",
      0., 1000.);

  /** TDatabase::ParamDB->PE_NR */
  db.add("diffusion_coefficient", 1.,
      "Diffusion coefficient: a factor in front of the diffusion term.",
      0., 1.);

  /** TDatabase::ParamDB->PERMEABILITY */
  db.add("permeability", (double) 1.,
      "permeability coefficient: a factor in front of the resistance term.",
      (double) 0., (double) 1000000.);
  /** TDatabase::ParamDB->VISCOSITY */
  db.add("viscosity", (double) 1.,
      "viscosity coefficient: a factor in front of the Laplacian or the resistance term.",
      (double) 0., (double) 1000000.);
  /** TDatabase::ParamDB->EFFECTIVE_VISCOSITY */
  db.add("effective_viscosity", (double) 1.,
      "effective_viscosity coefficient: a factor in front of the Laplacian term.",
      (double) 0., (double) 1000000.);


  // NEW LB 10.10.18
  db.add("inverse_permeability", (double) 0.,
      "viscosity/permeability: a factor in front of the resistance term.",
      (double) 0., (double) 1000000.);

  ///@todo boundary-related parameters could be moved soon to a dedicated class
  // Neumann BC
  db.add("n_neumann_bd",0u,
   "Number of boundary components where Neumann BC are imposed"
   " This variable is used as a preliminary flag. The number of Neumann BC"
   " must match the length of the parameter neumann_id.",0u,9u);

  db.add("neumann_id", {1u,1u,1u,1u,1u,1u,1u,1u,1u,1u},
   "Component ID of Neumann boundaries. On these edges the terms (hat p,v.n)_E"
   " will be assembled. The values -hat p- are given via the parameters "
   " neumann_values. "
   " Attention: at the moment max. 10 boundaries are allowed. This can be "
   " easily changed in the file " __FILE__,0u,999u);

  db.add("neumann_value", {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
   " Neumann values to be assigned to the respetive boundary."
   " Warning: the length must be the same as the one given in neumann_boundary_id");

  // Nitsche BC
  db.add("n_nitsche_bd",0u,
   " Number of boundary components where Nitsche BC are imposed"
   " This variable is used as a preliminary flag. The number of Nitsche BC"
   " must match the length of the parameter nitsche_id",0u,9u);

  db.add("nitsche_id", {1u,1u,1u,1u,1u,1u,1u,1u,1u,1u},
   "Component ID of boundaries where the Nitsche method is used to "
   " impose (weak) Dirichlet BC."
   " Attention: at the moment max. 10 boundaries are allowed. This can be "
   " easily changed in the file " __FILE__,0u,999u);

  db.add("nitsche_penalty", {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
   " Nitsche penalty (gamma) to be imposed on each boundary");


  db.add("symmetric_nitsche_u",-1,
   " Coefficient for the term -(u,mu dv/dn) "
   " Symmetric version (symmetric_nitsche_u = 1): assemble -(mu du/dn,v)-(u,mu dv/dn), "
   " non-sym (symmetric_nitsche_u = -1): assemble -(mu du/dn,v)+(u,mu dv/dn)");

  db.add("symmetric_nitsche_p",-1,
   " Coefficient for the term (u.n,q) "
   " Symmetric version (symmetric_nitsche_p = 1): assemble (p,v.n)+(u.n,q), "
   " non-sym (symmetric_nitsche_u = -1): assemble (p,v.n)-(u.n,q)");

  return db;
}

  Example3D::Example3D(const ParameterDatabase & db)
: example_database(Example3D::default_example_database()), exact_solution(),
  boundary_conditions(), boundary_data(), problem_coefficients(nullptr)
{
  this->example_database.merge(db, false);
}

Example3D::Example3D(std::vector <DoubleFunct3D*> exact,
    std::vector <BoundCondFunct3D*> bc,
    std::vector <BoundValueFunct3D*> bd, CoeffFct3D coeffs)
: example_database(Example3D::default_example_database()), exact_solution(exact),
  boundary_conditions(bc), boundary_data(bd), problem_coefficients(coeffs)
{ 
}
