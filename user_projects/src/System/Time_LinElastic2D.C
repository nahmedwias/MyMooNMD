#include <Time_LinElastic2D.h>
#include <Database.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include "../../../include/AssembleRoutines/Assemble2D.h"
#include "../../../include/AssembleRoutines/LocalAssembling2D.h"

/**************************************************************************** */
ParameterDatabase get_default_Time_LinElastic2D_parameters()
{
  Output::print<5>("creating a default Time_LinElastic2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TLinElastic2D parameter database");

  //set up a nonlinit_database and merge
  ParameterDatabase nl_db = ParameterDatabase::default_nonlinit_database();
  db.merge(nl_db,true);

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  // a default time database
  ParameterDatabase time_db = ParameterDatabase::default_time_database();
  db.merge(time_db,true);

  // add parameters specific to LinElastic2D
  db.add("lamecoeff_lambda", 0.5,
         "Lame coefficient lambda.", 0.,1000.);
  db.add("lamecoeff_mu", 0.5,
         "Lame coefficient mu.", 0., 1000.);
  db.add("rho_material", 7800.,
         "Density of material.", 0., 10000. );

  return db;
}

/**************************************************************************** */
Time_LinElastic2D::System_per_grid::System_per_grid(const Example_TimeLinElastic2D& example,
                                                    TCollection& coll)
: fe_space_(&coll, (char*)"space", (char*)"time_linelastic2d space", example.get_bc(0),
            TDatabase::ParamDB->ANSATZ_ORDER, nullptr),
  stiffness_matrix_({&fe_space_,&fe_space_}),
  mass_matrix_({&fe_space_,&fe_space_}),
  rhs_(this->stiffness_matrix_,true),
  solution_(this->stiffness_matrix_,false),
  u_(&fe_space_,(char*)"u", (char*)"u", solution_.block(0),
     solution_.length(0), 2),
  velocity_(this->stiffness_matrix_,false),
  v_(&fe_space_,(char*)"v", (char*)"v", velocity_.block(0),
     velocity_.length(0), 2),
  acceleration_(this->stiffness_matrix_,false),
  a_(&fe_space_,(char*)"a", (char*)"a", acceleration_.block(0),
     acceleration_.length(0), 2),
  vector_lambda_(solution_),
  lambda_(&fe_space_,(char*)"lambda", (char*)"lambda",
          vector_lambda_.block(0), vector_lambda_.length(0)),
  vector_mu_(solution_),
  mu_(&fe_space_,(char*)"mu", (char*)"mu",
       vector_mu_.block(0), vector_mu_.length(0)),
 vector_rho_(solution_),
 rho_(&fe_space_,(char*)"rho", (char*)"rho",
            vector_rho_.block(0), vector_rho_.length(0))

{
  // K has 4 blocks, M has 2 diagonal blocks
  stiffness_matrix_ = BlockFEMatrix::LinElastic2D(fe_space_);
  mass_matrix_ = BlockFEMatrix::Mass_LinElastic2D(fe_space_);

  cout << "CONSTRUCTOR OF SYSTEM PER GRID OK!!" << endl;
}

/**************************************************************************** */
Time_LinElastic2D::Time_LinElastic2D(const TDomain& domain,
                                     const ParameterDatabase& param_db,
                                     int reference_id)
: Time_LinElastic2D(domain, param_db, Example_TimeLinElastic2D(param_db), reference_id)
{
  cout << "CONSTRUCTOR1 OK!" << endl;
}

/**************************************************************************** */
Time_LinElastic2D::Time_LinElastic2D(const TDomain& domain,
                                     const ParameterDatabase& param_db,
                                     const Example_TimeLinElastic2D& ex,
                                     int reference_id)
: db_(get_default_Time_LinElastic2D_parameters()),
  example_(ex), solver_(param_db), outputwriter_(param_db),
  systems_()
{
  db_.merge(param_db,false);

  /* Construction of the systems per Grid */
  bool usingMultigrid = this->solver_.is_using_multigrid();
  if(!usingMultigrid)
  {
    // create the collection of cells from the domain (finest grid)
    TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
    // create finite element space and function, a matrix, rhs, and solution
    this->systems_.emplace_back(this->example_, *coll);

    // Initialization of Lame coefficients
    this->systems_.front().vector_lambda_= db_["lamecoeff_lambda"];
    this->systems_.front().vector_mu_    = db_["lamecoeff_mu"];
    this->systems_.front().vector_rho_   = db_["rho_material"];

    // interpolate initial displacements
    TFEFunction2D * u1 = this->systems_.front().u_.GetComponent(0);
    TFEFunction2D * u2 = this->systems_.front().u_.GetComponent(1);
    TFEFunction2D * v1 = this->systems_.front().v_.GetComponent(0);
    TFEFunction2D * v2 = this->systems_.front().v_.GetComponent(1);
    u1->Interpolate(example_.get_initial_cond(0));
    u2->Interpolate(example_.get_initial_cond(1));
    v1->Interpolate(example_.get_initial_cond(2));
    v2->Interpolate(example_.get_initial_cond(3));
    // add the fe function to the output object.
    outputwriter_.add_fe_vector_function(&this->systems_.front().u_);
    outputwriter_.add_fe_vector_function(&this->systems_.front().v_);
  }
  else
  {
    ErrThrow("MULTIGRID NOT IMPLEMENTED YET!!!");
  }

  // print out the information (cells, dofs, etc)
  this->output_problem_size_info();

  //  this->systems_.front().u_.GetComponent(0)->WriteSol("./","initial_values");
  cout << "CONSTRUCTOR2 OK!" << endl;
}

/**************************************************************************** */
void Time_LinElastic2D::assemble_initial_time()
{
  for(auto &s : this->systems_)
  {
    /* First, Local assembling */
    /* the first two fe_functions are u1, u2, just in case..
     * it is not used normally in the local assembling object
     * the three other fe_functions are lambda, mu  and rho
     * functions. They can be constant or variable
     * */
    TFEFunction2D * fe_functions[5] = {s.u_.GetComponent(0),
                                       s.u_.GetComponent(1),
                                       &s.lambda_, &s.mu_, &s.rho_};

    LocalAssembling2D local_assembling(TLinElastic2D_WholeSystem,
                                       fe_functions,
                                       this->example_.get_coeffs());


    /* Second, input information for Assemble2D */
    size_t n_fe_space = 1;
    const TFESpace2D * space = &s.fe_space_;
    size_t n_square_matrices = 6;
    TSquareMatrix2D *sqMatrices[6]{nullptr}; // maximum number of square matrices
    size_t nRhs = 2;
    double *RHSs[2] = {s.rhs_.block(0), s.rhs_.block(1)};
    const TFESpace2D *fe_rhs[2] = {space, space};

    BoundCondFunct2D * boundary_conditions[2] = {
      space->GetBoundCondition(), space->GetBoundCondition() };
    std::array<BoundValueFunct2D*, 2> non_const_bound_values;
    non_const_bound_values[0] = example_.get_bd()[0];
    non_const_bound_values[1] = example_.get_bd()[1];

    //set matrices
     std::vector<std::shared_ptr<FEMatrix>> stiffness_blocks =
         s.stiffness_matrix_.get_blocks_uniquely();
     std::vector<std::shared_ptr<FEMatrix>> mass_blocks =
         s.mass_matrix_.get_blocks_uniquely();
     sqMatrices[0] =reinterpret_cast<TSquareMatrix2D*>(stiffness_blocks.at(0).get());//K11
     sqMatrices[1] =reinterpret_cast<TSquareMatrix2D*>(stiffness_blocks.at(1).get());//K12
     sqMatrices[2] =reinterpret_cast<TSquareMatrix2D*>(stiffness_blocks.at(2).get());//K21
     sqMatrices[3] =reinterpret_cast<TSquareMatrix2D*>(stiffness_blocks.at(3).get());//K22
     sqMatrices[4] =reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());//M11
     sqMatrices[5] =reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(1).get());//M22

     /* Three, reset or preprocess for Assemble2D */
     s.rhs_.reset();
//sqMatrices[0]->reset(); //sqMatrices[1]->reset(); sqMatrices[2]->reset();
//     sqMatrices[3]->reset(); sqMatrices[4]->reset(); sqMatrices[5]->reset();

     /* Four, Assemble2D */
     Assemble2D(n_fe_space, &space, n_square_matrices,
                sqMatrices, 0, NULL, nRhs, RHSs,
                fe_rhs, boundary_conditions,
                non_const_bound_values.data(), local_assembling);
  }
//  this->get_stiff_matrix().get_blocks().at(0)->Print("A11");
//  this->get_stiff_matrix().get_blocks().at(1)->Print("A12");
//  this->get_stiff_matrix().get_blocks().at(2)->Print("A21");
//  this->get_stiff_matrix().get_blocks().at(3)->Print("A22");
//  this->get_stiff_matrix().print_matrix_info("stiffness");
//  this->get_mass_matrix().print_matrix_info("mass");
  cout << "END OF ASSEMBLE INITIAL TIME" << endl;
}

/**************************************************************************** */
void Time_LinElastic2D::solve_initial_acceleration()
{

}

/**************************************************************************** */
void Time_LinElastic2D::output_problem_size_info() const
{
  int n_u = this->get_fe_space().GetN_DegreesOfFreedom();
  int n_u_active = this->get_fe_space().GetN_ActiveDegrees();
  int n_dof = 2 * n_u; // total number of degrees of freedom

  double h_min, h_max;
  TCollection * coll = this->get_fe_space().GetCollection();
  coll->GetHminHmax(&h_min, &h_max);
  Output::stat("LinElastic2D", "Mesh data and problem size");
  Output::dash("cells              :  ", setw(10), coll->GetN_Cells());
  Output::dash("h (min, max)       :  ", setw(10), h_min, setw(10), " ", h_max);
  Output::dash("dof u_fe_space     :  ", setw(10), n_u );
  Output::dash("dof active         :  ", setw(10), n_u_active);
  Output::dash("dof all            :  ", setw(10), n_dof);
}

