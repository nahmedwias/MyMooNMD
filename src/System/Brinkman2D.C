#include <Brinkman2D.h>
#include <MainUtilities.h> // get velocity and pressure space
#include <LinAlg.h> // DDot
#include <DirectSolver.h>
#include <Assembler4.h>
#include <Assemble2D.h>
#include <BoundaryAssembling2D.h>
#include <Database.h>
#include <ParameterDatabase.h>
#include <memory>
#include <BlockFEMatrix.h>
#include <sstream>
#include <Vertex.h>
#include <AuxParam2D.h>

#include <ParMooN_repository_info.h>

const std::string path = parmoon::source_directory;
const std::string path_to_repo = path + "/data/input_files/";

 void Coefficient_Function(double *in, double *out) 
  {
  // coordinates:  x at in[0], y at in[1]
  // value of conductivity at in[2]
   out[0] = in[2]; 
}

// Create a Brinkman specific database
ParameterDatabase Brinkman2D::get_default_Brinkman2D_parameters()
{
  Output::print<5>("creating a default Brinkman2D parameter database");

  ParameterDatabase brinkman2d_db = ParameterDatabase::parmoon_default_database();
  brinkman2d_db.set_name("Brinkman2D parameter database");

  brinkman2d_db.add("Galerkin_type", "symmetric Galerkin formulation", 
      "This string enables to choose between symmetry and non-symmetry of the standard Galerkin scheme."
      "This might be irrelevant for direct solvers as long as no addiditional terms (stabilization) are apparent.",
      {"symmetric Galerkin formulation", "nonsymmetric Galerkin formulation"});

  brinkman2d_db.add("PkPk_stab", false,
      "Use a Brinkman specific assembling routine corresponding to a "
      "residual-based equal-order stabilization for the Brinkman problem."
      "This only works in two space dimensions and is currently only " 
      "meaningfull for the finite element pair P1/P1."
      "Usually this is used in the main program.",
      {true,false});

  brinkman2d_db.add("EqualOrder_PressureStab_type", "symmetric GLS", 
      "This string sets the type of residual-based momentum stabilization," 
      "usually used for the stabilization of P_k/P_k finite elements, for the" 
      "Brinkman problem in 2D. For stability of the method, the chosen "
      "equal-order stabilization has to fit the other chosen stabilizations and "
      "variants of the standard Galerkin scheme (symmetry or non-symmetry).",
      {"symmetric GLS", "nonsymmetric GLS"});

  brinkman2d_db.add("GradDiv_stab", false,
      "Use an assembling routine for the stabilization of the "
      "divergence constraint.This only works in two space dimensions and"
      " is meaningfull for the finite elemnt space P1/P1 only."
      "Usually this is used in the main program.",
      {true,false});

  brinkman2d_db.add("equal_order_stab_scaling", "by h_T",
      "This string enables to switch between the prefactor h_T^2/(mueff + sigma h_T^2) and h_T^2/(mueff + sigma L_0^2) for some characteristic length of the domain.",
      {"by h_T", "by L_0"});

  brinkman2d_db.add("equal_order_stab_weight_PkPk", (double) 0., "", (double) -1000, (double) 1000 );

  brinkman2d_db.add("refinement_n_initial_steps", (size_t) 2.0 , "", (size_t) 0, (size_t) 10000);

  brinkman2d_db.add("corner_stab_weight", (double)  0.0, "This quantity is the weight of the corner stabilization used for Nitsche corners of the domain", (double) -1000.0 , (double) 1000.0 );

  brinkman2d_db.add("coefficient_function_type", (int) 0 ,
  		"Set the parameter equal to 0 (default value) if the coefficients are constant, resp. explicitly set in the example file. "
		  "If you want to use a coefficient function that is spatially varying, then adjust coeffs via parameters "
		  "in the respective example file and choose: "
  		"coefficient_function_type: 1 -- if the coefficient function will be analytically defined in "
  		"ParMooN_Brinkman2D.C; "
  		"coefficient_function_type: 2 -- if you want to use a coefficient function that is defined by "
  		"a .mesh-file combined with a file describing a corresponding FEFunction2D;"
  		"coefficient_function_type: 3 -- if a mesh and a fitting TFEFunction2D (appropriate for ReadSol()) file are "
  		"provided and only the simulation on the finest grid is of interest. There is a default path to a "
  		"coefficient_function given in the main() (corresponding to example 10, "
  		"boundary_file: ../../../ParMooN/data/mesh/geothermal2d.PRM, "
  		"geo_file: ../../../ParMooN/data/mesh/geothermal2d.mesh, level 0), "
  		"which has to be changed manually to another path if desired.",
			(int) 0 ,(int) 3 );

  brinkman2d_db.add("use_source_sink_function", false, "This activates the use of a function on the rhs "
  		"of the divergence constraint, serving as source or sink. Thus, setting this parameter to true, "
  		"has an influence on the local assembling. Up to now, an approximation of the Dirac distribution "
  		"can be interpolated by setting 'coefficient_function_type' equal to 1. and then be used in the "
  		"assembling. Therefore coeff[9] has to be set equal to parameters[0] in the respective example file "
  		"(see e.g. Geothermal_Energy_Brinkman2D.h).", {true,false});

  brinkman2d_db.add("read_coefficient_function_directory", path_to_repo + "default_written_coefficient_function.Sol" ,
  		"This allows to use coefficient_function_type 2 and 3. "
  		"The File has to fit with the mesh (refinement level). A default file  is contained in input_files/ .");

  brinkman2d_db.add("write_coefficient_function_directory", "." , "");

  brinkman2d_db.add("write_velocity1_directory", "." , "This allows to save the computed Brinkmna2d velocity in "
  		"x-direction in a file for later use ( used in coefficient_function_type 3. ");

  brinkman2d_db.add("write_velocity2_directory", "." , "This allows to save the computed Brinkmna2d velocity in "
  		"y-direction in a file for later use ( used in coefficient_function_type 3. ");

  /* // Possible candidates for own database (maybe also a boundary assembling database, or into the assembling 2d database)
     brinkman2d_db.add("s1", 0.0,
     "Use an assembling routine corresponding to a residual-based "
     "equal-order stabilization for the Brinkman problem."
     "This only works in two space "
     "dimensions and is meaningfull for the finite elemnt spaces P1/P1 and P2/P2 only."
     "Usually this is used in the main program.",
     {-1.0, 0.0, 1.0});

     brinkman2d_db.add("s2", 0.0,
     "Use an assembling routine corresponding to a residual-based "
     "equal-order stabilization for the Brinkman problem."
     "This only works in two space "
     "dimensions and is meaningfull for the finite elemnt spaces P1/P1 and P2/P2 only."
     "Usually this is used in the main program.",
     {-1.0, 0.0, 1.0});

     brinkman2d_db.add("equal_order_stab_weight_PkPk", 0.0,
     "Use an assembling routine corresponding to a residual-based "
     "equal-order stabilization for the Brinkman problem."
     "This only works in two space "
     "dimensions and is meaningfull for the finite elemnt spaces P1/P1 and P2/P2 only."
     "Usually this is used in the main program.",-1000.0,1000.0);
   */
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  brinkman2d_db.merge(out_db, true); // merge the database 'db' with the database out_db

  //Output::print("!!!!!!!!!!Solver type: ", brinkman2d_db["solver_type"]);

  return brinkman2d_db;
}



/** ************************************************************************ */
// Constructor of a Brinkman2D::System_per_grid object
Brinkman2D::System_per_grid::System_per_grid (const Example_Brinkman2D& example,
    TCollection& coll, std::pair<int,int> velocity_pressure_orders,
    Brinkman2D::Matrix type)
  : velocity_space(&coll, "u", "Brinkman velocity", example.get_bc(0),
      velocity_pressure_orders.first, nullptr),
  pressure_space(&coll, "p", "Brinkman pressure", example.get_bc(2),
      velocity_pressure_orders.second, nullptr)
{
  // Note that at the moment only Matrix::Type14 is allowed
  // (see the constructor the of a Brinkman2D object below)
  switch (type)
  {
    case Brinkman2D::Matrix::Type1:
      matrix = BlockFEMatrix::NSE2D_Type1(velocity_space, pressure_space);
      break;
    case Brinkman2D::Matrix::Type2:
      matrix = BlockFEMatrix::NSE2D_Type2(velocity_space, pressure_space);
      break;
    case Brinkman2D::Matrix::Type3:
      matrix = BlockFEMatrix::NSE2D_Type3(velocity_space, pressure_space);
      break;
    case Brinkman2D::Matrix::Type4:
      matrix = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);
      break;
    case Brinkman2D::Matrix::Type14:
      matrix = BlockFEMatrix::NSE2D_Type14(velocity_space, pressure_space);
      break;
    default:
      ErrThrow("Unknown Brinkman type given to constructor of Brinkman2D::System_per_grid.");
  }

  rhs = BlockVector(matrix, true);
  solution = BlockVector(matrix, false);

  u = TFEVectFunct2D(&velocity_space, "u", "u", solution.block(0),
      solution.length(0), 2);
  p = TFEFunction2D(&pressure_space, "p", "p", solution.block(2),
      solution.length(2));
}

/** ************************************************************************ */
Brinkman2D::Brinkman2D(const TDomain& domain, const ParameterDatabase& param_db,
    int reference_id)
: Brinkman2D(domain, param_db, Example_Brinkman2D(param_db),reference_id)
{
  // note that the way we construct the example above will produce a memory
  // leak, but that class is small.
  // FIXME Find a workaround - we do not like memory leaks at all,
  // because they pollute our valgrind tests.
}

/** ************************************************************************ */
Brinkman2D::Brinkman2D(const TDomain & domain, const ParameterDatabase& param_db,
    const Example_Brinkman2D & e,
    unsigned int reference_id)
: brinkman2d_db(Brinkman2D::get_default_Brinkman2D_parameters()), solver(param_db),
  outputWriter(param_db),
  systems(), example(e), defect(), oldResiduals(), 
  initial_residual(1e10), errors()
{
  this->brinkman2d_db.merge(param_db,true);   
  // brinkman2d_db.info(true); 

  // TODO LB 13.11.17    check_and_set_assemble_type(brinkman2d_db["EqualOrder_PressureStab_type"], brinkman2d_db["brinkman_assemble_option"]);

  // set the velocity and preesure spaces
  std::pair <int,int> velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE,
      TDatabase::ParamDB->PRESSURE_SPACE);

  // this function outputs the velocity and pressure order
  this->get_velocity_pressure_orders(velocity_pressure_orders);

  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);

  // Here the constructor of Brinkman2D::System_per_grid is hidden in the so-called double-ended queue 'systems'. This reduces e.g. the amount of copies.
  // We always use Matrix Type 14 such that we can handle pressure stabilizations (Matrix C).
  // TODO LB: implement case distinction for different matrix types according to input parameters.
  //if (PkPk_stab)
  this->systems.emplace_back(example, *coll, velocity_pressure_orders, Brinkman2D::Matrix::Type14);
  //else 
  //this->systems.emplace_back(example, *coll, velocity_pressure_orders, Brinkman2D::Matrix::Type4);

  // the defect has the same structure as the rhs (and as the solution)
  this->defect.copy_structure(this->systems.front().rhs);

  outputWriter.add_fe_vector_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());

  // print out some information
  int n_u = this->get_velocity_space().GetN_DegreesOfFreedom();
  int n_u_active = this->get_velocity_space().GetN_ActiveDegrees();
  int n_p = this->get_pressure_space().GetN_DegreesOfFreedom();
  int n_dof = 2 * n_u + n_p; // total number of degrees of freedom

  double h_min, h_max;
  coll->GetHminHmax(&h_min, &h_max);
  Output::print("------------------------------------------------------");
  Output::print<1>("N_Cells            : ", setw(10), coll->GetN_Cells());
  Output::print<1>("h (min,max)        : ", setw(10), h_min, " ", setw(12), h_max);
  Output::print<1>("dof velocity       : ", setw(10), 2* n_u);
  Output::print<1>("dof velocity active: ", setw(10), 2* n_u_active);
  Output::print<1>("dof pressure       : ", setw(10), n_p);
  Output::print<1>("dof all            : ", setw(10), n_dof);
  Output::print("------------------------------------------------------");
}

/** ************************************************************************ */
Brinkman2D::~Brinkman2D()
{
}

/** ************************************************************************ */
void Brinkman2D::get_velocity_pressure_orders(std::pair <int,int>
    &velocity_pressure_orders)
{
  Output::print<1>("velocity space", setw(10), TDatabase::ParamDB->VELOCITY_SPACE);
  Output::print<1>("pressure space", setw(10), TDatabase::ParamDB->PRESSURE_SPACE);
}

/** ************************************************************************ */
void Brinkman2D::set_parameters()
{
}

/** ************************************************************************ */
void Brinkman2D::check_input_parameters()
{
  // if (brinkman2d_db["EqualOrder_PressureStab_type"].is("symmetric GLS") && brinkman2d_db["Galerkin_type"].is("nonsymmetric Galerkin formulation") || brinkman2d_db["EqualOrder_PressureStab_type"].is("nonsymmetric GLS") && brinkman2d_db["Galerkin_type"].is("symmetric Galerkin formulation")  )
  // { 
  //   Output::print("Warning, the stabilization type does not fit perfectly to the sign of the divergence constraint. This might cause instability.");
  //   //getch();
  //   //system("pause");
  // }
}

/** ************************************************************************ */
void Brinkman2D::assemble(size_t level, TFEFunction2D* coefficient_function)

{
  //Valgrind test start
  //std::vector<int> testvector(5,1);
  //Output::print("testvector.size()", testvector.size());
  //Output::print("testvector[4]", testvector[4]);
  //Output::print("testvector[7]", testvector[7]);
  //Output::print("testvector.at(4)", testvector.at(4));
  //Output::print("testvector.at(7)", testvector.at(7));
  //Valgrind test end

const TFESpace2D *coefficient_function_space;

  for(System_per_grid& s : this->systems)
  {
    s.rhs.reset(); //right hand side reset (TODO: is that necessary?)

    const TFESpace2D *v_space = &s.velocity_space;
    const TFESpace2D *p_space = &s.pressure_space;
    std::array<BoundValueFunct2D*, 3> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    non_const_bound_values[2] = example.get_bd()[2];

    //same for all: the local assembling object
    TFEFunction2D *fe_functions[4] =
    { s.u.GetComponent(0), s.u.GetComponent(1), &s.p, NULL };

    if (coefficient_function)
    {
      coefficient_function_space = coefficient_function->GetFESpace2D();
      fe_functions[3] = coefficient_function;
    }


    if( brinkman2d_db["use_source_sink_function"].is(true))
    {
    	if ( brinkman2d_db["coefficient_function_type"].is(0) )
    	{
    		Error("The chosen parameters 'use_source_sink_function' and 'coefficient_function_type' do not harmonize. "
    				"The source/sink function would not be used in brinkman2d::assemble(). Change the parameter"
    				" 'coefficient_function_type' accordingly!");
    	}
    	else
    	{
    		TDatabase::ParamDB->SOURCE_SINK_FUNCTION = "true";
    	}
    }

    if ( brinkman2d_db["PkPk_stab"].is(false) && brinkman2d_db["Galerkin_type"].is("nonsymmetric Galerkin formulation") )
    {
      TDatabase::ParamDB->SIGN_MATRIX_BI = 1;
    }
    else if ( brinkman2d_db["PkPk_stab"].is(false) && brinkman2d_db["Galerkin_type"].is("symmetric Galerkin formulation") )
    {
      TDatabase::ParamDB->SIGN_MATRIX_BI = -1;
    }
    else if ( brinkman2d_db["PkPk_stab"].is(true) )
    {
      if ( brinkman2d_db["Galerkin_type"].is("nonsymmetric Galerkin formulation") && brinkman2d_db["EqualOrder_PressureStab_type"].is("nonsymmetric GLS") )
      {
        TDatabase::ParamDB->SIGN_MATRIX_BI = 1;
      }
      else if ( brinkman2d_db["PkPk_stab"].is(true)  && brinkman2d_db["Galerkin_type"].is("symmetric Galerkin formulation")  && brinkman2d_db["EqualOrder_PressureStab_type"].is("symmetric GLS") )
      {
        TDatabase::ParamDB->SIGN_MATRIX_BI = -1;
      }
      else if (brinkman2d_db["EqualOrder_PressureStab_type"].is("symmetric GLS") && brinkman2d_db["Galerkin_type"].is("nonsymmetric Galerkin formulation"))
      {
        Output::print("WARNING, the stabilization type does not fit perfectly to the sign of the divergence constraint. This might cause instability. Therefore, both were set to be symmetric. Are you sure you want to continue?");
        TDatabase::ParamDB->SIGN_MATRIX_BI = -1;
        //double tmp;
        //std::cin>>tmp;
      }
      else if (brinkman2d_db["EqualOrder_PressureStab_type"].is("nonsymmetric GLS") && brinkman2d_db["Galerkin_type"].is("symmetric Galerkin formulation"))
      {
        Output::print("WARNING, the stabilization type does not fit perfectly to the sign of the divergence constraint. This might cause instability. Therefore, both were set to be nonsymmetric. Are you sure you want to continue?");
        TDatabase::ParamDB->SIGN_MATRIX_BI = 1;
      }
      else
      { 
        Output::print("WARNING: Please specify the discrete formulation you wish to use via the" 
            " parameters EqualOrder_PressureStab_type and Galerkin_type. ");
      }
    } 

    LocalAssembling_type type;
    // use a list of LocalAssembling2D objects
    std::vector< std::shared_ptr <LocalAssembling2D >> la_list;

    type = LocalAssembling_type::Brinkman2D_Galerkin1;
    std::shared_ptr <LocalAssembling2D> la(
      new LocalAssembling2D(this->brinkman2d_db, type,fe_functions,
                            this->example.get_coeffs()));

       if (coefficient_function)
    {
      // modify la such that it includes the TFEFunction2D coefficient_function
      la->setBeginParameter({0});
      la->setN_ParamFct(1);
      la->setParameterFct({Coefficient_Function});
      la->setN_FeValues(1);
      la->setFeValueFctIndex({3});
      la->setFeValueMultiIndex({D00});
    }

    la_list.push_back(la);
    Output::print<>("The ", brinkman2d_db["Galerkin_type"].value_as_string(), " has been used.");

    if (brinkman2d_db["PkPk_stab"].is(true) && TDatabase::ParamDB->VELOCITY_SPACE >= 1)
    {
      if (brinkman2d_db["equal_order_stab_scaling"].is("by h_T"))
      {
        TDatabase::ParamDB->l_T = 1;
      }
      else if (brinkman2d_db["equal_order_stab_scaling"].is("by L_0"))
      {
        TDatabase::ParamDB->l_T = -1;
      }
      type = LocalAssembling_type::Brinkman2D_Galerkin1ResidualStabP1;
      std::shared_ptr <LocalAssembling2D> la_P1P1stab(
        new LocalAssembling2D(this->brinkman2d_db, type, fe_functions,
                              this->example.get_coeffs()));

if (coefficient_function)
    {
      // modify la such that it includes the TFEFunction2D coefficient_function
      la_P1P1stab->setBeginParameter({0});
      la_P1P1stab->setN_ParamFct(1);
      la_P1P1stab->setParameterFct({Coefficient_Function});
      la_P1P1stab->setN_FeValues(1);
      la_P1P1stab->setFeValueFctIndex({3});
      la_P1P1stab->setFeValueMultiIndex({D00});
    }
     
      la_list.push_back(la_P1P1stab);
      Output::print<>("The ", brinkman2d_db["EqualOrder_PressureStab_type"].value_as_string(), " stabilization for P1/P1 has been applied.");
    }
    if(brinkman2d_db["PkPk_stab"].is(true) && TDatabase::ParamDB->VELOCITY_SPACE > 1)
    {
      if (brinkman2d_db["equal_order_stab_scaling"].is("by h_T"))
      {
        TDatabase::ParamDB->l_T = 1;
      }
      else if (brinkman2d_db["equal_order_stab_scaling"].is("by L_0"))
      {
        TDatabase::ParamDB->l_T = -1;
      }
      type = LocalAssembling_type::Brinkman2D_Galerkin1ResidualStabP2;
      std::shared_ptr <LocalAssembling2D> la_P2P2stab(
        new LocalAssembling2D(this->brinkman2d_db, type, fe_functions,
                              this->example.get_coeffs()));
      
    if (coefficient_function)
    {
      // modify la such that it includes the TFEFunction2D coefficient_function
      la_P2P2stab->setBeginParameter({0});
      la_P2P2stab->setN_ParamFct(1);
      la_P2P2stab->setParameterFct({Coefficient_Function});
      la_P2P2stab->setN_FeValues(1);
      la_P2P2stab->setFeValueFctIndex({3});
      la_P2P2stab->setFeValueMultiIndex({D00});
    }

la_list.push_back(la_P2P2stab);
      Output::print<>("The ", brinkman2d_db["EqualOrder_PressureStab_type"].value_as_string(), " stabilization for Pk/Pk, k>=2,  has been applied.");

    }

    if (brinkman2d_db["GradDiv_stab"].is(true))
    {
      Output::print("Grad-Div stabilization has been applied.");
      type = LocalAssembling_type::Brinkman2D_GradDivStabilization;
      std::shared_ptr <LocalAssembling2D> la_graddiv(
        new LocalAssembling2D(this->brinkman2d_db, type, fe_functions,
                              this->example.get_coeffs()));
   
        if (coefficient_function)
    {
      // modify la such that it includes the TFEFunction2D coefficient_function
      la_graddiv->setBeginParameter({0});
      la_graddiv->setN_ParamFct(1);
      la_graddiv->setParameterFct({Coefficient_Function});
      la_graddiv->setN_FeValues(1);
      la_graddiv->setFeValueFctIndex({3});
      la_graddiv->setFeValueMultiIndex({D00});
    }

   la_list.push_back(la_graddiv);
    }

    /* //--------------------------------------------------------------------------------------------------
    //Old Assemble routine
    declare the variables which Assemble2D needs and each brinkman type has to fill
    size_t N_FESpaces = 2;

    const TFESpace2D *fespmat[2] = {v_space, p_space};
    TSquareMatrix2D *sq_matrices[5]{nullptr}; // it's five pointers maximum (Type14)

    TMatrix2D *rect_matrices[4]{nullptr}; // it's four pointers maximum (Types 2, 4, 14)

    BoundCondFunct2D * boundary_conditions[3] = {
    v_space->GetBoundCondition(), v_space->GetBoundCondition(),
    p_space->GetBoundCondition() };

    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely();
    // Note: We use only Type 14 for Brinkman (for now)
    // call the assemble method with the information that has been patched together
    // //old Assemble2D.C function
    size_t n_sq_mat = 5;
    sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
    sq_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
    sq_matrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
    sq_matrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
    sq_matrices[4] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());

    size_t n_rect_mat = 4;
    rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); // first the lying B blocks
    rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
    rect_matrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); // than the standing B blocks
    rect_matrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());

    double *RHSs[3] = {s.rhs.block(0), s.rhs.block(1), nullptr}; // third place gets only filled
    RHSs[2] = s.rhs.block(2); // NSE type 14 includes pressure rhs
    const TFESpace2D *fesprhs[3] = {v_space, v_space, nullptr};  // if NSE type is 4 or 14
    fesprhs[2]  = p_space;
    size_t N_Rhs = 3; // is 3 if NSE type is 4 or 14 (else it is 2)

    Assemble2D(N_FESpaces, fespmat, n_sq_mat, sq_matrices,
    n_rect_mat, rect_matrices, N_Rhs, RHSs, fesprhs,
    boundary_conditions, non_const_bound_values.data(), la);

    //--------------------------------------------------------------------------------------------------
     */

    // Brinkmann-specific choices
    std::vector<const TFESpace2D*> spaces_for_matrix;
    spaces_for_matrix.resize(2);

    if (coefficient_function)
    { spaces_for_matrix.resize(3);
      spaces_for_matrix[2] = coefficient_function_space;
    }

    spaces_for_matrix[0] = v_space;
    spaces_for_matrix[1] = p_space;

    std::vector<const TFESpace2D*> spaces_for_rhs;
    spaces_for_rhs.resize(3);
    spaces_for_rhs[0] = v_space;
    spaces_for_rhs[1] = v_space;
    spaces_for_rhs[2] = p_space;

    Assembler4 Ass;
    Ass.Assemble2D(s.matrix,s.rhs,
        spaces_for_matrix,spaces_for_rhs,
        example,la_list);
    //===========================================================================//
    // Weakly Imposing Essential Boundary Conditions - Boundary Integrals


    for (int k = 0; k < TDatabase::ParamDB->n_neumann_boundary; k++)
    {
      BoundaryAssembling2D::rhs_g_v_n(s.rhs, v_space,
          nullptr,                                                  // g = 1
          TDatabase::ParamDB->neumann_boundary_id[k],            // boundary component
          -1.*TDatabase::ParamDB->neumann_boundary_value[k]);       // mult
    }

    for (int k = 0; k < TDatabase::ParamDB->n_g_v_boundary; k++)
    {
      BoundaryAssembling2D::rhs_uD_v(s.rhs, v_space,
          this->example.get_bd(0),                                 // access to U1BoundValue in the current example
          this->example.get_bd(1),                                 // access to U2BoundValue in the current example
          TDatabase::ParamDB->g_v_boundary_id[k],                  // boundary component
          TDatabase::ParamDB->g_v_boundary_value[k],               // mult
          false);                                                  // rescale local integral by edge values
    }

    // test: (in the true case, this function should be assembled on "DIRICHLET" boundaries
    // BoundaryAssembling2D::matrix_v_n_v_n(s.matrix,v_space,1,1.);

    for (int k = 0; k < TDatabase::ParamDB->n_unvn_boundary; k++)
    {
      BoundaryAssembling2D::matrix_u_n_v_n(s.matrix, v_space,
          TDatabase::ParamDB->unvn_boundary_id[k],          // boundary component
          TDatabase::ParamDB->unvn_boundary_value[k],      // mult
          false);
    }

    for (int k = 0; k < TDatabase::ParamDB->n_gradunv_boundary; k++)
    {
      BoundaryAssembling2D::matrix_gradu_n_v(s.matrix, v_space,
          TDatabase::ParamDB->gradunv_boundary_id[k],     // boundary component
          TDatabase::ParamDB->gradunv_boundary_value[k]); // mult
    }

    for (int k = 0; k < TDatabase::ParamDB->n_u_v_boundary; k++)
    {
      BoundaryAssembling2D::matrix_u_v(s.matrix, v_space,
          TDatabase::ParamDB->u_v_boundary_id[k],               // boundary component
          TDatabase::ParamDB->u_v_boundary_value[k],            // mult
          false);                                               // rescale local integral by edge values
    }

    for (int k = 0; k < TDatabase::ParamDB->n_p_v_n_boundary; k++)
    {
      BoundaryAssembling2D::matrix_p_v_n(s.matrix, v_space, p_space,
          TDatabase::ParamDB->p_v_n_boundary_id[k],           // boundary component
          TDatabase::ParamDB->p_v_n_boundary_value[k]);       // mult
    }


    //===========================================================================//
    // Nitsche Combination - Weak Essential Boundary Conditions
    for (int k = 0; k < TDatabase::ParamDB->n_nitsche_boundary; k++)
    {
      if (Output::getVerbosity() == 5)
      {   std::string Nitsche_variant;
        std::string Penalty;
        if (TDatabase::ParamDB->s1 == 1 && TDatabase::ParamDB->s2 == 1 )
          Nitsche_variant = "Symmetric";
        else if (TDatabase::ParamDB->s1 == -1 && TDatabase::ParamDB->s2 == -1 )
          Nitsche_variant = "Non-symmetric";
        else if (TDatabase::ParamDB->s1 != TDatabase::ParamDB->s2)
          Nitsche_variant = "Unbalanced";
        else if (TDatabase::ParamDB->nitsche_penalty[k] == 0)
          Penalty = "Penalty-free, ";
        Output::print<5>(Penalty, Nitsche_variant ," Nitsche approach on boundary with ID ", TDatabase::ParamDB->nitsche_boundary_id[k],
            " with Nitsche parameter: ", TDatabase::ParamDB->nitsche_penalty[k], " .");
      }


      BoundaryAssembling2D::matrix_gradu_n_v(s.matrix, v_space,
          TDatabase::ParamDB->nitsche_boundary_id[k],           // boundary component
          -1. * (double)brinkman2d_db["effective_viscosity"]); // mult  ; if all is in t^2, then mult = -1.*t*t is appropriate with t^2 matches brinkman2d_db["effective_viscosity"]

      BoundaryAssembling2D::matrix_gradv_n_u(s.matrix, v_space,
          TDatabase::ParamDB->nitsche_boundary_id[k],           // boundary component
          -1. * TDatabase::ParamDB->s1 * (double) brinkman2d_db["effective_viscosity"]); // mult ; if all is in t^2, then mult =-1.*TDatabase::ParamDB->s1*t*t is appropriate

      BoundaryAssembling2D::rhs_gradv_n_uD(s.rhs, v_space,
          this->example.get_bd(0),                                 // access to U1BoundValue in the example,
          this->example.get_bd(1),                                 // access to U2BoundValue in the example,
          TDatabase::ParamDB->nitsche_boundary_id[k],    // boundary component
          -1. * TDatabase::ParamDB->s1 * (double) brinkman2d_db["effective_viscosity"]); // mult ; if all is in t^2, then mult =-1.*TDatabase::ParamDB->s1*t*t is appropriate

      BoundaryAssembling2D::matrix_p_v_n(s.matrix, v_space, p_space,
          TDatabase::ParamDB->nitsche_boundary_id[k],         // boundary component
          1.);                                                // mult

      BoundaryAssembling2D::matrix_q_u_n(s.matrix, v_space, p_space,
          TDatabase::ParamDB->nitsche_boundary_id[k],         // boundary component
          1. * TDatabase::ParamDB->s2);

      BoundaryAssembling2D::rhs_q_uD_n(s.rhs, v_space, p_space,
          this->example.get_bd(0),                                 // access to U1BoundValue in the example,
          this->example.get_bd(1),                                 // access to U2BoundValue in the example,
          TDatabase::ParamDB->nitsche_boundary_id[k],         // boundary component
          1. * TDatabase::ParamDB->s2);

      BoundaryAssembling2D::matrix_u_v(s.matrix, v_space,
          TDatabase::ParamDB->nitsche_boundary_id[k],           // boundary component
          //t*t*TDatabase::ParamDB->nitsche_penalty[k],   //t*TDatabase::ParamDB->nitsche_penalty[k],             // mult
          TDatabase::ParamDB->nitsche_penalty[k] * (double) brinkman2d_db["effective_viscosity"], //mueff*TDatabase::ParamDB->nitsche_penalty[k],             // mult
          true);                                                // rescale local integral by edge values

      BoundaryAssembling2D::rhs_uD_v(s.rhs, v_space,
          this->example.get_bd(0),                                 // access to U1BoundValue in the example,
          this->example.get_bd(1),                                 // access to U2BoundValue in the example,
          TDatabase::ParamDB->nitsche_boundary_id[k],              // boundary component
          //t*t*TDatabase::ParamDB->nitsche_penalty[k],                // mult
          TDatabase::ParamDB->nitsche_penalty[k]* (double) brinkman2d_db["effective_viscosity"], //mult
          true);            // rescale local integral by edge values 

      // LB TEST 19.01.18
      BoundaryAssembling2D::matrix_u_n_v_n(s.matrix, v_space, 
          TDatabase::ParamDB->nitsche_boundary_id[k],           // boundary component↲
          TDatabase::ParamDB->nitsche_penalty[k]* ((double) brinkman2d_db["viscosity"]/(double) brinkman2d_db["permeability"]) * TDatabase::ParamDB->L_0 * TDatabase::ParamDB->L_0, true);
      //LB TEST 19.01.18^
    }
    //New LB 15.01.18
    BoundaryAssembling2D::matrix_cornerjump_u_n_cornerjump_v_n(s.matrix, v_space, // TDatabase::ParamDB->nitsche_boundary_id[k], 
        1, // nBoundaryParts
        (double) brinkman2d_db["corner_stab_weight"] * ((double) brinkman2d_db["effective_viscosity"] + ((double) brinkman2d_db["viscosity"]/(double) brinkman2d_db["permeability"]) * ( TDatabase::ParamDB->L_0 * TDatabase::ParamDB->L_0 )));

    //NewLB 15.01.18^

    //--------------------------------------------------------------------------------------------------

    // copy Dirichlet values from right hand side into solution
    s.solution.copy_nonactive(s.rhs);

    // create an output file containing the whole FE matrix. This can be read into Matlab using the Matlab function mmread.m
    //s.matrix.get_combined_matrix()->write("Brinkman2D_Matrix_mmread_output");
    // Output::print("Creating output file with FE Matrix");


    // create an output file containing the whole FE matrix. This can be read into Matlab using the Matlab function mmread.m
    /*    std::stringstream matrix_name;
      matrix_name << "Brinkman2D_Matrix_mmread_output_NSGLS"  << TDatabase::ParamDB->equal_order_stab_weight_PkPk << "_GradDivStab"<< TDatabase::ParamDB->grad_div_stab_weight << "_mueff" << std::scientific << setprecision(2) << brinkman2d_db["effective_viscosity"] << "_K" << setprecision(2) << brinkman2d_db["permeability"] << "_level" << level;
       s.matrix.get_combined_matrix()->write(matrix_name.str());
     */  Output::print("Creating output file with FE Matrix");

    // TODO Maybe we have to explicitely set non-actives in non-diagonal blocks
    // to zero here, that was done in former code, but maybe we can move it to the solver part


    // tidy up; This is necessary since GetComponent() is used for the first two entries and this function creates new memory ('new') which has to be deleted manually.
    delete fe_functions[0];
    delete fe_functions[1];
  }
}

/** ************************************************************************ */
bool Brinkman2D::stopIt(unsigned int iteration_counter)
{
  // compute the residuals with the current matrix and solution
  this->computeNormsOfResiduals();
  // the current norm of the residual
  const double normOfResidual = this->getFullResidual();
  // store initial residual, so later we can print the overall reduction
  if(iteration_counter == 0)
  {
    initial_residual = normOfResidual;
  }
  // the residual from 10 iterations ago
  const double oldNormOfResidual = this->oldResiduals.front().fullResidual;

  const unsigned int Max_It = brinkman2d_db["nonlinloop_maxit"];
  const double convergence_speed = brinkman2d_db["nonlinloop_slowfactor"];
  bool slow_conv = false;

  if( normOfResidual >= convergence_speed*oldNormOfResidual )
  {
    slow_conv = true;
  }

  double limit = brinkman2d_db["nonlinloop_epsilon"];
  if ( brinkman2d_db["nonlinloop_scale_epsilon_with_size"] )
  {
    limit *= sqrt(this->get_size());
    Output::print<1>("stopping tolerance for nonlinear iteration ", limit);
  }

  // check if the iteration has converged, or reached the maximum number of
  // iterations or if convergence is too slow. Then return true otherwise false
  if( ( normOfResidual<=limit) || (iteration_counter==Max_It ) || ( slow_conv ) )
  {
    if( slow_conv )
      Output::print<1>(" SLOW !!! ", normOfResidual/oldNormOfResidual);
    // stop iteration
    Output::print<1>(" ITE : ", setw(4), iteration_counter, setprecision(8),
        " RES : ", normOfResidual, " Reduction : ",
        normOfResidual/initial_residual);
    return true;
  }
  else
    return false;
}


/** ************************************************************************ */
void Brinkman2D::computeNormsOfResiduals()
{
  System_per_grid& s = this->systems.front();
  unsigned int n_u_dof = s.solution.length(0);
  unsigned int n_p_dof = s.solution.length(2);

  // copy rhs to defect
  this->defect = s.rhs;
  s.matrix.apply_scaled_add(s.solution, defect,-1.);

  if( s.matrix.pressure_projection_enabled() )
  {
    IntoL20FEFunction(&defect[2*n_u_dof], n_p_dof, &this->get_pressure_space(),
        TDatabase::ParamDB->VELOCITY_SPACE,
        TDatabase::ParamDB->PRESSURE_SPACE);
  }

  // square of the norms of the residual components
  double impuls_Residual = Ddot(2*n_u_dof, &this->defect[0],&this->defect[0]);
  double mass_residual = Ddot(n_p_dof, &this->defect[2*n_u_dof],
      &this->defect[2*n_u_dof]);

  Residuals currentResiduals(impuls_Residual, mass_residual);
  oldResiduals.add(currentResiduals);
}

/** ************************************************************************ */
void Brinkman2D::solve()
{
  System_per_grid& s = this->systems.front();


  /* if (brinkman2d_db["preconditioner"].is("augmented_Lagrangian_based"))
     {
     this->solver.solve_augmented(s.matrix, s.rhs, s.solution);
     }
     else
     {*/
  this->solver.solve(s.matrix, s.rhs, s.solution);
  //}

  //s.rhs.write("righthandside_periodic_test_after");
  //s.matrix.get_combined_matrix()->write("matrix_periodic_test_after");


  // /// @todo consider storing an object of DirectSolver in this class
  //DirectSolver direct_solver(s.matrix,
  //		       DirectSolver::DirectSolverTypes::umfpack);
  // direct_solver.solve(s.rhs, s.solution);
  // 
  // //Output::print("coefficients of solution for basis", s.solution);

  if(s.matrix.pressure_projection_enabled())
    s.p.project_into_L20();

  bool periodic_boundary = true;
if (periodic_boundary)
	 s.p.project_into_L20();
}

/** ************************************************************************ */
void Brinkman2D::output(int level, int i)
{
  bool no_output = !brinkman2d_db["output_write_vtk"] && !brinkman2d_db["output_compute_errors"];
  if( no_output )
    return;

  System_per_grid& s = this->systems.front();
  this->u1 = s.u.GetComponent(0);
  this->u2 = s.u.GetComponent(1);

  /* For Thermohydraulic_Brinkman2D */
  u1->WriteSol(brinkman2d_db["write_velocity1_directory"], "Brinkman_ux");
  u2->WriteSol(brinkman2d_db["write_velocity2_directory"], "Brinkman_uy");


  //----------------------------------------------------------------------------------
  // The following output is made for the geometry channel.mesh or channel_simple.mesh
  if (Output::getVerbosity() == 5)
  {Output::print("The values the solution u1, u2 takes at the line y=1 are saved in u_values_atline.txt.");
    std::ostringstream oss;
        oss << "u_values_LineCut_2Neumann2NSPFNitsche_ExponentialFlow_Perm_" << (double) brinkman2d_db["permeability"] << "EffVisc" << (double) brinkman2d_db["effective_viscosity"] << "Visc" << (double) brinkman2d_db["viscosity"] << "_level_" << level << "scaling_by_L0.txt";

    std::string var = oss.str();
    std::ofstream velfile(var);

    ///   std::ofstream velfile("u_values_atline.txt");
    double values_u1[3];
    double values_u2[3];
    for ( int k=0; k<(4*20)+1; k++ )
    {
      //double x = k*(0.05/4);
      //double y = 0.25; //1.;
      double y = k*(0.05/4);
      double x = 0.25; //1.;

      u1->FindGradient(x, y, values_u1);
      u2->FindGradient(x, y, values_u2);
      velfile << x << " " << y << " " << values_u1[0] << " " << values_u2[0] << endl;
    }
    velfile.close();
  }
  //----------------------------------------------------------------------------------

  //u1->Interpolate(example.exact_solution.GetComponent(0));

  // print the value of the largest and smallest entry in the finite element
  // vector
  if( (size_t) brinkman2d_db["verbosity"]> 1 )
  {
    u1->PrintMinMax();
    u2->PrintMinMax();
    s.p.PrintMinMax();
  }

  outputWriter.add_fe_function(&s.p);
  outputWriter.add_fe_vector_function(&s.u);
  outputWriter.write();

  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(brinkman2d_db["output_compute_errors"])
  {
    std::array<double, 8> errors_temporary;
    errors_temporary.fill(0);
    errors_temporary[1] = s.u.GetL2NormDivergenceError(example.get_exact(0), example.get_exact(1));

    TAuxParam2D aux_error;
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};

    double err[5];
    const TFESpace2D *velocity_space = &this->get_velocity_space();

    // errors in first velocity component
    auto newfunction = [](TBaseCell* cell, int ID){return cell->IsBoundaryCell(ID);};

    u1->GetErrors(example.get_exact(0), 3, AllDerivatives, 2, L2H1Errors,
        nullptr, &aux_error, 1, &velocity_space, err, newfunction);
    // errors in second velocity component
    u2->GetErrors(example.get_exact(1), 3, AllDerivatives, 2, L2H1Errors,
        nullptr, &aux_error, 1, &velocity_space, err + 2, newfunction);


    //cout << "example.get_exact(0)" << example.get_exact(0) << endl;

    errors_temporary.at(0) = sqrt(err[0]*err[0] + err[2]*err[2]);     // err[0]=L2-Error in u1, err[2]=L2-Error in u2
    errors_temporary.at(2) = sqrt(err[1]*err[1] + err[3]*err[3]);     // err[1]=H1-Error in u1, err[3]=H1-Error in u2

    // compute the L2-norm of the velocity error at the Nitsche boundaries
    double boundary_error_l2_u1[1];
    double boundary_error_l2_u2[1];
    u1->GetL2BoundaryError(example.get_bd(0),
        &aux_error, 1, &velocity_space, 
        boundary_error_l2_u1);
    u2->GetL2BoundaryError(example.get_bd(1),
        &aux_error, 1, &velocity_space, 
        boundary_error_l2_u2);

    double boundary_error_l2 = sqrt( boundary_error_l2_u1[0] + boundary_error_l2_u2[0] );
    
    // compute the L2-norm of the normal velocity error at the Nitsche boundaries
    double un_boundary_error = s.u.GetL2NormNormalComponentError(example.get_bd(0), example.get_bd(1));
    
    // compute the L2-norm of the pressure errors
    TAuxParam2D aux;
    const TFESpace2D * pointer_to_p_space = &s.pressure_space;
    s.p.GetErrors(example.get_exact(2), 3, AllDerivatives, 2, L2H1Errors,
        nullptr, &aux, 1, &pointer_to_p_space,
 //example.get_coeffs(), &aux, 1, &pointer_to_p_space,
        errors_temporary.data() + 3);
    
    errors_temporary.at(5) = boundary_error_l2;
    errors_temporary.at(6) = un_boundary_error;


    // Here, the error in the L^2-norm at the boundary is computed.
    //Therefore, a FEMatrix containing the necessary terms and the approximate xolution Blockvector have to be used
    const TFESpace2D * v_space = &s.velocity_space;
    BlockFEMatrix Boundary_matrix_u, Boundary_matrix_u_n;
    Boundary_matrix_u = BlockFEMatrix::NSE2D_Type14(s.velocity_space, s.pressure_space);
    Boundary_matrix_u_n = BlockFEMatrix::NSE2D_Type14(s.velocity_space, s.pressure_space);

    for (int k = 0; k < TDatabase::ParamDB->n_nitsche_boundary; k++)
    {
      BoundaryAssembling2D::matrix_u_v(Boundary_matrix_u, v_space,
          TDatabase::ParamDB->nitsche_boundary_id[k],           // boundary component
          1,                                                    // scale factor
          false);                                               // rescale local integral by edge values

      BoundaryAssembling2D::matrix_u_n_v_n(Boundary_matrix_u_n, v_space, 
          TDatabase::ParamDB->nitsche_boundary_id[k],           // boundary component↲
          1,                                                    // scale factor
          false);                                               // rescale local integral by edge values
    }

    std::vector< double > boundary_errors;
    boundary_errors.resize(2,0);
    // boundary_errors[0]=u^T*A*u;
    BlockVector Au, Au_n;
    Au = BlockVector(Boundary_matrix_u, true);
    Au_n = BlockVector(Boundary_matrix_u_n, true);

    Boundary_matrix_u.print_matrix_info("Boundary_matrix_u Infos");
    Boundary_matrix_u.apply(s.solution, Au);
    boundary_errors[0]=sqrt(dot(s.solution, Au));
    Boundary_matrix_u_n.print_matrix_info("Boundary_matrix_u_n Infos");
    Boundary_matrix_u_n.apply(s.solution, Au_n);
    boundary_errors[1]=sqrt(dot(s.solution, Au_n));

    //Output::print("The value of the boundary norm \sum  || u ||_E is: ", boundary_errors[0]);
    //Output::print("The value of the boundary norm \sum  || u.n ||_E is: is: ", boundary_errors[1]);

    Output::print("------------------------------------------------------");
    Output::print<1>(" L2(u):      ", setprecision(14), errors_temporary[0]);
    Output::print<1>(" L2(div(u)): ", setprecision(14), errors_temporary[1]);
    Output::print<1>(" H1-semi(u): ", setprecision(14), errors_temporary[2]);
    Output::print<1>(" L2(p):      ", setprecision(14), errors_temporary[3]);
    Output::print<1>(" H1-semi(p): ", setprecision(14), errors_temporary[4]);
    //Output::print<1>(" L2(u)_boundary OLD: ", setprecision(14), boundary_errors[0]);
    //New
    Output::print<1>(" L2(u)_boundary: ", setprecision(14), boundary_error_l2);
    //Output::print<1>(" L2(u.n)_boundary OLD: ", setprecision(14), boundary_errors[1]);
    //New
    Output::print<1>(" L2(u.n)_boundary: ", setprecision(14), un_boundary_error);
    Output::print("------------------------------------------------------");

    // copy: This is necessary to save the calculated errors (errors_temporary) in the member variable errors.
    //       The vector errors_temporary is deleted with the end of the current scope and cannot be accessed from outside.
    std::copy(errors_temporary.begin(), errors_temporary.end()-1, this->errors.begin());
 }
  delete u1;
  delete u2;

   // create an output file containing the whole FE matrix. This can be read into Matlab using the Matlab function mmread.m
 std::stringstream matrix_name;
 matrix_name << "Brinkman2D_Matrix_mmread_output_NSGLS"  << TDatabase::ParamDB->equal_order_stab_weight_PkPk << "_GradDivStab"<< TDatabase::ParamDB->grad_div_stab_weight << "_CornerStab"<< brinkman2d_db["corner_stab_weight"] << "_mueff" << std::scientific << setprecision(2) << brinkman2d_db["effective_viscosity"] << "_K" << setprecision(2) << brinkman2d_db["permeability"] <<"_L0_" <<  TDatabase::ParamDB->L_0 << "_level" << level;
 s.matrix.get_combined_matrix()->write(matrix_name.str());


}


/////////////// Routines for periodic boundary conditions /////////////////
/** ************************************************************************ */
void Brinkman2D::findPeriodicDOFs()
{
	// threshold for two doubles to be equal
	const double eps = 1e-8;

	// left and right boundary as periodic boundary
	const double x_left = 0.0;
	const double x_right = 2.0;
	if( ! brinkman2d_db["example"].is(11) )
		ErrThrow("Only the Riverbed example is supported concerning periodic boundaries");

	const TFESpace2D * fespace = &this->get_velocity_space();
	const TCollection* coll = fespace->GetCollection();

	int n_u_active = this->get_velocity_space().GetN_ActiveDegrees();

	//const int n_dofs = fespace->GetN_DegreesOfFreedom();
	int N_Cells = coll->GetN_Cells();
	// first loop over cells
	for(int cell = 0; cell < N_Cells; cell++)
	{
		TBaseCell *cell1 = coll->GetCell(cell);
		// check if face on boundary
		for(int j1 = 0; j1 < cell1->GetN_Edges(); j1++)
		{
			const TJoint *joint1 = cell1->GetJoint(j1);
			//not on boundary
			if(joint1->GetType() != BoundaryEdge)
				continue;

			const int n_Vert = cell1->GetN_Vertices();
			double x11, x12, y11, y12;
			// compute coordinates of vertices
			cell1->GetVertex(j1)->GetCoords(x11, y11);
			cell1->GetVertex((j1 + 1) % n_Vert)->GetCoords(x12, y12);

			//Todo: change implementation such that periodic dofs in one and the same cell are kept as well


			// check if vertex x = x_left
			bool left = false;
			if(fabs(x11 - x_left) > eps || fabs(x12 - x_left) > eps)
			{ // one vertex does not lie on the left boundary
				if(fabs(x11 - x_right) > eps || fabs(x12 - x_right) > eps)
					continue; // one vertex does not lie on the right boundary
			}
			else
				{
				left = true;
				}

			const FE2D FEid1 = fespace->GetFE2D(cell, cell1);
			const TFE2D * FE1 = TFEDatabase2D::GetFE2D(FEid1);

			// global indices of all degrees of freedom in this cell
			const int *globalDOF1 = fespace->GetGlobalDOF(cell);
			// local degrees of freedom which correspond to this edge
			const int* localDOF1 = FE1->GetFEDesc2D()->GetJointDOF(j1);
			const int N_localDOF1 = FE1->GetFEDesc2D()->GetN_JointDOF();

			// midpoint of edge (specific for periodicity in x direction)
			const double y1 = (y11 + y12) / 2;

			// find cell which should be coupled to this cell
			// inner loop over the cells
			for(int i2_cell = cell + 1; i2_cell < N_Cells; i2_cell++)
			{
				TBaseCell *cell2 = coll->GetCell(i2_cell);
				// check if face on boundary
				for(int j2 = 0; j2 < cell2->GetN_Edges(); j2++)
				{
					const TJoint *joint2 = cell2->GetJoint(j2);
					//not on boundary
					if(joint2->GetType() != BoundaryEdge)
						continue;

					double x21, x22, y21, y22;
					cell2->GetVertex(j2)->GetCoords(x21, y21);
					cell2->GetVertex((j2 + 1) % n_Vert)->GetCoords(x22, y22);

					if(fabs(x21 - (left ? x_right : x_left)) > eps  || fabs(x22 - (left ? x_right : x_left)) > eps)
						continue; // one vertex does not lie on the correct boundary
					// the two edges are on the correct boundary parts
					// check if their midpoints have the same y-coordinate
					const double y2 = (y21 + y22) / 2;
					if(fabs(y1 - y2) > eps)
						continue;

					// found two edges will should be identified

					const FE2D FEid2 = fespace->GetFE2D(i2_cell, cell2);
					const TFE2D * FE2 = TFEDatabase2D::GetFE2D(FEid2);

					// global indices of all degrees of freedom in this cell
					const int *globalDOF2 = fespace->GetGlobalDOF(i2_cell);
					// local degrees of freedom which correspond to this edge
					const int* localDOF2 = FE2->GetFEDesc2D()->GetJointDOF(j2);
					const int N_localDOF2 = FE2->GetFEDesc2D()->GetN_JointDOF();

					if(FEid1 != FEid2)
					{
						ErrThrow("Error in making periodic boundary. ",
								"Two different finite elements");
					}
					if(N_localDOF1 != N_localDOF2)
					{
						ErrThrow("Error in making periodic boundary. ",
								"Different numbers of dofs on the periodic boundary");
					}

					/*if(TDatabase::ParamDB->SC_VERBOSE > 2)
						OutPut( " creating a vertical periodic boundary at y=(" << y21 << "," << y22 << ")\n");
								*/

					for(int edge_dof = 0; edge_dof < N_localDOF1; edge_dof++)
					{
						// due to counterclockwise numbering in each cell we have to go
						// through one edge the opposite way:
						const int dof1 = globalDOF1[localDOF1[N_localDOF1 - 1 - edge_dof]];
						const int dof2 = globalDOF2[localDOF2[edge_dof]];



						if(dof1 >= n_u_active && dof2 >= n_u_active)
							continue;
						else if ((dof1 >= n_u_active && dof2 < n_u_active) || (dof1 < n_u_active && dof2 >= n_u_active))
						{
							ErrThrow("Error in findPeriodicDOFs() - active dofs.");
						}

						if (left)
						{
							periodic_dofs[dof1] = dof2;
						// LB DEBUG 03.09.18	periodic_dofs[dof1 + n_dofs] = dof2 + n_dofs; // second component  WOZU? LB
							//OutPut(" dofs " << dof1 << "\t" << dof2 << endl);
							//OutPut("      " << dof1+n_dofs << "\t" << dof2+n_dofs << endl)
						 /* cout<< "x11: "<< x11 << ", y11: "<< y11 << endl;
						  cout<< "x12: "<< x12 << ", y12: "<< y12 << endl;
						  cout<< "x21: "<< x21 << ", y21: "<< y21 << endl;
						  cout<< "x22: "<< x22 << ", y22: "<< y22 << endl;
				      cout<< "dof1: "<< dof1 << ", dof2: "<<  dof2<< endl;
				      cout<< "dof1 + n_dofs: "<< dof1 + n_dofs << ", dof2 + n_dofs: "<<  dof2 + n_dofs<< endl;
					  	*/
						}
			//LB Debug 28.08.18
				else
						{
							periodic_dofs[dof2] = dof1;
							// LB DEBUG 03.09.18		periodic_dofs[dof2 + n_dofs] = dof1 + n_dofs; // second component  WOZU? LB
							//OutPut(" dofs " << dof2 << "\t" << dof1 << endl);
							//OutPut("      " << dof2+n_dofs << "\t" << dof1+n_dofs << endl)
							/*cout<< "ELSE: "<< endl;
							 cout<< "x11: "<< x11 << ", y11: "<< y11 << endl;
													  cout<< "x12: "<< x12 << ", y12: "<< y12 << endl;
													  cout<< "x21: "<< x21 << ", y21: "<< y21 << endl;
													  cout<< "x22: "<< x22 << ", y22: "<< y22 << endl;
						  cout<< "dof1: "<< dof1 << ", dof2: "<<  dof2<< endl;
						  cout<< "dof1 + n_dofs: "<< dof1 + n_dofs << ", dof2 + n_dofs: "<<  dof2 + n_dofs<< endl;
						*/
						}
					}
				}
			}
		}
	}
	OutPut( "There are " << periodic_dofs.size()/2 << " periodic Brinkman2D degrees of freedom\n");
}


/** ************************************************************************ */
void Brinkman2D::makePeriodicBoundary(std::shared_ptr<TMatrix> mat,
		bool stokesMat, bool p )
{
	if(periodic_dofs.empty())
		ErrThrow("called Brinkman2D::makePeriodicBoundary with map ",
				"'periodic_dofs' not yet set");

	int do_once = 0;

	if(!mat)
	{
		auto blocks = this->get_matrix().get_blocks_uniquely();
		makePeriodicBoundary(blocks[0], true, true);  // A11
		do_once +=1 ;
		makePeriodicBoundary(blocks[1], false, true); // A12
		makePeriodicBoundary(blocks[2], false, true); // B1T
		makePeriodicBoundary(blocks[3], false, true); // A21
		makePeriodicBoundary(blocks[4], true, true);  // A22
		makePeriodicBoundary(blocks[5], false, true); // B2T
		// TODO add remaining blocks for Brinkman problem
		//makePeriodicBoundary(blocks[6], true, false); // B1
		//makePeriodicBoundary(blocks[7], true, false); // B2
		//makePeriodicBoundary(blocks[8], false, true); // C
		return;
	}

	int const * const rowPtr = mat->GetRowPtr();
	int const * const colPtr = mat->GetKCol();
	double const * const entries = mat->GetEntries();
	const int n_rows = mat->GetN_Rows();
	// this map will be passed to "mat->changeRows(new_rows)"·
	std::map<int, std::map<int, double> > new_rows;

//	for(int k = 0; k < this->get_rhs().length(0) + 1; k++ )
//	cout<<"this->get_rhs()[k]:  " << this->get_rhs()[k] << endl;
//	cout<<"DONE!!!!!!!!!!" << endl;

	for(std::map<int, int>::iterator ii = periodic_dofs.begin();
			ii != periodic_dofs.end() && ii->first < n_rows; ++ii)
	{
		//OutPut( ii->first << ": " << ii->second << endl);
		// the row with number "ii->second" is added to the row number "ii->first".
		// Then row number "ii->second" is replaced by a row with two entries, 1·
		// on the diagonal and -1 on the "ii->first"-th entry. Also the right hand·
		// side is changed to be zero

		// loop over all entries in row "ii->first" of A-matrices
		for(int i = rowPtr[ii->first]; i < rowPtr[1 + ii->first]; i++)
		{
			if(entries[i] != 0.0 || p)
				(new_rows[ii->first])[colPtr[i]] += entries[i];
			//cout<< "add row " <<  ii->first << endl;

//			cout<<"this->get_rhs()[i]:  " << this->get_rhs()[i] << endl;
		} // up to here this row would simply be copied.
		// loop over all entries in row "ii->second" of A-matrices
		for(int i = rowPtr[ii->second]; i < rowPtr[1 + ii->second]; i++)
		{
			if(entries[i] != 0.0 || p)
				{
				(new_rows[ii->first])[colPtr[i]] += entries[i];
//				cout<<"this->get_rhs()[i]:  " << this->get_rhs()[i] << endl;
				//cout<< "and row " <<  ii->second << endl;
				}
		}

//New LB 28.08.18 start
		//cout<<"ii->first:  " << ii->first << ", ii->second: " << ii->second<< endl;
		//cout<<"OLD: this->get_rhs()[ii->first]:  " << this->get_rhs()[ii->first] << endl;
		if (stokesMat && p)
		{
			if (do_once == 0)
			{
		this->get_rhs()[ii->first] += this->get_rhs()[ii->second];
	 		this->get_rhs()[ii->first + this->get_rhs().length(0)] += this->get_rhs()[ii->second + this->get_rhs().length(0)];
			}
		}
	//	cout<<"JOOOOOOOOOOO" << endl;
	//	cout<<"NEW: this->get_rhs()[ii->first]:  " << this->get_rhs()[ii->first] << endl;
  //New LB 28.08.18 end


		if(stokesMat)
		{
			//LB Debug 30.08.18 new (new_rows[ii->second]); // set entire row to zero
			(new_rows[ii->second])[ii->second] = 1.0; // diagonal
			(new_rows[ii->second])[ii->first] = -1.0; // coupling
			// here two entries for the right hand side are handled.
			// One would be enough, but which one depends on which matrix this is
			// (A11 or A22).


//cout<< "1 -1 pair in row: "<< ii->second << endl;

			this->get_rhs()[ii->second] = 0;
	 		this->get_rhs()[ii->second + this->get_rhs().length(0)] = 0;
		}
		else //LB Debug 30.08.18 old if(!p || n_rows != mat->GetN_Columns())
			{
			(new_rows[ii->second]); // set entire row to zero
			//cout<< "zero row in row "<< ii->second << endl;
			}
		/* //LB Debug 30.08.18 old
		 * else
		{
			(new_rows[ii->second])[ii->second] = 0.0; // diagonal
			(new_rows[ii->second])[ii->first] = 0.0; // coupling
			cout<< "zero entries" << endl;
		}*/
	}
	mat->changeRows(new_rows);
}

/** ************************************************************************ */
void Brinkman2D::checkPeriodicDOFs()
{
	std::map<int, int>::iterator it;
	double const * const vals = this->get_solution().get_entries();
	for(it = periodic_dofs.begin(); it != periodic_dofs.end(); ++it)
	{
		double diff = vals[it->first] - vals[it->second];
		if(fabs(diff) > 1e-6)
			Error("WRONG Brinkman2D here " << it->first << "\t" << it->second << "\t" << setprecision(10) << diff << endl);
	}
}


/////////////// ///////////////  ///////////////  ///////////////

/** ************************************************************************ */
double Brinkman2D::getL2VelocityError() const
{
  return this->errors[0];
}

/** ************************************************************************ */
double Brinkman2D::getL2DivergenceError() const
{
  return this->errors[1];
}

/** ************************************************************************ */
double Brinkman2D::getH1SemiVelocityError() const
{
  return this->errors[2];
}

/** ************************************************************************ */
double Brinkman2D::getL2PressureError() const
{
  return this->errors[3];
}

/** ************************************************************************ */
double Brinkman2D::getL2BoundaryError() const
{
  return this->errors[5];
}
/** ************************************************************************ */
double Brinkman2D::getL2NormNormalComponentError() const
{
  return this->errors[6];
}

/** ************************************************************************ */
double Brinkman2D::getH1SemiPressureError() const
{
  return this->errors[4];
}


/** ************************************************************************ */
std::array< double, int(8) > Brinkman2D::get_errors()
{
  return errors;
}

/** ************************************************************************ */
const Brinkman2D::Residuals& Brinkman2D::getResiduals() const
{
  return this->oldResiduals.back();
}

/** ************************************************************************ */
double Brinkman2D::getImpulsResidual() const
{
  return this->oldResiduals.back().impulsResidual;
}

/** ************************************************************************ */
double Brinkman2D::getMassResidual() const
{
  return this->oldResiduals.back().massResidual;
}

/** ************************************************************************ */
double Brinkman2D::getFullResidual() const
{
  return this->oldResiduals.back().fullResidual;
}

/** ************************************************************************ */
  Brinkman2D::Residuals::Residuals()
: impulsResidual(1e10), massResidual(1e10), fullResidual(1e10)
{}

/** ************************************************************************ */
  Brinkman2D::Residuals::Residuals(double imR, double maR)
: impulsResidual(sqrt(imR)), massResidual(sqrt(maR)),
  fullResidual(sqrt(imR + maR))
{}

/** ************************************************************************ */
std::ostream& operator<<(std::ostream& s, const Brinkman2D::Residuals& n)
{
  s << setw(14) << n.impulsResidual << "\t" << setw(14)
    << n.massResidual << "\t" << setw(14) << n.fullResidual;
  return s;
}


