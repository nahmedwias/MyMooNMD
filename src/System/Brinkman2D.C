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


// Create a Brinkman specific database
ParameterDatabase Brinkman2D::get_default_Brinkman2D_parameters()
{
    Output::print<5>("creating a default Brinkman2D parameter database");

    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.set_name("Brinkman2D parameter database");

    db.add("Galerkin_type", "symmetric Galerkin formulation", 
           "This string enables to choose between symmetry and non-symmetry of the standard Galerkin scheme."
           "This might be irrelevant for direct solvers as long as no addiditional terms (stabilization) are apparent.",
           {"symmetric Galerkin formulation", "nonsymmetric Galerkin formulation"});

    db.add("PkPk_stab", false,
           "Use a Brinkman specific assembling routine corresponding to a "
           "residual-based equal-order stabilization for the Brinkman problem."
           "This only works in two space dimensions and is currently only " 
           "meaningfull for the finite element pair P1/P1."
           "Usually this is used in the main program.",
           {true,false});

     db.add("EqualOrder_PressureStab_type", "symmetric GLS", 
           "This string sets the type of residual-based momentum stabilization," 
           "usually used for the stabilization of P_k/P_k finite elements, for the" 
           "Brinkman problem in 2D. For stability of the method, the chosen "
           "equal-order stabilization has to fit the other chosen stabilizations and "
           "variants of the standard Galerkin scheme (symmetry or non-symmetry).",
           {"symmetric GLS", "nonsymmetric GLS"});

     db.add("GradDiv_stab", false,
            "Use an assembling routine for the stabilization of the "
            "divergence constraint.This only works in two space dimensions and"
            " is meaningfull for the finite elemnt space P1/P1 only."
            "Usually this is used in the main program.",
            {true,false});
    db.add("equal_order_stab_scaling", "by h_T",
           "This string enables to switch between the prefactor h_T^2/(mueff + sigma h_T^2) and h_T^2/(mueff + sigma L_0^2) for some characteristic length of the domain.",
            {"by h_T", "by L_0"});
      db.add("equal_order_stab_weight_PkPk", (double) 0., "", (double) -1000, (double) 1000 );

      db.add("refinement_n_initial_steps", (size_t) 2.0 , "", (size_t) 0, (size_t) 10000);


/* // Possible candidates for own database (maybe also a boundary assembling database, or into the assembling 2d database)
    db.add("s1", 0.0,
           "Use an assembling routine corresponding to a residual-based "
           "equal-order stabilization for the Brinkman problem."
           "This only works in two space "
           "dimensions and is meaningfull for the finite elemnt spaces P1/P1 and P2/P2 only."
           "Usually this is used in the main program.",
           {-1.0, 0.0, 1.0});

    db.add("s2", 0.0,
           "Use an assembling routine corresponding to a residual-based "
           "equal-order stabilization for the Brinkman problem."
           "This only works in two space "
           "dimensions and is meaningfull for the finite elemnt spaces P1/P1 and P2/P2 only."
           "Usually this is used in the main program.",
           {-1.0, 0.0, 1.0});

    db.add("equal_order_stab_weight_PkPk", 0.0,
           "Use an assembling routine corresponding to a residual-based "
           "equal-order stabilization for the Brinkman problem."
           "This only works in two space "
           "dimensions and is meaningfull for the finite elemnt spaces P1/P1 and P2/P2 only."
           "Usually this is used in the main program.",-1000.0,1000.0);
*/
    ParameterDatabase out_db = ParameterDatabase::default_output_database();
    db.merge(out_db, true); // merge the database 'db' with the database out_db
    
    return db;
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
: db(Brinkman2D::get_default_Brinkman2D_parameters()), solver(param_db),
  outputWriter(param_db),
  systems(), example(e), defect(), oldResiduals(), 
  initial_residual(1e10), errors()
{
    this->db.merge(param_db,true);   
   // db.info(true); 

// TODO LB 13.11.17    check_and_set_assemble_type(db["EqualOrder_PressureStab_type"], db["brinkman_assemble_option"]);

   // set the velocity and preesure spaces
   std::pair <int,int> velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE,
                                                TDatabase::ParamDB->PRESSURE_SPACE);

    // this function outputs the velocity and pressure order
    this->get_velocity_pressure_orders(velocity_pressure_orders);
    
    // create the collection of cells from the domain (finest grid)
    TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
    
    // Here the constructor of Brinkman2D::System_per_grid is hidden in the so-called double-ended queue 'systems'. This reduces e.g. the amount of copies.
    // We always use Matrix Type 14 such that we can handle pressure stabilizations (Matrix C).
    // TODO LB: implement case distiction for diffrerent matrix types according to input parameters.
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
 // if (db["EqualOrder_PressureStab_type"].is("symmetric GLS") && db["Galerkin_type"].is("nonsymmetric Galerkin formulation") || db["EqualOrder_PressureStab_type"].is("nonsymmetric GLS") && db["Galerkin_type"].is("symmetric Galerkin formulation")  )
 // { 
 //   Output::print("Warning, the stabilization type does not fit perfectly to the sign of the divergence constraint. This might cause instability.");
 //   //getch();
 //   //system("pause");
 // }
}

/** ************************************************************************ */
void Brinkman2D::assemble()
{
  //Valgrind test start
  //std::vector<int> testvector(5,1);
  //Output::print("testvector.size()", testvector.size());
  //Output::print("testvector[4]", testvector[4]);
  //Output::print("testvector[7]", testvector[7]);
  //Output::print("testvector.at(4)", testvector.at(4));
  //Output::print("testvector.at(7)", testvector.at(7));
  //Valgrind test end

  for(System_per_grid& s : this->systems)
  {
    s.rhs.reset(); //right hand side reset (TODO: is that necessary?)

    const TFESpace2D * v_space = &s.velocity_space;
    const TFESpace2D * p_space = &s.pressure_space;
    std::array<BoundValueFunct2D*, 3> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    non_const_bound_values[2] = example.get_bd()[2];

    //same for all: the local asembling object
    TFEFunction2D *fe_functions[3] =
    { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };

    if ( db["PkPk_stab"].is(false) && db["Galerkin_type"].is("nonsymmetric Galerkin formulation") )
    {
      TDatabase::ParamDB->SIGN_MATRIX_BI = 1;
    }
    else if ( db["PkPk_stab"].is(false) && db["Galerkin_type"].is("symmetric Galerkin formulation") )
    {
      TDatabase::ParamDB->SIGN_MATRIX_BI = -1;
    }
    else if ( db["PkPk_stab"].is(true) )
    {
      if ( db["Galerkin_type"].is("nonsymmetric Galerkin formulation") && db["EqualOrder_PressureStab_type"].is("nonsymmetric GLS") )
      {
        TDatabase::ParamDB->SIGN_MATRIX_BI = 1;
      }
      else if ( db["PkPk_stab"].is(true)  && db["Galerkin_type"].is("symmetric Galerkin formulation")  && db["EqualOrder_PressureStab_type"].is("symmetric GLS") )
      {
        TDatabase::ParamDB->SIGN_MATRIX_BI = -1;
      }
      else if (db["EqualOrder_PressureStab_type"].is("symmetric GLS") && db["Galerkin_type"].is("nonsymmetric Galerkin formulation"))
      {
        Output::print("WARNING, the stabilization type does not fit perfectly to the sign of the divergence constraint. This might cause instability. Therefore, both were set to be symmetric. Are you sure you want to continue?");
         TDatabase::ParamDB->SIGN_MATRIX_BI = -1;
       //double tmp;
        //std::cin>>tmp;
      }
      else if (db["EqualOrder_PressureStab_type"].is("nonsymmetric GLS") && db["Galerkin_type"].is("symmetric Galerkin formulation"))
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

    LocalAssembling2D_type type;
    // use a list of LocalAssembling2D objects
    std::vector< std::shared_ptr <LocalAssembling2D >> la_list;

    type = Brinkman2D_Galerkin1;
    std::shared_ptr <LocalAssembling2D> la(new LocalAssembling2D(type, fe_functions,
          this->example.get_coeffs()));
    la_list.push_back(la);
    Output::print<>("The ", db["Galerkin_type"].value_as_string(), " has been used.");


    if (db["PkPk_stab"].is(true) && TDatabase::ParamDB->VELOCITY_SPACE == 1)
    {
      if (db["equal_order_stab_scaling"].is("by h_T"))
        {
        TDatabase::ParamDB->l_T = 1;
        }
      else if (db["equal_order_stab_scaling"].is("by L_0"))
        {
        TDatabase::ParamDB->l_T = -1;
        }
      type = Brinkman2D_Galerkin1ResidualStabP1;
      std::shared_ptr <LocalAssembling2D> la_P1P1stab(new LocalAssembling2D(type, fe_functions,
            this->example.get_coeffs()));
      la_list.push_back(la_P1P1stab);
      Output::print<>("The ", db["EqualOrder_PressureStab_type"].value_as_string(), " stabilization for P1/P1 has been applied.");
    }
    else if(db["PkPk_stab"].is(true) && TDatabase::ParamDB->VELOCITY_SPACE != 1)
    {
      Output::print<>("Up to now, the GLS stabilization is only implemented correctly for P1/P1, this does not fit the input data.");
    }

    if (db["GradDiv_stab"].is(true))
    {
      Output::print("Grad-Div stabilization has been applied.");
      type = Brinkman2D_GradDivStabilization;
      std::shared_ptr <LocalAssembling2D> la_graddiv(new LocalAssembling2D(type, fe_functions,
            this->example.get_coeffs()));
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

    BoundaryAssembling2D bi;

    for (int k = 0; k < TDatabase::ParamDB->n_neumann_boundary; k++)
    {
      bi.rhs_g_v_n(s.rhs, v_space,
          NULL,                                                  // g = 1
          TDatabase::ParamDB->neumann_boundary_id[k],            // boundary component
          -1.*TDatabase::ParamDB->neumann_boundary_value[k]);       // mult
    }

    for (int k = 0; k < TDatabase::ParamDB->n_g_v_boundary; k++)
    {
      bi.rhs_uD_v(s.rhs, v_space,
          this->example.get_bd(0),                                 // access to U1BoundValue in the current example
          this->example.get_bd(1),                                 // access to U2BoundValue in the current example
          TDatabase::ParamDB->g_v_boundary_id[k],                  // boundary component
          TDatabase::ParamDB->g_v_boundary_value[k],               // mult
          false);                                                  // rescale local integral by edge values
    }

    // test: (in the true case, this function should be assembled on "DIRICHLET" boundaries
    // bi.matrix_v_n_v_n(s.matrix,v_space,1,1.);

    for (int k = 0; k < TDatabase::ParamDB->n_unvn_boundary; k++)
    {
      bi.matrix_u_n_v_n(s.matrix, v_space,
          TDatabase::ParamDB->unvn_boundary_id[k],          // boundary component
          TDatabase::ParamDB->unvn_boundary_value[k]);      // mult
    }

    for (int k = 0; k < TDatabase::ParamDB->n_gradunv_boundary; k++)
    {
      bi.matrix_gradu_n_v(s.matrix, v_space,
          TDatabase::ParamDB->gradunv_boundary_id[k],     // boundary component
          TDatabase::ParamDB->gradunv_boundary_value[k]); // mult
    }

    for (int k = 0; k < TDatabase::ParamDB->n_u_v_boundary; k++)
    {
      bi.matrix_u_v(s.matrix, v_space,
          TDatabase::ParamDB->u_v_boundary_id[k],               // boundary component
          TDatabase::ParamDB->u_v_boundary_value[k],            // mult
          false);                                               // rescale local integral by edge values
    }

    for (int k = 0; k < TDatabase::ParamDB->n_p_v_n_boundary; k++)
    {
      bi.matrix_p_v_n(s.matrix, v_space, p_space,
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

      bi.matrix_gradu_n_v(s.matrix, v_space,
          TDatabase::ParamDB->nitsche_boundary_id[k],           // boundary component
          -1. * TDatabase::ParamDB->EFFECTIVE_VISCOSITY);      // mult  ; if all is in t^2, then mult = -1.*t*t is appropriate with t^2 matches TDatabase::ReadParam->EFFECTIVE_VISCOSITY

      bi.matrix_gradv_n_u(s.matrix, v_space,
          TDatabase::ParamDB->nitsche_boundary_id[k],           // boundary component
          -1. * TDatabase::ParamDB->s1 * TDatabase::ParamDB->EFFECTIVE_VISCOSITY);                   // mult ; if all is in t^2, then mult =-1.*TDatabase::ParamDB->s1*t*t is appropriate

      bi.rhs_gradv_n_uD(s.rhs, v_space,
          this->example.get_bd(0),                                 // access to U1BoundValue in the example,
          this->example.get_bd(1),                                 // access to U2BoundValue in the example,
          TDatabase::ParamDB->nitsche_boundary_id[k],    // boundary component
          -1. * TDatabase::ParamDB->s1 * TDatabase::ParamDB->EFFECTIVE_VISCOSITY);   // mult ; if all is in t^2, then mult =-1.*TDatabase::ParamDB->s1*t*t is appropriate

      bi.matrix_p_v_n(s.matrix, v_space, p_space,
          TDatabase::ParamDB->nitsche_boundary_id[k],         // boundary component
          1.);                                                // mult

      bi.matrix_q_u_n(s.matrix, v_space, p_space,
          TDatabase::ParamDB->nitsche_boundary_id[k],         // boundary component
          1. * TDatabase::ParamDB->s2);

      bi.rhs_q_uD_n(s.rhs, v_space, p_space,
          this->example.get_bd(0),                                 // access to U1BoundValue in the example,
          this->example.get_bd(1),                                 // access to U2BoundValue in the example,
          TDatabase::ParamDB->nitsche_boundary_id[k],         // boundary component
          1. * TDatabase::ParamDB->s2);

      bi.matrix_u_v(s.matrix, v_space,
          TDatabase::ParamDB->nitsche_boundary_id[k],           // boundary component
          //t*t*TDatabase::ParamDB->nitsche_penalty[k],   //t*TDatabase::ParamDB->nitsche_penalty[k],             // mult
          TDatabase::ParamDB->nitsche_penalty[k],   //t*TDatabase::ParamDB->nitsche_penalty[k],             // mult
          true);                                                // rescale local integral by edge values

      bi.rhs_uD_v(s.rhs, v_space,
          this->example.get_bd(0),                                 // access to U1BoundValue in the example,
          this->example.get_bd(1),                                 // access to U2BoundValue in the example,
          TDatabase::ParamDB->nitsche_boundary_id[k],              // boundary component
          //t*t*TDatabase::ParamDB->nitsche_penalty[k],                // mult
          TDatabase::ParamDB->nitsche_penalty[k],                // mult
          true);            // rescale local integral by edge values
    }


    //--------------------------------------------------------------------------------------------------

    // copy Dirichlet values from right hand side into solution
    s.solution.copy_nonactive(s.rhs);

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
    
    const unsigned int Max_It = db["nonlinloop_maxit"];
    const double convergence_speed = db["nonlinloop_slowfactor"];
    bool slow_conv = false;
    
    if( normOfResidual >= convergence_speed*oldNormOfResidual )
    {
        slow_conv = true;
    }
    
    double limit = db["nonlinloop_epsilon"];
    if ( db["nonlinloop_scale_epsilon_with_size"] )
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
   this->solver.solve(s.matrix, s.rhs, s.solution);

   // /// @todo consider storing an object of DirectSolver in this class
   //DirectSolver direct_solver(s.matrix,
	//		       DirectSolver::DirectSolverTypes::umfpack);
   // direct_solver.solve(s.rhs, s.solution);
   // 
   // //Output::print("coefficients of solution for basis", s.solution);
    
    if(s.matrix.pressure_projection_enabled())
      s.p.project_into_L20();
}

/** ************************************************************************ */
void Brinkman2D::output(int i)
{
  bool no_output = !db["output_write_vtk"] && !db["output_compute_errors"];
  if( no_output )
    return;

  System_per_grid& s = this->systems.front();
  TFEFunction2D* u1 = s.u.GetComponent(0);
  TFEFunction2D* u2 = s.u.GetComponent(1);

  //----------------------------------------------------------------------------------
  // The following output is made for the geometry channel.mesh or channel_simple.mesh
  if (Output::getVerbosity() == 5)
  {Output::print("The values the solution u1, u2 takes at the line y=1 are saved in u_values_atline.txt.");
    std::ofstream velfile("u_values_atline.txt");
    double values_u1[3];
    double values_u2[3];
    for ( int k=0; k<(4*20)+1; k++ )
    {
      double x = k*(0.05/4);
      double y=1.;
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
  if( (size_t)db["verbosity"]> 1 )
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
  if(db["output_compute_errors"])
  {
    std::array<double, 6> errors_temporary;
    errors_temporary.fill(0);
    errors_temporary[1] = s.u.GetL2NormDivergence();

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

    errors_temporary.at(0) = sqrt(err[0]*err[0] + err[2]*err[2]);     // err[0]=L2-Error in u1, err[2]=L2-Error in u2
    errors_temporary.at(2) = sqrt(err[1]*err[1] + err[3]*err[3]);     // err[1]=H1-Error in u1, err[3]=H1-Error in u2

    TAuxParam2D aux;
    const TFESpace2D * pointer_to_p_space = &s.pressure_space;
    s.p.GetErrors(example.get_exact(2), 3, AllDerivatives, 2, L2H1Errors,
        example.get_coeffs(), &aux, 1, &pointer_to_p_space,
        errors_temporary.data() + 3);

    /* LB 14.11. OLD  (alternative computation of the pressure errors)
       TAuxParam2D aux_error;
    // errors in pressure
    const TFESpace2D *pressure_space = &this->get_pressure_space();
    s.p.GetErrors(example.get_exact(2), 3, AllDerivatives, 2, L2H1Errors,
    nullptr, &aux_error, 1, &pressure_space, err, newfunction);

    errors_temporary.at(2) = err[0];    // L2(p): errors[2]
    errors_temporary.at(3) = err[1];    // H1-semi(p):  errors[3]
    //        double tsquare=(TDatabase::ParamDB->VISCOSITY/TDatabase::ParamDB->EFFECTIVE_VISCOSITY)*TDatabase::ParamDB->PERMEABILITY;
    //        errors_temporary.at(4) = sqrt(tsquare * errors_temporary[1]*errors_temporary[1] + errors_temporary[0]*errors_temporary[0]);
    //        Output::print<1>("t^2 * H1-semi(u) + L2(u): ", setprecision(14),  errors_temporary[4]);
    // Output::print<1>( setprecision(5),  "$",errors_temporary[0], "$ & $", errors_temporary[1], "$ & $", errors_temporary[2], "$ & $", errors_temporary[3], "$"& $", errors_temporary[4], "$");
     */ //LB 14.11.17 OLD ^


    Output::print("------------------------------------------------------");
    Output::print<1>(" L2(u):      ", setprecision(14), errors_temporary[0]);
    Output::print<1>(" L2(div(u)): ", setprecision(14), errors_temporary[1]);
    Output::print<1>(" H1-semi(u): ", setprecision(14), errors_temporary[2]);
    Output::print<1>(" L2(p):      ", setprecision(14), errors_temporary[3]);
    Output::print<1>(" H1-semi(p): ", setprecision(14), errors_temporary[4]);
    Output::print("------------------------------------------------------");

    // copy: This is necessary to save the calculated errors (errors_temporary) in the member variable errors.
    //       The vector errors_temporary is deleted with the end of the current scope and cannot be accessed from outside.
    std::copy(errors_temporary.begin(), errors_temporary.end()-1, this->errors.begin());
  }
  delete u1;
  delete u2;
}

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
double Brinkman2D::getH1SemiPressureError() const
{
    return this->errors[4];
}


/** ************************************************************************ */
std::array< double, int(6) > Brinkman2D::get_errors()
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

