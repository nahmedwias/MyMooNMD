#include <Brinkman2D.h>
#include <MainUtilities.h> // get velocity and pressure space
#include <LinAlg.h> // DDot
#include <DirectSolver.h>
#include <Assembler4.h>
#include <Assemble2D.h>
#include <BoundaryAssembling2D.h>

ParameterDatabase get_default_Brinkman2D_parameters()
{
    Output::print<5>("creating a default Brinkman2D parameter database");
    
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.set_name("Brinkman2D parameter database");
    
    ParameterDatabase out_db = ParameterDatabase::default_output_database();
    db.merge(out_db, true);
    
    db.add("P1P1_stab", false,
           "Use an assembling routine corresponding to a residual-based "
           "equal-order stabilization for the Brinkman problem."
           "This only works in two space "
           "dimensions and is meaningfull for the finite elemnt space P1/P1 only."
           "Usually this is used in the main program.",
           {true,false});
    
    db.add("P2P2_stab", false,
           "Use an assembling routine corresponding to a residual-based "
           "equal-order stabilization for the Brinkman problem."
           "This only works in two space "
           "dimensions and is meaningfull for the finite elemnt spaces P1/P1 and P2/P2 only."
           "Usually this is used in the main program.",
           {true,false});
    
//    db.add("equal_order_stab_weight", 0,
//           "Use an assembling routine corresponding to a residual-based "
//           "equal-order stabilization for the Brinkman problem."
//           "This only works in two space "
//           "dimensions and is meaningfull for the finite elemnt spaces P1/P1 and P2/P2 only."
//           "Usually this is used in the main program.",-1000,1000);
    return db;
}



/** ************************************************************************ */
Brinkman2D::System_per_grid::System_per_grid (const Example_Brinkman2D& example,
                                              TCollection& coll, std::pair<int,int> velocity_pressure_orders,
                                              Brinkman2D::Matrix type)
: velocity_space(&coll, (char*)"u", (char*)"Brinkman velocity", example.get_bc(0),
                 velocity_pressure_orders.first, nullptr),
pressure_space(&coll, (char*)"p", (char*)"Brinkman pressure", example.get_bc(2),
               velocity_pressure_orders.second, nullptr)
// TODO CB: Building the matrix here and rebuilding later is due to the
// highly non-functional class TFEVectFunction2D (and TFEFunction2D,
// which do neither provide default constructors nor working copy assignments.)
,matrix({&velocity_space, &velocity_space, &pressure_space}),
rhs(matrix, true),
solution(matrix, false),
u(&velocity_space, (char*)"u", (char*)"u", solution.block(0),
  solution.length(0), 2),
p(&pressure_space, (char*)"p", (char*)"p", solution.block(2),
  solution.length(2))
{
    // rebuild the matrix due to NSE type. We must be sure, that the rhs and solution
    // vector which we built above do fit the new matrix, too!
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
: db(get_default_Brinkman2D_parameters()), solver(param_db),
  outputWriter(param_db),
  systems(), example(e), defect(), oldResiduals(), 
initial_residual(1e10), errors()
{
    this->db.merge(param_db,false);
    
    // set the argument to false for more detailed output on the console
//    Output::print<>("HIIIIIIEEEEERRRRRRRRRR");
//    db.info(true);
    
    std::pair <int,int> velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE,
                                                 TDatabase::ParamDB->PRESSURE_SPACE);
    
    // set the velocity and preesure spaces
    // this function returns a pair which consists of velocity and pressure order
    this->get_velocity_pressure_orders(velocity_pressure_orders);
    
    // create the collection of cells from the domain (finest grid)
    TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
    
    // we always use Matrix Type 14
    this->systems.emplace_back(example, *coll, velocity_pressure_orders, Brinkman2D::Matrix::Type14);
    
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
    Output::print<1>("N_Cells            : ", setw(10), coll->GetN_Cells());
    Output::print<1>("h (min,max)        : ", setw(10), h_min, " ", setw(12), h_max);
    Output::print<1>("dof velocity       : ", setw(10), 2* n_u);
    Output::print<1>("dof velocity active: ", setw(10), 2* n_u_active);
    Output::print<1>("dof pressure       : ", setw(10), n_p);
    Output::print<1>("dof all            : ", setw(10), n_dof);
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
void Brinkman2D::assemble()
{
    for(System_per_grid& s : this->systems)
    {
        s.rhs.reset(); //right hand side reset (TODO: is that necessary?)
        
        const TFESpace2D * v_space = &s.velocity_space;
        const TFESpace2D * p_space = &s.pressure_space;
        
        // declare the variables which Assemble2D needs and each brinkman type has to fill
        size_t N_FESpaces = 2;
        
        const TFESpace2D *fespmat[2] = {v_space, p_space};
        TSquareMatrix2D *sq_matrices[5]{nullptr}; // it's five pointers maximum (Type14)
        
        TMatrix2D *rect_matrices[4]{nullptr}; // it's four pointers maximum (Types 2, 4, 14)

        BoundCondFunct2D * boundary_conditions[3] = {
            v_space->GetBoundCondition(), v_space->GetBoundCondition(),
            p_space->GetBoundCondition() };

        std::array<BoundValueFunct2D*, 3> non_const_bound_values;
        non_const_bound_values[0] = example.get_bd()[0];
        non_const_bound_values[1] = example.get_bd()[1];
        non_const_bound_values[2] = example.get_bd()[2];
        
        //same for all: the local asembling object
        TFEFunction2D *fe_functions[3] =
        { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
        
        LocalAssembling2D_type type;
        
        if (db["P2P2_stab"].is(true))
        {type=Brinkman2D_Galerkin1ResidualStabP2;
        Output::print<>("P2P2 Stabilization");}
        else if (db["P1P1_stab"].is(true))
        {type=Brinkman2D_Galerkin1ResidualStabP1;
            Output::print<>("P1P1 Stabilization");}
        else
        {type=Brinkman2D_Galerkin1;}
        
        LocalAssembling2D la(type, fe_functions,
                             this->example.get_coeffs());
        
        std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely();
        
        // Note: We use only Type 14 for Brinkman (for now)
        //--------------------------------------------------------------------------------------------------
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
	
	
	/*Assemble2D(N_FESpaces, fespmat, n_sq_mat, sq_matrices,
		   n_rect_mat, rect_matrices, N_Rhs, RHSs, fesprhs,
		   boundary_conditions, non_const_bound_values.data(), la);
	*/
        //--------------------------------------------------------------------------------------------------

	// use a list of LocalAssembling2D
	std::vector< LocalAssembling2D* > la_list;
	la_list.push_back(&la);
	
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

        //--------------------------------------------------------------------------------------------------
        // Weakly Imposing Boundary Conditions - Boundary Integrals
        
        BoundaryAssembling2D bi;
        
        for (int k=0;k<TDatabase::ParamDB->n_neumann_boundary;k++)
        {
            bi.rhs_g_v_n(s.rhs,
                         v_space,
                         NULL,                                                  // g = 1
                         TDatabase::ParamDB->neumann_boundary_id[k],            // boundary component
                         -TDatabase::ParamDB->neumann_boundary_value[k]);       // mult
        }
        
        for (int k=0;k<TDatabase::ParamDB->n_g_v_boundary;k++)
        {
            bi.rhs_g_v(s.rhs,
                       v_space,
                       this->example.get_bd(0),                                 // access to U1BoundValue in the current example
                       this->example.get_bd(1),                                 // access to U2BoundValue in the current example
                       TDatabase::ParamDB->g_v_boundary_id[k],                  // boundary component
                       TDatabase::ParamDB->g_v_boundary_value[k],               // mult
                       false);                                                  // rescale local integral by edge values
        }
        
        // test (in the true case, this function should be assembled on "DIRICHLET" boundaries
        // bi.matrix_v_n_v_n(s.matrix,v_space,1,1.);
        
        for (int k=0;k<TDatabase::ParamDB->n_unvn_boundary;k++)
        {
            bi.matrix_v_n_v_n(s.matrix,
                              v_space,
                              TDatabase::ParamDB->unvn_boundary_id[k],          // boundary component
                              TDatabase::ParamDB->unvn_boundary_value[k]);      // mult
        }
        
        for (int k=0;k<TDatabase::ParamDB->n_gradunv_boundary;k++)
        {
            bi.matrix_gradv_n_v(s.matrix,
                                v_space,
                                TDatabase::ParamDB->gradunv_boundary_id[k],     // boundary component
                                TDatabase::ParamDB->gradunv_boundary_value[k]); // mult
        }
        
        for (int k=0;k<TDatabase::ParamDB->n_u_v_boundary;k++)
        {
            bi.matrix_u_v(s.matrix,
                          v_space,
                          TDatabase::ParamDB->u_v_boundary_id[k],               // boundary component
                          TDatabase::ParamDB->u_v_boundary_value[k],            // mult
                          false);                                               // rescale local integral by edge values
        }
        
        for (int k=0;k<TDatabase::ParamDB->n_p_v_n_boundary;k++)
            
        {
            bi.matrix_p_v_n(s.matrix,
                            v_space,
                            p_space,
                            TDatabase::ParamDB->p_v_n_boundary_id[k],           // boundary component
                            TDatabase::ParamDB->p_v_n_boundary_value[k]);       // mult
        }
        
        // Nitsche combination - weak Dirichlet
        for (int k=0;k<TDatabase::ParamDB->n_nitsche_boundary;k++)
        {
            double K = TDatabase::ParamDB->PERMEABILITY;
            double nu = TDatabase::ParamDB->VISCOSITY;
            double nu_eff = TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
            double t = fabs(sqrt((nu_eff/nu)*K));
            
            bi.matrix_gradv_n_v(s.matrix,
                                v_space,
                                TDatabase::ParamDB->nitsche_boundary_id[k],     // boundary component
                                -TDatabase::ParamDB->EFFECTIVE_VISCOSITY);      // mult
            
            bi.matrix_p_v_n(s.matrix,
                            v_space,
                            p_space,
                            TDatabase::ParamDB->nitsche_boundary_id[k],         // boundary component
                            1.);                                                // mult
            
            bi.matrix_u_v(s.matrix,
                          v_space,
                          TDatabase::ParamDB->nitsche_boundary_id[k],           // boundary component
                          t*TDatabase::ParamDB->nitsche_penalty[k],             // mult
                          true);                                                // rescale local integral by edge values
            
            bi.rhs_g_v(s.rhs,
                       v_space,
                       this->example.get_bd(0),                                 // access to U1BoundValue in the example,
                       this->example.get_bd(1),                                 // access to U2BoundValue in the example,
                       TDatabase::ParamDB->nitsche_boundary_id[k],              // boundary component
                       t*TDatabase::ParamDB->nitsche_penalty[k],                // mult
                       true);                                                   // rescale local integral by edge values
        }
        
        //--------------------------------------------------------------------------------------------------
        
        // copy Dirichlet values from right hand side into solution
        s.solution.copy_nonactive(s.rhs);
        
        // TODO Maybe we have to explicitely set non-actives in non-diagonal blocks
        // to zero here, that was done in former code, but maybe we can move it to the solver part
        
        // tidy up
        delete fe_functions[0];
        delete fe_functions[1];
	Output::print("end-2");
    }
    Output::print("end-2");
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
    
    if( TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE )
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
    
    /// @todo consider storing an object of DirectSolver in this class
    DirectSolver direct_solver(s.matrix,
			       DirectSolver::DirectSolverTypes::umfpack);
    direct_solver.solve(s.rhs, s.solution);
    
    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
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
    /*  // write solution to a vtk file
     if(db["output_write_vtk"])
     {
     outputWriter.write(i);
     }
     */
    
    // measure errors to known solution
    // If an exact solution is not known, it is usually set to be zero, so that
    // in such a case here only integrals of the solution are computed.
    if( db["output_compute_errors"] )
    {
        double err[4];
        TAuxParam2D aux_error;
        MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
        const TFESpace2D *velocity_space = &this->get_velocity_space();
        const TFESpace2D *pressure_space = &this->get_pressure_space();
        
        
        // errors in first velocity component
        u1->GetErrors(example.get_exact(0), 3, AllDerivatives, 2, L2H1Errors,
                      nullptr, &aux_error, 1, &velocity_space, err);
        // errors in second velocity component
        u2->GetErrors(example.get_exact(1), 3, AllDerivatives, 2, L2H1Errors,
                      nullptr, &aux_error, 1, &velocity_space, err + 2);
        
        errors.at(0) = sqrt(err[0]*err[0] + err[2]*err[2]);     // err[0]=L2-Error in u1, err[2]=L2-Error in u2
        errors.at(1) = sqrt(err[1]*err[1] + err[3]*err[3]);     // err[1]=H1-Error in u1, err[3]=H1-Error in u2
        Output::print<1>("L2(u)     : ", setprecision(14), errors[0]);
        Output::print<1>("H1-semi(u): ", setprecision(14), errors[1]);
        
        // errors in pressure
        s.p.GetErrors(example.get_exact(2), 3, AllDerivatives, 2, L2H1Errors,
                      nullptr, &aux_error, 1, &pressure_space, err);
        
        
        errors.at(2) = err[0];
        errors.at(3) = err[1];
        Output::print<1>("L2(p)     : ", setprecision(14), errors[2]);
        Output::print<1>("H1-semi(p): ", setprecision(14),  errors[3]);
        
        
    } // if(TDatabase::ParamDB->MEASURE_ERRORS)
    delete u1;
    delete u2;
}

/** ************************************************************************ */
double Brinkman2D::getL2VelocityError() const
{
    return this->errors[0];
}

/** ************************************************************************ */
//double Brinkman2D::getL2DivergenceError() const
//{
//    return this->errors[1];
//}

/** ************************************************************************ */
double Brinkman2D::getH1SemiVelocityError() const
{
    return this->errors[1];
}

/** ************************************************************************ */
double Brinkman2D::getL2PressureError() const
{
    return this->errors[2];
}

/** ************************************************************************ */
double Brinkman2D::getH1SemiPressureError() const
{
    return this->errors[3];
}


/** ************************************************************************ */
std::array< double, int(4) > Brinkman2D::get_errors()
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

