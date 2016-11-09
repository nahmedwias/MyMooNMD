#include <Brinkman3D.h>
#include <MainUtilities.h> // get velocity and pressure space
#include <LinAlg.h> // DDot
#include <DirectSolver.h>
///#include <Assembler4.h>
#include <Assemble3D.h>
#include <Database.h>
#include <LocalAssembling3D.h>
#include <Output3D.h>
#include <sys/stat.h>

#ifdef _MPI
#include "mpi.h"
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif
#include <PETScSolver.h>

#include <BoundaryAssembling3D.h>



ParameterDatabase get_default_Brinkman3D_parameters()
{
    Output::print<5>("Creating a default Brinkman3D parameter database");
    
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.set_name("Brinkman3D parameter database");
    
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
Brinkman3D::System_per_grid::System_per_grid (const Example_Brinkman3D& example,
                                              TCollection& coll,
                                              std::pair<int,int> velocity_pressure_orders,
                                              Brinkman3D::Matrix type)
#ifdef _MPI
, int maxSubDomainPerDof
#endif
: velocity_space(&coll, (char*)"u", (char*)"Brinkman3D velocity", example.get_bc(0), //bd cond at 0 is x velo bc
                 velocity_pressure_orders.first),
pressure_space(&coll, (char*)"p", (char*)"Brinkman3D pressure", example.get_bc(3), //bd condition at 3 is pressure bc
               velocity_pressure_orders.second)
// TODO CB: Building the matrix here and rebuilding later is due to the
// highly non-functional class TFEVectFunction3D (and TFEFunction3D,
// which do neither provide default constructors nor working copy assignments.)
,matrix({&velocity_space, &velocity_space,&velocity_space, &pressure_space}),
rhs(matrix, true),
solution(matrix, false),
u(&velocity_space, (char*)"u", (char*)"u", solution.block(0),
  solution.length(0), 3),
p(&pressure_space, (char*)"p", (char*)"p", solution.block(3),
  solution.length(3))

{
    // rebuild the matrix due to NSE type. We must be sure, that the rhs and solution
    // vector which we built above do fit the new matrix, too!
    switch (type)
    {
        case Brinkman3D::Matrix::Type1:
            matrix = BlockFEMatrix::NSE3D_Type1(velocity_space, pressure_space);
            break;
        case Brinkman3D::Matrix::Type2:
            matrix = BlockFEMatrix::NSE3D_Type2(velocity_space, pressure_space);
            break;
        case Brinkman3D::Matrix::Type3:
            matrix = BlockFEMatrix::NSE3D_Type3(velocity_space, pressure_space);
            break;
        case Brinkman3D::Matrix::Type4:
            matrix = BlockFEMatrix::NSE3D_Type4(velocity_space, pressure_space);
            break;
        case Brinkman3D::Matrix::Type14:
            matrix = BlockFEMatrix::NSE3D_Type14(velocity_space, pressure_space);
            break;
        default:
            ErrThrow("Unknown Brinkman type given to constructor of Brinkman3D::System_per_grid.");
    }
    
#ifdef _MPI
    
    velocity_space.initialize_parallel(maxSubDomainPerDof);
    pressure_space.initialize_parallel(maxSubDomainPerDof);
    
    //print some information
    velocity_space.get_communicator().print_info();
    pressure_space.get_communicator().print_info();
    
#endif
    
}

/** ************************************************************************ */
void Brinkman3D::output_problem_size_info() const
{
    int my_rank = 0;
    
#ifndef _MPI
    const TFESpace3D & velocity_space = this->systems.front().velocity_space;
    const TFESpace3D & pressure_space = this->systems.front().pressure_space;
    
    size_t nDofu  = velocity_space.GetN_DegreesOfFreedom();
    size_t nDofp  = pressure_space.GetN_DegreesOfFreedom();
    size_t nTotal = 3*nDofu + nDofp;
    size_t nActive= 3*velocity_space.GetActiveBound();
    
    TCollection* coll = velocity_space.GetCollection();
    
    double hmin, hmax;
    coll->GetHminHmax(&hmin, &hmax);
#else
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    auto velocity_comm = systems.front().velocity_space.get_communicator();
    auto pressure_comm = systems.front().pressure_space.get_communicator();
    int nDofu  = velocity_comm.get_n_global_dof();
    int nDofp  = pressure_comm.get_n_global_dof();
    int nTotal = 3*nDofu + nDofp;
#endif
    
    if(my_rank ==0)
    {
        Output::stat("Brinkman3D", "Mesh data and problem size on finest grid");
#ifndef _MPI
        Output::dash("N_Cells      :  ", setw(10), coll->GetN_Cells());
        Output::dash("h(min, max)  :  ", setw(10), hmin, setw(10), " ", hmax);
#endif
        Output::dash("ndof velocity:  ", setw(10), 3*nDofu );
        Output::dash("ndof pressure:  ", setw(10), nDofp);
        Output::dash("ndof total   :  ", setw(10), nTotal );
#ifndef _MPI
        Output::dash("nActive      :  ", setw(10), nActive);
#endif
    }
}



/** ************************************************************************ */
Brinkman3D::Brinkman3D(const TDomain & domain,
                       const ParameterDatabase& param_db,
                       const Example_Brinkman3D& e,
                       int reference_id
#ifdef _MPI
                       , int maxSubDomainPerDof
#endif
)
: systems(), example(e), db(get_default_Brinkman3D_parameters()),
solver(param_db), defect(), norms_of_old_residuals(), initial_residual(1e10),
errors()

{
    this->db.merge(param_db,false);
    
    // set the argument to false for more detailed output on the console
    // db.info(true);
    
    std::pair <int,int> velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE,
                                                 TDatabase::ParamDB->PRESSURE_SPACE);
    
    // set the velocity and preesure spaces
    // this function returns a pair which consists of velocity and pressure order
    this->get_velocity_pressure_orders(velocity_pressure_orders);
    
//    // create the collection of cells from the domain (finest grid collection)
//    TCollection *coll = collections.front();
    
    // create the collection of cells from the domain (finest grid)
    TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
    
    // Output the mesh
    coll->writeMesh("test.mesh");
    //exit(1);
    
    // create finite element space and function, matrix, rhs, and solution
    switch(TDatabase::ParamDB->NSTYPE)
    {
        case 4:
#ifdef _MPI
            this->systems.emplace_back(e,
                                       *coll,
                                       velocity_pressure_orders,
                                       Brinkman3D::Matrix::Type4,
                                       maxSubDomainPerDof);
#else
            this->systems.emplace_back(e,
                                       *coll,
                                       velocity_pressure_orders,
                                       Brinkman3D::Matrix::Type4);
#endif
            break;
        case 14:
#ifdef _MPI
            this->systems.emplace_back(e,
                                       *coll,
                                       velocity_pressure_orders,
                                       Brinkman3D::Matrix::Type14,
                                       maxSubDomainPerDof);
#else
            this->systems.emplace_back(e,
                                       *coll,
                                       velocity_pressure_orders,
                                       Brinkman3D::Matrix::Type14);
#endif
            break;
    }
    
    
    output_problem_size_info();
    
    // the defect has the same structure as the rhs (and as the solution)
    this->defect.copy_structure(this->systems.front().rhs);
    
}


/** ************************************************************************ */
Brinkman3D::Brinkman3D(std::list<TCollection* > collections,
                       const ParameterDatabase& param_db,
                       const Example_Brinkman3D& e
#ifdef _MPI
                       , int maxSubDomainPerDof
#endif
)
: systems(), example(e), db(get_default_Brinkman3D_parameters()),
solver(param_db), defect(), norms_of_old_residuals(), initial_residual(1e10),
errors()

{
    this->db.merge(param_db,false);
    
    // set the argument to false for more detailed output on the console
    // db.info(true);
    
    std::pair <int,int> velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE,
                                                 TDatabase::ParamDB->PRESSURE_SPACE);
    
    // set the velocity and preesure spaces
    // this function returns a pair which consists of velocity and pressure order
    this->get_velocity_pressure_orders(velocity_pressure_orders);
    
    // create the collection of cells from the domain (finest grid collection)
    TCollection *coll = collections.front();
    
    // create finite element space and function, matrix, rhs, and solution
    switch(TDatabase::ParamDB->NSTYPE)
    {
        case 4:
#ifdef _MPI
            this->systems.emplace_back(example,
                                       *coll,
                                       velocity_pressure_orders,
                                       Brinkman3D::Matrix::Type4,
                                       maxSubDomainPerDof);
#else
            this->systems.emplace_back(example,
                                       *coll,
                                       velocity_pressure_orders,
                                       Brinkman3D::Matrix::Type4);
#endif
            break;
        case 14:
#ifdef _MPI
            this->systems.emplace_back(example,
                                       *coll,
                                       velocity_pressure_orders,
                                       Brinkman3D::Matrix::Type14,
                                       maxSubDomainPerDof);
#else
            this->systems.emplace_back(example,
                                       *coll,
                                       velocity_pressure_orders,
                                       Brinkman3D::Matrix::Type14);
#endif
    }
    
    output_problem_size_info();
    
    // the defect has the same structure as the rhs (and as the solution)
    this->defect.copy_structure(this->systems.front().rhs);
    
}

/** ************************************************************************ */
Brinkman3D::~Brinkman3D()
{
}

/** ************************************************************************ */
void Brinkman3D::get_velocity_pressure_orders(std::pair <int,int>
                                              &velocity_pressure_orders)
{
    Output::print<1>("velocity space", setw(10), TDatabase::ParamDB->VELOCITY_SPACE);
    Output::print<1>("pressure space", setw(10), TDatabase::ParamDB->PRESSURE_SPACE);
}

/** ************************************************************************ */
void Brinkman3D::set_parameters()
{
}

/** ************************************************************************ */
void Brinkman3D::assemble()
{
    for(System_per_grid& s : this->systems)
    {
        const TFESpace3D *v_space = &s.velocity_space;
        const TFESpace3D *p_space = &s.pressure_space;
        
        const TFESpace3D *fespmat[2] = {v_space, p_space};
        TSquareMatrix3D *sq_matrices[10]{nullptr}; // maximum number of square matrices (Type14)
        TMatrix3D *rect_matrices[6]{nullptr}; // it's four pointers maximum (Types 2, 4, 14)
        
        size_t N_FESpaces = 2; // spaces used for assembling matrices
        
        //        LocalAssembling3D_type type;
        
        //        if (db["P2P2_stab"].is(true))
        //        {type=Brinkman3D_Galerkin1ResidualStabP2;
        //        Output::print<>("P2P2 Stabilization");}
        //        else if (db["P1P1_stab"].is(true))
        //        {type=Brinkman3D_Galerkin1ResidualStabP1;
        //            Output::print<>("P1P1 Stabilization");}
        //        else
        //        {type=Brinkman3D_Galerkin1;}
        //
        //        LocalAssembling3D la(type, fe_functions,
        //                             this->example.get_coeffs());
        
        std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely();
        
        size_t n_sq_mat,n_rect_mat,N_Rhs;
        //right hand side must be adapted
        N_Rhs= 4;
        std::vector<double*> RHSs(N_Rhs); // right hand side
        //std::vector<double*> RHSs[4] = {s.rhs.block(0), s.rhs.block(1), s.rhs.block(2), nullptr}; // third place gets only filled
        s.rhs.reset();
        RHSs[0]=s.rhs.block(0);
        RHSs[1]=s.rhs.block(1);
        RHSs[2]=s.rhs.block(2);
        RHSs[3]=nullptr; //will be reset for type 4 and 14
        // Note: For now, we use only Type 4 and 14 for Brinkman
        //--------------------------------------------------------------------------------------------------
        // call the assemble method with the information that has been patched together
        // //old Assemble3D.C function
        switch(TDatabase::ParamDB->NSTYPE)
        {
            case 4:
                n_sq_mat = 9;
                sq_matrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
                sq_matrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks[1].get());
                sq_matrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks[2].get());
                sq_matrices[3]=reinterpret_cast<TSquareMatrix3D*>(blocks[4].get());
                sq_matrices[4]=reinterpret_cast<TSquareMatrix3D*>(blocks[5].get());
                sq_matrices[5]=reinterpret_cast<TSquareMatrix3D*>(blocks[6].get());
                sq_matrices[6]=reinterpret_cast<TSquareMatrix3D*>(blocks[8].get());
                sq_matrices[7]=reinterpret_cast<TSquareMatrix3D*>(blocks[9].get());
                sq_matrices[8]=reinterpret_cast<TSquareMatrix3D*>(blocks[10].get());
                
                n_rect_mat = 6;
                rect_matrices[0]=reinterpret_cast<TMatrix3D*>(blocks[12].get()); //first the lying B blocks
                rect_matrices[1]=reinterpret_cast<TMatrix3D*>(blocks[13].get());
                rect_matrices[2]=reinterpret_cast<TMatrix3D*>(blocks[14].get());
                rect_matrices[3]=reinterpret_cast<TMatrix3D*>(blocks[3].get()); //than the standing B blocks
                rect_matrices[4]=reinterpret_cast<TMatrix3D*>(blocks[7].get());
                rect_matrices[5]=reinterpret_cast<TMatrix3D*>(blocks[11].get());
                
                RHSs[3] = s.rhs.block(3); // NSE type 14 includes pressure rhs
                break;
            case 14:
                n_sq_mat = 10;
                sq_matrices[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
                sq_matrices[1] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
                sq_matrices[2] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
                sq_matrices[3] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
                sq_matrices[4] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
                sq_matrices[5] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
                sq_matrices[6] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
                sq_matrices[7] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());
                sq_matrices[8] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(10).get());
                sq_matrices[9] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(15).get());
                
                double *C=sq_matrices[9]->GetEntries();
                Output::print("Matrix C:");
                Output::print(C);
                
                n_rect_mat = 6;
                rect_matrices[0] = reinterpret_cast<TMatrix3D*>(blocks.at(12).get()); // first the lying B blocks
                rect_matrices[1] = reinterpret_cast<TMatrix3D*>(blocks.at(13).get());
                rect_matrices[2] = reinterpret_cast<TMatrix3D*>(blocks.at(14).get());
                rect_matrices[3] = reinterpret_cast<TMatrix3D*>(blocks.at(3).get());// than the standing B blocks
                rect_matrices[4] = reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
                rect_matrices[5] = reinterpret_cast<TMatrix3D*>(blocks.at(11).get());
                
                N_Rhs = 4; // is 4 if NSE type is 4 or 14 (else it is 3)
                RHSs[3] = s.rhs.block(3); // NSE type 14 includes pressure rhs
                break;
        }
        
        const TFESpace3D *fesprhs[4] = {v_space, v_space, v_space, p_space};  // if NSE type is 4 or 14
        
//???        for(unsigned int i=0; i<n_sq_mat; i++)
//            sq_matrices[i]->reset();
//        for(unsigned int i=0; i<n_rect_mat;i++)
//            rect_matrices[i]->reset();
        
        // boundary conditions and boundary values
        BoundCondFunct3D * boundary_conditions[4]={
            fespmat[0]->getBoundCondition(), fespmat[0]->getBoundCondition(),
            fespmat[0]->getBoundCondition(), fespmat[1]->getBoundCondition() };
        
        std::array<BoundValueFunct3D*, 4> non_const_bound_values;
        non_const_bound_values[0]=example.get_bd()[0];
        non_const_bound_values[1]=example.get_bd()[1];
        non_const_bound_values[2]=example.get_bd()[2];
        non_const_bound_values[3]=example.get_bd()[3];
        
        // finite element function used for nonlinear term
        std::vector<TFEFunction3D*> feFunction(4, nullptr);
        // finite element functions
        feFunction[0]=s.u.GetComponent(0);
        feFunction[1]=s.u.GetComponent(1);
        feFunction[2]=s.u.GetComponent(2);
        feFunction[3]=&s.p;
        
        // local assembling object
        const LocalAssembling3D la(LocalAssembling3D_type::Brinkman3D_Galerkin1,
                                   feFunction.data(),
                                   example.get_coeffs());
        
        // assemble now the matrices and right hand side
        Assemble3D(N_FESpaces, fespmat,
                   n_sq_mat, sq_matrices,
                   n_rect_mat, rect_matrices,
                   N_Rhs, RHSs.data(), fesprhs,
                   boundary_conditions,  non_const_bound_values.data(), la);
        
        sq_matrices[0]->write("a11.out");
        sq_matrices[4]->write("a22.out");
        sq_matrices[8]->write("a33.out");
        
        rect_matrices[0]->write("b1.out");
        rect_matrices[1]->write("b2.out");
        rect_matrices[2]->write("b3.out");
        rect_matrices[3]->write("b1t.out");
        rect_matrices[4]->write("b2t.out");
        rect_matrices[5]->write("b3t.out");

        //delete the temorary feFunctions gained by GetComponent
        for(int i = 0; i<3; ++i)
            delete feFunction[i];
        
    //////}// endfor auto grid
    
    //copy non-actives from rhs to solution on finest grid
    this->systems.front().solution.copy_nonactive(systems.front().rhs);
    
    /** When we call copy_nonactive in MPI-case, we have to remember the following:
     * it can happen that some slave ACTTIVE DoFs are placed in the block of
     * NON-ACTIVE DoFs (because they are at the interface between processors).
     * Doing copy_nonactive changes then the value of these DOFs,although they are
     * actually active.
     * That's why we have to update the values so that the vector becomes consistent again.
     * This is done here.
     */
#ifdef _MPI
    double *u1 = this->systems.front().solution.block(0);
    double *u2 = this->systems.front().solution.block(1);
    double *u3 = this->systems.front().solution.block(2);
    double *p  = this->systems.front().solution.block(3);
    this->systems.front().velocitySpace.get_communicator().consistency_update(u1, 3);
    this->systems.front().velocitySpace.get_communicator().consistency_update(u2, 3);
    this->systems.front().velocitySpace.get_communicator().consistency_update(u3, 3);
    this->systems.front().pressureSpace.get_communicator().consistency_update(p, 3);
#endif
    
    
    //--------------------------------------------------------------------------------------------------
    
    //        Assembler4 Ass(type, fe_functions, this->example.get_coeffs());
    // Brinkmann-specific choices
    //        std::vector<const TFESpace3D*> spaces_for_matrix;
    //        spaces_for_matrix.resize(2);
    //        spaces_for_matrix[0] = v_space;
    //        spaces_for_matrix[1] = p_space;
    //
    //        std::vector<const TFESpace3D*> spaces_for_rhs;
    //        spaces_for_rhs.resize(3);
    //        spaces_for_rhs[0] = v_space;
    //        spaces_for_rhs[1] = v_space;
    //        spaces_for_rhs[2] = p_space;
    //
    //        Ass.Assemble2D(s.matrix,s.rhs,
    //	             spaces_for_matrix,spaces_for_rhs,
    //	             example);
    // Weakly Imposing Boundary Conditions - Boundary Integrals
    

    
    TCollection* coll = v_space->GetCollection();
    std::vector<TBaseCell*> allCells;
    for (int i=0 ; i < coll->GetN_Cells(); i++)
      {
	ListOfAllCells[i]= coll->GetCell(i);
      }
    
    BoundaryAssembling3D bi;
    for (int k=0;k<TDatabase::ParamDB->n_neumann_boundary;k++)
    {
      bi.rhs_g_v_n(s.rhs,v_space,NULL,
		   allCells,            
		   -TDatabase::ParamDB->neumann_boundary_value[k]);
    }
    
    //        for (int k=0;k<TDatabase::ParamDB->n_g_v_boundary;k++)
    //        {
    //            bi.rhs_g_v(s.rhs,
    //                       v_space,
    //                       this->example.get_bd(0),                                 // access to U1BoundValue in the current example
    //                       this->example.get_bd(1),                                 // access to U2BoundValue in the current example
    //                       TDatabase::ParamDB->g_v_boundary_id[k],                  // boundary component
    //                       TDatabase::ParamDB->g_v_boundary_value[k],               // mult
    //                       false);                                                  // rescale local integral by edge values
    //        }
    //
    //        // test (in the true case, this function should be assembled on "DIRICHLET" boundaries
    //        // bi.matrix_v_n_v_n(s.matrix,v_space,1,1.);
    //
    //        for (int k=0;k<TDatabase::ParamDB->n_unvn_boundary;k++)
    //        {
    //            bi.matrix_v_n_v_n(s.matrix,
    //                              v_space,
    //                              TDatabase::ParamDB->unvn_boundary_id[k],          // boundary component
    //                              TDatabase::ParamDB->unvn_boundary_value[k]);      // mult
    //        }
    //
    //        for (int k=0;k<TDatabase::ParamDB->n_gradunv_boundary;k++)
    //        {
    //            bi.matrix_gradv_n_v(s.matrix,
    //                                v_space,
    //                                TDatabase::ParamDB->gradunv_boundary_id[k],     // boundary component
    //                                TDatabase::ParamDB->gradunv_boundary_value[k]); // mult
    //        }
    //
    //        for (int k=0;k<TDatabase::ParamDB->n_u_v_boundary;k++)
    //        {
    //            bi.matrix_u_v(s.matrix,
    //                          v_space,
    //                          TDatabase::ParamDB->u_v_boundary_id[k],               // boundary component
    //                          TDatabase::ParamDB->u_v_boundary_value[k],            // mult
    //                          false);                                               // rescale local integral by edge values
    //        }
    //
    //        for (int k=0;k<TDatabase::ParamDB->n_p_v_n_boundary;k++)
    //
    //        {
    //            bi.matrix_p_v_n(s.matrix,
    //                            v_space,
    //                            p_space,
    //                            TDatabase::ParamDB->p_v_n_boundary_id[k],           // boundary component
    //                            TDatabase::ParamDB->p_v_n_boundary_value[k]);       // mult
    //        }
    //
    //        // Nitsche combination - weak Dirichlet
    //        for (int k=0;k<TDatabase::ParamDB->n_nitsche_boundary;k++)
    //        {
    //            double K = TDatabase::ParamDB->PERMEABILITY;
    //            double nu = TDatabase::ParamDB->VISCOSITY;
    //            double nu_eff = TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
    //            double t = fabs(sqrt((nu_eff/nu)*K));
    //
    //            bi.matrix_gradv_n_v(s.matrix,
    //                                v_space,
    //                                TDatabase::ParamDB->nitsche_boundary_id[k],     // boundary component
    //                                -TDatabase::ParamDB->EFFECTIVE_VISCOSITY);      // mult
    //
    //            bi.matrix_p_v_n(s.matrix,
    //                            v_space,
    //                            p_space,
    //                            TDatabase::ParamDB->nitsche_boundary_id[k],         // boundary component
    //                            1.);                                                 // mult
    //
    //            bi.matrix_u_v(s.matrix,
    //                          v_space,
    //                          TDatabase::ParamDB->nitsche_boundary_id[k],           // boundary component
    //                          t*TDatabase::ParamDB->nitsche_penalty[k],               // mult
    //                          true);                                                // rescale local integral by edge values
    //
    //            bi.rhs_g_v(s.rhs,
    //                       v_space,
    //                       this->example.get_bd(0),                                 // access to U1BoundValue in the example,
    //                       this->example.get_bd(1),                                 // access to U2BoundValue in the example,
    //                       TDatabase::ParamDB->nitsche_boundary_id[k],              // boundary component
    //                       t*TDatabase::ParamDB->nitsche_penalty[k],                  // mult
    //                       true);                                                   // rescale local integral by edge values
    //        }
    
    //--------------------------------------------------------------------------------------------------
    }
}

/** ************************************************************************ */
bool Brinkman3D::stopIt(unsigned int iteration_counter)
{
}


/** ************************************************************************ */
void Brinkman3D::compute_norm_of_residual()
{
    //    System_per_grid& s = this->systems.front();
    //    unsigned int n_u_dof = s.solution.length(0);
    //    unsigned int n_p_dof = s.solution.length(3);
    //
    //    // copy rhs to defect
    //#ifdef _MPI
    //    //MPI: solution in consistency level 3 (TODO: maybe this is superfluous here
    //    // (because solution might be in level 3 consistency already)!)
    //    auto comms = s.matrix.get_communicators();
    //    for (size_t bl = 0; bl < comms.size() ;++bl)
    //    {
    //        comms[bl]->consistency_update(s.solution.block(bl), 3);
    //    }
    //#endif
    //
    //    this->defect = s.rhs;
    //    s.matrix.apply_scaled_add(s.solution, defect,-1.);
    //
    //    if( TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE )
    //    {
    //            IntoL20Vector3D(&defect[3*n_u_dof], n_p_dof,
    //                            TDatabase::ParamDB->PRESSURE_SPACE);
    //    }
    //
    //    // This is the calculation of the residual, given the defect.
    //    BlockVector defect_impuls({n_u_dof,n_u_dof,n_u_dof});
    //    BlockVector defect_mass(n_p_dof);
    //    //copy the entries (BlockVector offers no functionality to do this more nicely)
    //    for(size_t i = 0; i<3*n_u_dof ;++i)
    //        defect_impuls.get_entries()[i] = defect.get_entries()[i];
    //        for(size_t i =0 ; i<n_p_dof ; ++i)
    //            defect_mass.get_entries()[i] = defect.get_entries()[3*n_u_dof + i];
    //
    //#ifdef _MPI
    //            double impuls_residual_square = defect_impuls.norm_global({comms[0],comms[1],comms[2]});
    //            impuls_residual_square *= impuls_residual_square;
    //            double mass_residual_square = defect_mass.norm_global({comms[3]});
    //            mass_residual_square *= mass_residual_square;
    //#else
    //            double impuls_residual_square = defect_impuls.norm();
    //            impuls_residual_square *= impuls_residual_square;
    //            double mass_residual_square = defect_mass.norm();
    //            mass_residual_square *= mass_residual_square;
    //#endif
    //
    //            Residuals current_residuals(impuls_residual_square, mass_residual_square);
    //            norms_of_old_residuals.add(current_residuals);
}


/** ************************************************************************ */
void Brinkman3D::solve()
{
    System_per_grid& s = this->systems.front();
    
    /** @todo consider storing an object of DirectSolver in this class
     * DirectSolver direct_solver(s.matrix,
     * DirectSolver::DirectSolverTypes::umfpack);
     * direct_solver.solve(s.rhs, s.solution); */
    
    
    // solving:
#ifndef _MPI
    this->solver.solve(s.matrix, s.rhs, s.solution);
#endif
#ifdef _MPI
    if(this->solver.get_db()["solver_type"].is("direct"))
    {
        //set up a MUMPS wrapper
        MumpsWrapper mumps_wrapper(s.matrix_);
        
        //kick off the solving process
        mumps_wrapper.solve(s.rhs_, s.solution_);
    }
    else
        this->solver.solve(s.matrix, s.rhs, s.solution); // same as sequential
#endif
    
    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    {
        s.p.project_into_L20();
    }
}

/** ************************************************************************ */
void Brinkman3D::solve_with_Petsc(ParameterDatabase parmoon_db)
{
    System_per_grid& s = this->systems.front();
    
    /** @todo consider storing an object of DirectSolver in this class
     * DirectSolver direct_solver(s.matrix,
     * DirectSolver::DirectSolverTypes::umfpack);
     * direct_solver.solve(s.rhs, s.solution); */
    
//AENDERUNG
    PETScSolver PETScObject(s.matrix, parmoon_db);
    
    PETScObject.solve(s.rhs, s.solution);
//---------
    
////AENDERUNG
//    // solving:
//#ifndef _MPI
//    this->solver.solve(s.matrix, s.rhs, s.solution);
//#endif
//#ifdef _MPI
//    if(this->solver.get_db()["solver_type"].is("direct"))
//    {
//        //set up a MUMPS wrapper
//        MumpsWrapper mumps_wrapper(s.matrix);
//        
//        //kick off the solving process
//        mumps_wrapper.solve(s.rhs_, s.solution);
//    }
//    else
//        this->solver.solve(s.matrix, s.rhs, s.solution); // same as sequential
//#endif
////---------
    
    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    {
        s.p.project_into_L20();
    }
}



/** ************************************************************************ */
void Brinkman3D::output(int i)
{
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
    
    bool no_output = !db["output_write_vtk"] && !db["output_compute_errors"];
    if( no_output )
        return;
    
    System_per_grid& s = this->systems.front();
    TFEFunction3D* u1 = s.u.GetComponent(0);
    TFEFunction3D* u2 = s.u.GetComponent(1);
    TFEFunction3D* u3 = s.u.GetComponent(2);
    
    // print the value of the largest and smallest entry in the finite element vector
    if( (size_t)db["verbosity"]> 1 )
    {
        u1->PrintMinMax(std::string("u1"));
        u2->PrintMinMax(std::string("u2"));
        u3->PrintMinMax(std::string("u3"));
        s.p.PrintMinMax(std::string("p"));
    }
    
    // write solution to a vtk file
    if(db["output_write_vtk"])
    {
        // last argument in the following is domain, but is never used in this class
        TOutput3D Output(5, 5, 2, 1, NULL);
        Output.AddFEFunction(&s.p);
        Output.AddFEVectFunct(&s.u);
#ifdef _MPI
        char SubID[] = "";
        if(my_rank == 0)
            mkdir(db["output_vtk_directory"], 0777);
        std::string dir = db["output_vtk_directory"];
        std::string base = db["output_basename"];
        Output.Write_ParVTK(MPI_COMM_WORLD, 0, SubID, dir, base);
#else
        // Create output directory, if not already existing.
        mkdir(db["output_vtk_directory"], 0777);
        std::string filename = this->db["output_vtk_directory"];
        filename += "/" + this->db["output_basename"].value_as_string();
        
        if(i >= 0)
            filename += "_" + std::to_string(i);
        filename += ".vtk";
        Output.WriteVtk(filename.c_str());
#endif
    }
    
    // measure errors to known solution
    // If an exact solution is not known, it is usually set to be zero, so that
    // in such a case here only integrals of the solution are computed.
    if( db["output_compute_errors"] )
    {
        double err_u1[4]; // of these arrays only the two first entries are used,
        double err_u2[4]; // but the evil GetErrors() will corrupt memory if these
        double err_u3[4]; // have not at least size 4
        double err_p[4];
        
        TAuxParam3D aux(1, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
        MultiIndex3D nsAllDerivs[4] = {D000, D100, D010, D001};
        const TFESpace3D *velocity_space = &this->get_velocity_space();
        const TFESpace3D *pressure_space = &this->get_pressure_space();
        
        // errors in first velocity component
        u1->GetErrors(example.get_exact(0), 4, nsAllDerivs, 2,
                      L2H1Errors, nullptr, &aux, 1, &velocity_space, err_u1);
        // errors in second velocity component
        u2->GetErrors(example.get_exact(1), 4, nsAllDerivs, 2,
                      L2H1Errors, nullptr, &aux, 1, &velocity_space, err_u2);
        // errors in third velocity component
        u3->GetErrors(example.get_exact(2), 4, nsAllDerivs, 2,
                      L2H1Errors, nullptr, &aux, 1, &velocity_space, err_u3);
        // errors in pressure
        s.p.GetErrors(example.get_exact(3), 4, nsAllDerivs, 2, L2H1Errors,
                      nullptr, &aux, 1, &pressure_space, err_p);
        
#ifdef _MPI
        double err_red[8]; //memory for global (across all processes) error
        double err_send[8]; //fill send buffer
        err_send[0]=err_u1[0];
        err_send[1]=err_u1[1];
        err_send[2]=err_u2[0];
        err_send[3]=err_u2[1];
        err_send[4]=err_u3[0];
        err_send[5]=err_u3[1];
        err_send[6]=err_p[0];
        err_send[7]=err_p[1];
        
        MPI_Allreduce(err_send, err_red, 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for(i=0;i<8;i++)
        {//MPI: sqrt was skipped in GetErrors function - do it here globally!
            err_red[i] = sqrt(err_red[i]);
        }
        //fill the reduced errors back where they belong
        err_u1[0] = err_red[0];
        err_u1[1] = err_red[1];
        err_u2[0] = err_red[2];
        err_u2[1] = err_red[3];
        err_u3[0] = err_red[4];
        err_u3[1] = err_red[5];
        err_p[0] = err_red[6];
        err_p[1] = err_red[7];
#else
        int my_rank =0;
#endif
        
        errors.at(0) = sqrt(err_u1[0]*err_u1[0] + err_u2[0]*err_u2[0] + err_u3[0]*err_u3[0]);//L2
        errors.at(1) = sqrt(err_u1[1]*err_u1[1] + err_u2[1]*err_u2[1] + err_u3[1]*err_u3[1]);//H1-semi
        errors.at(2) = err_p[0];
        errors.at(3) = err_p[1];
        
        //print errors
        if(my_rank == 0)
        {
            Output::stat("Brinkman3D", "Measured errors");
            Output::dash("L2(u)     : ", setprecision(10), errors.at(0));
            Output::dash("H1-semi(u): ", setprecision(10), errors.at(1));
            Output::dash("L2(p)     : ", setprecision(10), errors.at(2));
            Output::dash("H1-semi(p): ", setprecision(10), errors.at(3));
        }
    }// if(this->db["compute_errors"])
    delete u1;
    delete u2;
    delete u3;
    
    //do postprocessing step depending on what the example implements
    example.do_post_processing(*this);
}


/** ************************************************************************ */

TFEFunction3D* Brinkman3D::get_velocity_component(int i)
{
    if(i==0)
        return this->systems.front().u.GetComponent(0);
    else if(i==1)
        return this->systems.front().u.GetComponent(1);
    else  if(i==2)
        return this->systems.front().u.GetComponent(2);
    else
        throw std::runtime_error("There are only three velocity components!");
}

const Residuals& Brinkman3D::get_residuals() const
{
    return norms_of_old_residuals.back();
}

double Brinkman3D::get_impuls_residual() const
{
    return norms_of_old_residuals.back().impulsResidual;
}

double Brinkman3D::get_mass_residual() const
{
    return norms_of_old_residuals.back().massResidual;
}

double Brinkman3D::get_full_residual() const
{
    return norms_of_old_residuals.back().fullResidual;
}

std::array<double, int(4)> Brinkman3D::get_errors() const
{
    return errors;
}

