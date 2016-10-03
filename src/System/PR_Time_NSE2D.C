#include <PR_Time_NSE2D.h>
#include <Database.h>
#include <LocalAssembling2D.h>
#include <SquareMatrix2D.h>
#include <Assemble2D.h>

using namespace std;

/// used for the matrix vector multiplications
extern "C" void dgemm_(char *TRANSA, char *TRANSB, int *m, int *n, int *k,
           double *alpha, double *A, int *lda, double *B,
           int *ldb, double *beta, double *C, int *ldc);

ParameterDatabase get_default_PR_TNSE2D_parameters()
{
  Output::print<5>("creating a default Pressure Robust TNSE2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default TNSE2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("Pressure Robust TNSE2D parameter database");

  //Time_NSE2D requires a nonlinear iteration, set up a nonlinit_database and merge
  ParameterDatabase nl_db = ParameterDatabase::default_nonlinit_database();
  db.merge(nl_db,true);

  // a default output database - needed here as long as there's no
  // class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  return db;
}
PR_Time_NSE2D::system_per_grid::system_per_grid(const TFESpace2D &velsp, const TFESpace2D &prsp,
                          const Example_TimeNSE2D &ex, int order)
    : projSpace_(velsp.GetCollection(), (char*)"u", (char*)"projection space",
                 ex.get_bc(2), order, nullptr),
      projMatrix_({&velsp}, {&projSpace_, &projSpace_, &prsp}),
      modifedMassMatrix_({&velsp, &velsp})
{
  projMatrix_.print_coloring_pattern("projection matrix");
  modifedMassMatrix_.print_coloring_pattern("modified mass matrix");

  if(TDatabase::ParamDB->NSTYPE !=4 && TDatabase::ParamDB->DISCTYPE==300)
  {
    ErrThrow("pressure robust FEM only implemented for NSTYPE 4");
  }

  projMatrix_=BlockFEMatrix::Projection_NSE2D(velsp,projSpace_,prsp);

  modifedMassMatrix_=BlockFEMatrix::BlockFEMatrix::NSE2D_Type4(velsp, prsp);
}


PR_Time_NSE2D::PR_Time_NSE2D(const TDomain &domain,
                             const ParameterDatabase &param_db,
                             const Example_TimeNSE2D &ex, int ref_id)
    : Time_NSE2D(domain, param_db, ex, ref_id) ,
     db_(get_default_PR_TNSE2D_parameters()), systems_(),
     example_(ex), solver_(param_db)
{
  db_.merge(param_db);
  this->set_projection_space();
  int order=TDatabase::ParamDB->PROJECTION_SPACE;

  if(TDatabase::ParamDB->DISCTYPE == RECONSTRUCTION)
  {
    if(TDatabase::ParamDB->NSTYPE != 4)
    {
      ErrThrow("TDatabase::ParamDB->NSTYPE = ", 4 ,
               " is only supported in Pressre-Robust FEM.");
    }
  }
  bool usingMultigrid = this->solver.is_using_multigrid();

  if(!usingMultigrid)
  {
    // create the collection of cells from the domain (finest grid)
    this->systems_.emplace_back(this->Time_NSE2D::get_velocity_space(),
                          this->Time_NSE2D::get_pressure_space(),
                          example_, order);

    TFEFunction2D * u1 = this->systems.front().u.GetComponent(0);
    TFEFunction2D * u2 = this->systems.front().u.GetComponent(1);

    u1->Interpolate(example_.get_initial_cond(0));
    u2->Interpolate(example_.get_initial_cond(1));
  }
  else
  {
    auto multigrid = this->solver_.get_multigrid();
    // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;
    size_t n_levels = multigrid->get_n_geometric_levels();
    int finest = domain.get_ref_level();
    int coarsest = finest - n_levels + 1;
    for (int grid_no = finest; grid_no >= coarsest; --grid_no)
    {
      systems_.emplace_back(this->Time_NSE2D::get_velocity_space(),
                            this->Time_NSE2D::get_pressure_space(),
                            example_, order);
      //prepare input argument for multigrid object
      matrices.push_front(&systems.back().matrix);
    }
    multigrid->initialize(matrices);
  }
  this->print_info();
}

void PR_Time_NSE2D::set_projection_space() const
{
  switch(TDatabase::ParamDB->VELOCITY_SPACE)
  {
     case 2: // BDM2
       TDatabase::ParamDB->PROJECTION_SPACE = 1012;
       break;
     case 22:
       TDatabase::ParamDB->PROJECTION_SPACE = 1012;
       break;
     case 3:
       TDatabase::ParamDB->PROJECTION_SPACE = 1013;
       break;
     case 4:
       TDatabase::ParamDB->PROJECTION_SPACE = 1014;
       break;
  }
  Output::print("projection space ", setw(10),
                TDatabase::ParamDB->PROJECTION_SPACE);
}

void PR_Time_NSE2D::print_info() const
{
  this->Time_NSE2D::output_problem_size_info();
  int n_dof_pr=this->get_projection_space().GetN_DegreesOfFreedom();
  int n_active=this->get_projection_space().GetActiveBound();
  Output::print("dof RT/BDM        : ", setw(10), n_dof_pr);
  Output::print("dof active RT/BDM : ", setw(10), n_active);

  if(TDatabase::ParamDB->DISCTYPE == RECONSTRUCTION)
    Output::print("Pressure Robust FEM");
  else
    Output::print("      Classical FEM");
}

void PR_Time_NSE2D::assemble_initial()
{
  /// base class systems
  System_per_grid& sb=this->Time_NSE2D::systems.front();
  this->Time_NSE2D::assemble_initial_time();
  if(TDatabase::ParamDB->DISCTYPE != RECONSTRUCTION)
  {
    return;
  }
  /// assemble the initial matrices and rhs side
  /// depending on the disctypes

  if(TDatabase::ParamDB->DISCTYPE==RECONSTRUCTION)
  {
    LocalAssembling2D_type la_type;
    if(db_["problem_type"].is(4))
      la_type=RECONSTR_TSTOKES;
    else if(db_["problem_type"].is(6))
      la_type=RECONSTR_TNSE;
    else
      ErrThrow("problem type ", db_["problem_type"], " is not supported");
    Output::print("assembling the initial matrices ");
    for(system_per_grid& sd : this->systems_)
    {
      TFEFunction2D* fef[3]={sb.u.GetComponent(0),sb.u.GetComponent(1),&sb.p};
      LocalAssembling2D la(la_type, fef, example_.get_coeffs());
      // setting the space for sign in GetSignOfThisDOF();
      size_t vel_space = TDatabase::ParamDB->VELOCITY_SPACE;
      TDatabase::ParamDB->VELOCITY_SPACE = TDatabase::ParamDB->PROJECTION_SPACE;
      // variables used to assemble
      size_t nfeSp=2;
      const TFESpace2D* vel_sp=&sb.velocity_space;
      const TFESpace2D* pro_sp=&sd.projSpace_;
      const TFESpace2D* pointer_to_space[2]={vel_sp, pro_sp};
      // matrices to be assembled
      size_t nsq_ass = 3;
      size_t nre_ass = 4;
      // row and column spaces
      std::vector<int> rspace={0,0,1,1,1,0,0}; // 0: velocity space
      std::vector<int> cspace={0,0,1,0,0,1,1}; // 1: projection space
      // prepare everything that will be stored
      size_t nsq_stor=8;
      TSquareMatrix2D* sq_mat_store[8]{nullptr};
      std::vector<std::shared_ptr<FEMatrix>> mm
             = sd.modifedMassMatrix_.get_blocks_uniquely();
      std::vector<std::shared_ptr<FEMatrix>> mAB
             = sb.matrix.get_blocks_uniquely();
      sq_mat_store[0]=reinterpret_cast<TSquareMatrix2D*>(mm.at(0).get());
      sq_mat_store[1]=reinterpret_cast<TSquareMatrix2D*>(mm.at(1).get());
      sq_mat_store[2]=reinterpret_cast<TSquareMatrix2D*>(mm.at(3).get());
      sq_mat_store[3]=reinterpret_cast<TSquareMatrix2D*>(mm.at(4).get());

      sq_mat_store[4]=reinterpret_cast<TSquareMatrix2D*>(mAB.at(0).get());
      sq_mat_store[5]=reinterpret_cast<TSquareMatrix2D*>(mAB.at(1).get());
      sq_mat_store[6]=reinterpret_cast<TSquareMatrix2D*>(mAB.at(3).get());
      sq_mat_store[7]=reinterpret_cast<TSquareMatrix2D*>(mAB.at(4).get());
      // reset the matrices
      for(size_t i=0; i<nsq_stor; i++)
        sq_mat_store[i]->reset();
      // rectangular matrices to be stored
      size_t nre_stor=0;
      TMatrix2D** re_mat_store=nullptr;

      // right hand side to be assembled
      size_t nrhs_ass = 1;
      // row space for the right hand side
      std::vector<int> rspace_rhs={1};

      const TFESpace2D* pre_sp=&sb.pressure_space;
      // right hand side to be stored
      size_t nrhs_sto=2;
      // pressure block is empty always
      sb.rhs.reset();
      double *rhs_stor[3]={sb.rhs.block(0), sb.rhs.block(1),sb.rhs.block(2)};
      const TFESpace2D* pointer_to_rhs_space[3]={vel_sp,vel_sp,pre_sp};

      // boundary conditions
      BoundCondFunct2D *bc[3] ={ vel_sp->GetBoundCondition(),
                                 vel_sp->GetBoundCondition(),
                                 pre_sp->GetBoundCondition() };
      // boundary values
      std::array<BoundValueFunct2D*, 3> bv;
      bv[0] = example_.get_bd()[0];
      bv[1] = example_.get_bd()[1];
      bv[2] = example_.get_bd()[2];
      // assemble the correspoding matrices and right hand side
      Assemble2D_MixedFEM(nfeSp, pointer_to_space,
                 nsq_ass, nre_ass, rspace, cspace, nrhs_ass, rspace_rhs,
                 nsq_stor, sq_mat_store, nre_stor, re_mat_store,
                 nrhs_sto, rhs_stor,pointer_to_rhs_space, la,
                 matrices_reconstruction, nullptr,
                 ProjectionMatricesTNSE2D, bc, bv.data());

      TDatabase::ParamDB->VELOCITY_SPACE = vel_space;
    }
    // copy the solution and rhs vectors for the next time steps
    this->oldsol_=sb.solution;
    this->oldrhs_=sb.rhs;
  }
}

void PR_Time_NSE2D::rhs_assemble()
{
  Output::print<5>("assembling the system right hand side ");
  // System_per_grid for the member access from the Time_NSE2D class
  System_per_grid& sb=this->Time_NSE2D::systems.front();
  // System_per_grid for the member access from current class
  system_per_grid& sd=this->systems_.front();
  if(TDatabase::ParamDB->DISCTYPE==RECONSTRUCTION)
  {
    TFEFunction2D* fef[3]={sb.u.GetComponent(0),sb.u.GetComponent(1), &sb.p};
    LocalAssembling2D larhs(RECONSTR_GALERKIN_Rhs, fef,example_.get_coeffs());
    size_t vsp=TDatabase::ParamDB->VELOCITY_SPACE;
    size_t psp=TDatabase::ParamDB->PROJECTION_SPACE;
    TDatabase::ParamDB->VELOCITY_SPACE=psp;
    size_t nfeSp=2;
    const TFESpace2D* vel_sp=&sb.velocity_space;
    const TFESpace2D* pro_sp=&sd.projSpace_;
    const TFESpace2D* pointer_to_space[2]={vel_sp, pro_sp};
    // matrices to be assembled
    size_t nsq_ass = 0;
    size_t nre_ass = 2;
    // row and column spaces
    std::vector<int> rspace={0,0}; // 0: velocity space
    std::vector<int> cspace={1,1}; // 1: projection space

    // prepare everything that will be stored
    size_t nsq_stor=0;
    TSquareMatrix2D** sq_mat_store={nullptr};
    // rectangular matrices to be stored
    size_t nre_stor=0;
    TMatrix2D** re_mat_store=nullptr;

    // right hand side to be assembled
    size_t nrhs_ass = 1;
    // row space for the right hand side
    std::vector<int> rspace_rhs={1};

    const TFESpace2D* pre_sp=&sb.pressure_space;
    // right hand side to be stored
    size_t nrhs_sto=2;
    // pressure block is empty always
    sb.rhs.reset();
    double *rhs_stor[3]={sb.rhs.block(0), sb.rhs.block(1),sb.rhs.block(2)};
    const TFESpace2D* pointer_to_rhs_space[3]={vel_sp,vel_sp,pre_sp};
    // boundary conditions
    BoundCondFunct2D *bc[3] ={ vel_sp->GetBoundCondition(),
                               vel_sp->GetBoundCondition(),
                               pre_sp->GetBoundCondition() };
    // boundary values
    std::array<BoundValueFunct2D*, 3> bv;
    bv[0] = example_.get_bd()[0];
    bv[1] = example_.get_bd()[1];
    bv[2] = example_.get_bd()[2];

    Assemble2D_MixedFEM(nfeSp, pointer_to_space,
             nsq_ass, nre_ass, rspace, cspace, nrhs_ass, rspace_rhs,
             nsq_stor, sq_mat_store, nre_stor, re_mat_store, nrhs_sto,
             rhs_stor,pointer_to_rhs_space, larhs, nullptr, MatVectMult,
             projection_matrices, bc, bv.data());

    TDatabase::ParamDB->VELOCITY_SPACE=vsp;
    sb.solution.copy_nonactive(sb.rhs);
    double t2=TDatabase::TimeDB->THETA2;
    double t3=TDatabase::TimeDB->THETA3;
    double t4=TDatabase::TimeDB->THETA4;
    double tau=TDatabase::TimeDB->TIMESTEPLENGTH;
    // scale the right hand side with current step length
    sb.rhs.scaleActive(tau*t4);
    // add rhs from previous time step
    if(t3 !=0)
    {
      sb.rhs.addScaledActive(this->oldrhs_, tau*t3);
      // save the right hand side at current time step
      // for the next time step as old rhs
      oldrhs_.addScaledActive(sb.rhs, -1/(t3*tau));
      oldrhs_.scaleActive(-t3/t4);
      oldrhs_.copy_nonactive(sb.rhs);
    }
    // mass matrix times the old sol to current rhs
    sd.modifedMassMatrix_.apply_scaled_submatrix(oldsol_, sb.rhs, 2, 2, 1.0);
    // A's time the old solution to the rhs
    if(t2 !=0)
      sb.matrix.apply_scaled_submatrix(oldsol_, sb.rhs, 2, 2, t2*tau);
    // scale BT's with the time step length
    for(System_per_grid& s : this->Time_NSE2D::systems)
    {
      if(tau != oldtau)
      {
        double factor=TDatabase::TimeDB->THETA1*tau;
        if(oldtau !=0.)
        {
          factor /= oldtau;
          Output::print("change in tau", oldtau, "->", tau);
        }
        const std::vector<std::vector<size_t>> cell_positions = {{0,2}, {1,2}};
        s.matrix.scale_blocks(factor, cell_positions);
        if(TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT > 0)
        {
          const std::vector<std::vector<size_t>> cell_positions_t = {{2,0}, {2,1}};
          s.matrix.scale_blocks(factor, cell_positions_t);
        }
      }
    }// endfor
    oldtau = tau;
    sb.rhs.copy_nonactive(sb.solution);
  }
  else
  {
    this->Time_NSE2D::assemble_rhs();
  }
  Output::print<5>("assembled the system right hand side ");
}

void PR_Time_NSE2D::system_assemble()
{
  double t1=TDatabase::TimeDB->THETA1;
  double tau=TDatabase::TimeDB->TIMESTEPLENGTH;
  int i=0;
  if(TDatabase::ParamDB->DISCTYPE==RECONSTRUCTION)
  {
    for(System_per_grid& sb : this->systems)
    {
      // scale the A-blocks with tau*t1
      sb.matrix.scale_blocks_actives(t1*tau, {{0,0},{0,1},{1,0},{1,1}});
      // M00 + tau*t1* A00
      const FEMatrix& mfm00 =
         *this->systems_[i].modifedMassMatrix_.get_blocks().at(0).get();
      sb.matrix.add_matrix_actives(mfm00, 1.0, {{0,0}}, {false});
      // M01 + tau*t1* A01
      const FEMatrix& mfm01 =
         *this->systems_[i].modifedMassMatrix_.get_blocks().at(1).get();
      sb.matrix.add_matrix_actives(mfm01, 1.0, {{0,1}}, {false});

      // M10 + tau*t1* A10
      const FEMatrix& mfm10 =
         *this->systems_[i].modifedMassMatrix_.get_blocks().at(3).get();
      sb.matrix.add_matrix_actives(mfm10, 1.0, {{1,0}}, {false});

      // M11 + tau*t1* A11
      const FEMatrix& mfm11 =
         *this->systems_[i].modifedMassMatrix_.get_blocks().at(4).get();
      sb.matrix.add_matrix_actives(mfm11, 1.0, {{1,1}}, {false});
      i++;
    }
  }
  else
  {
    this->Time_NSE2D::assemble_system();
  }
}

void PR_Time_NSE2D::system_solve()
{
  if(TDatabase::ParamDB->DISCTYPE != RECONSTRUCTION)
  {
    Output::print("solving the system for CL FEM ");
    this->Time_NSE2D::solve();
    return;
  }
  Output::print("solving the system for PR FEM ");
  System_per_grid& sb=this->systems.front();
  if(this->solver_.is_using_multigrid())
    ErrThrow("multigrid solver is not tested yet");
  solver_.solve(sb.matrix, sb.rhs, sb.solution);
  // subtract the mass matrix and rescale to get back the A's
  this->descale();

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    sb.p.project_into_L20();

  // copy solution for next time step
  this->oldsol_ = sb.solution;
}

void PR_Time_NSE2D::descale()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double t1=TDatabase::TimeDB->THETA1;
  int i=0;
  for(System_per_grid& s : this->systems)
  {
    const FEMatrix& mfm00
      = *this->systems_[i].modifedMassMatrix_.get_blocks().at(0).get();
    const FEMatrix& mfm01
      = *this->systems_[i].modifedMassMatrix_.get_blocks().at(1).get();
    const FEMatrix& mfm10
      = *this->systems_[i].modifedMassMatrix_.get_blocks().at(3).get();
    const FEMatrix& mfm11
      = *this->systems_[i].modifedMassMatrix_.get_blocks().at(4).get();

    s.matrix.add_matrix_actives(mfm00, -1.0, {{0,0}}, {false});
    s.matrix.add_matrix_actives(mfm01, -1.0, {{0,1}}, {false});
    s.matrix.add_matrix_actives(mfm10, -1.0, {{1,0}}, {false});
    s.matrix.add_matrix_actives(mfm11, -1.0, {{1,1}}, {false});
    s.matrix.scale_blocks_actives(1./(t1*tau), {{0,0}, {0,1}, {1, 0}, {1, 1}});
    i++;
  }
}
