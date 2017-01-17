#include <PR_NSE2D.h>
#include <Database.h>
#include <Multigrid.h>
#include <LocalAssembling2D.h>
#include <Assemble2D.h>


ParameterDatabase get_default_PR_NSE2D_parameters()
{
  Output::print<5>("creating a default NSE2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NSE2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("PR_NSE2D parameter database");

  //NSE2D requires a nonlinear iteration, set up a nonlinit_database and merge
  db.merge(ParameterDatabase::default_nonlinit_database());

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  //stokes case - reduce no nonlin its TODO remove global database dependency
  if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 3)
  {
    if (TDatabase::ParamDB->PRESSURE_SEPARATION==1)
    {
       db["nonlinloop_maxit"] = 1;
    }
    else
    {
      db["nonlinloop_maxit"] = 1;
    }
  }

  return db;
}

PR_NSE2D::system_per_grid::system_per_grid(const Example_NSE2D &ex_,
                                           TCollection &coll, size_t order)
    : projSpace_(&coll, (char*)"pr", (char*)"pr",ex_.get_bc(2),
                 order, nullptr)
{

}

PR_NSE2D::PR_NSE2D(const TDomain & domain, const ParameterDatabase& param_db,
                   const Example_NSE2D _example, unsigned int reference_id)
  : NSE2D(domain, param_db, _example, reference_id),
    system_(), example_(_example),db_(get_default_PR_NSE2D_parameters()),
    solver_(param_db)
{
  db_.merge(param_db, false);
  // set projection space according to the velocity space
  this->set_projection_space();

  TCollection*coll = domain.GetCollection(It_Finest,0);
  size_t order=TDatabase::ParamDB->PROJECTION_SPACE;
  system_.emplace_back(example_, *coll, order);
}

void PR_NSE2D::set_projection_space() const
{
  switch(TDatabase::ParamDB->VELOCITY_SPACE)
  {
    case 2: // BDM2
      TDatabase::ParamDB->PROJECTION_SPACE = 1012;
      break;
    case 22:
      TDatabase::ParamDB->PROJECTION_SPACE = 1001;
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

void PR_NSE2D::print_info() const
{
  int n_dof_pr=this->get_projection_space().GetN_DegreesOfFreedom();
  int n_active=this->get_projection_space().GetActiveBound();
  Output::print("dof RT/BDM        : ", setw(10), n_dof_pr);
  Output::print("dof active RT/BDM : ", setw(10), n_active);

  if(TDatabase::ParamDB->DISCTYPE == RECONSTRUCTION)
      Output::print("Pressure Robust FEM");
  else
      Output::print("      Classical FEM");
}

void PR_NSE2D::all_assemble()
{
  // assembling of the matrices and rhs at once
  this->NSE2D::assemble();

  if(TDatabase::ParamDB->DISCTYPE == GALERKIN)
  {
    Output::print("Classical       FEM is considered");
    return;
  }
  else if(TDatabase::ParamDB->DISCTYPE==RECONSTRUCTION)
  {
    Output::print("Pressure Robust FEM is considered");
  }
  else
  {
    Output::print("DISCTYPE: ", TDatabase::ParamDB->DISCTYPE, " is used");
    return;
  }
  // system_per_grid& sd = this->system_.front();
  for(System_per_grid& sb : this->systems)
  {
    TFEFunction2D *fef[3] = { sb.u.GetComponent(0), sb.u.GetComponent(1),&sb.p };
    size_t vspace_ = TDatabase::ParamDB->VELOCITY_SPACE;
    TDatabase::ParamDB->VELOCITY_SPACE = TDatabase::ParamDB->PROJECTION_SPACE;
    // local assemble
    LocalAssembling2D larhs(RECONSTR_GALERKIN_Rhs, fef, example_.get_coeffs());

    size_t nfesp=2;
    const TFESpace2D* projsp=&this->get_projection_space();
    const TFESpace2D* velosp=&this->NSE2D::get_velocity_space();
    const TFESpace2D* pointer_to_space[2] = {velosp, projsp};
    //TODO: do we really need to assemble square matrices???
    size_t nsqmatass=2;
    std::vector<int> rspace={0,0,1};// row space for assembling (local?) matrices
    std::vector<int> cspace={1,1,1};// column space for assembling (local?) matricces

    // rhs
    size_t nrhsass=1;
    std::vector<int> r_rhsspac={1}; // row space for assembling rhs
    // variable for storing matrices and rhs
    size_t nsqmatstore=0;
    TSquareMatrix2D** sqmatstore=nullptr;
    size_t nrematstore=0;
    TMatrix2D** rematstore=nullptr;

    size_t nrhsstore=2;
    const TFESpace2D* pressp=&this->NSE2D::get_pressure_space();
    const TFESpace2D* pointer_to_rhsspace[3]={velosp,velosp,pressp};
    sb.rhs.reset();
    double *rhsstore[3]={sb.rhs.block(0), sb.rhs.block(1),sb.rhs.block(2)};

    // boundary conditions and boundary values
    BoundCondFunct2D *bc[3] ={ velosp->GetBoundCondition(),
                               velosp->GetBoundCondition(),
                               pressp->GetBoundCondition() };
    // boundary values
    std::array<BoundValueFunct2D*, 3> bv;
    bv[0] = example_.get_bd()[0];
    bv[1] = example_.get_bd()[1];
    bv[2] = example_.get_bd()[2];
    Assemble2D_VectFE(nfesp,pointer_to_space,
                      nsqmatass, rspace, cspace,nrhsass,r_rhsspac,
                      nsqmatstore, sqmatstore, nrematstore,rematstore,
                      nrhsstore, rhsstore,pointer_to_rhsspace,
                      larhs, nullptr, MatVectMult,
                      projection_matrices,bc,bv.data());
    TDatabase::ParamDB->VELOCITY_SPACE = vspace_;
  }
}

void PR_NSE2D::nonlinear_assemble()
{
  // standard Galerkin case
  if(TDatabase::ParamDB->DISCTYPE == GALERKIN)
  {
    this->NSE2D::assemble_nonlinear_term();
    return;
  }
  
  for(System_per_grid& s : this->systems)
  {
    std::vector<TFEFunction2D*> fefct;
    fefct.push_back(s.u.GetComponent(0));
    fefct.push_back(s.u.GetComponent(1));
    fefct.push_back(&s.p);
    // local assembling object
    LocalAssembling2D la(RECONSTR_NLGALERKIN, fefct.data(), example_.get_coeffs());
    // setting the space for sign in GetSignOfThisDOF();
    size_t vel_space = TDatabase::ParamDB->VELOCITY_SPACE;
    TDatabase::ParamDB->VELOCITY_SPACE = TDatabase::ParamDB->PROJECTION_SPACE;
    // prepare everything which is needed for assembling matrices 
    size_t nfesp = 2;
    const TFESpace2D* pro_sp  = &this->get_projection_space();
    const TFESpace2D* vel_sp = &s.velocity_space;
    // pointers to the fespace
    const TFESpace2D* pointer_to_space[2] = { vel_sp, pro_sp };
    
    size_t nsq_ass=2;
    size_t nre_ass=4;
    // row space for assembling matrices
    std::vector<int> rspace ={0, 0, 1, 1, 0, 0}; 
    // cols space for assembling matrices
    std::vector<int> cspace ={0, 0, 0, 0, 1, 1}; 
    // square matrices to be stored
    size_t nsq_stor = 4;
    TSquareMatrix2D* sq_mat_store[4]{nullptr};
    // rectangular matrices to be stored
    size_t nre_store = 0;
    TMatrix2D** re_mat_store=nullptr;
     // right hand side to be assembled
    size_t nrhs_ass = 0;
    // row space for the right hand side
    std::vector<int> rspace_rhs={};
    size_t nrhs_sto=0;
    double **rhs_stor{nullptr};
    const TFESpace2D ** pointer_to_rhs_space{nullptr};
    
    // matrices
    std::vector<std::shared_ptr<FEMatrix>>mat
                        = s.matrix.get_blocks_uniquely();
    sq_mat_store[0] = reinterpret_cast<TSquareMatrix2D*>(mat.at(0).get());
    sq_mat_store[1] = reinterpret_cast<TSquareMatrix2D*>(mat.at(1).get());
    sq_mat_store[2] = reinterpret_cast<TSquareMatrix2D*>(mat.at(3).get());
    sq_mat_store[3] = reinterpret_cast<TSquareMatrix2D*>(mat.at(4).get());
    // boundary conditions and boundary values
    const TFESpace2D* pre_sp = &s.pressure_space;
    BoundCondFunct2D *bc[3]={vel_sp->GetBoundCondition(), 
                             vel_sp->GetBoundCondition(),
                             pre_sp->GetBoundCondition()};
    std::array<BoundValueFunct2D*, 3> bv;
    bv[0]=example_.get_bd(0);
    bv[1]=example_.get_bd(1);
    bv[2]=example_.get_bd(2);
    
    Assemble2D_MixedFEM(nfesp, pointer_to_space,
                 nsq_ass, nre_ass, rspace, cspace, nrhs_ass, rspace_rhs,
                 nsq_stor, sq_mat_store, nre_store, re_mat_store,
                 nrhs_sto, rhs_stor,pointer_to_rhs_space, la,
                 nonlinear_term_reconstruct, nullptr,
                 ProjectionMatricesNSE2D, bc, bv.data());
    TDatabase::ParamDB->VELOCITY_SPACE = vel_space;
  }
}
