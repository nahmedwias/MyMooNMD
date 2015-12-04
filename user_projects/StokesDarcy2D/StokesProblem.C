#include <algorithm>
#include <StokesProblem.h>
#include <solver.h>
#include <preconditioner.h>
#include <FEDatabase2D.h>

/** ************************************************************************ */
StokesProblem::StokesProblem(const TDomain & domain, 
                             std::shared_ptr<Example_NSE2D> ex, 
                             InterfaceCondition t, 
                             std::vector<const TInnerInterfaceJoint*>& in,
                             int c)
 : NSE2D(domain, *ex, 2),
   typeOf_bci(t), DrhsNSE(this->get_rhs()), darcyToStokes(nullptr),
   // actually the old solution is only needed if 
   // TDatabase::ParamDB->StoDa_solutionStrategy != 0, but we do this anyway
   solution_old(this->get_solution()), 
   fe_u_solution_old(&this->get_velocity_space(), (char*)"Stokes_u_old",
                     (char*)"Stokes_u_old", solution_old.block(0), 
                     solution_old.length(0), 2),
   fe_p_solution_old(&this->get_pressure_space(), (char*)"Stokes_p_old", 
                     (char*)"Stokes_p_old", solution_old.block(2), 
                     solution_old.length(2)),
   // actually the direct solution is only needed if 
   // TDatabase::ParamDB->StoDa_solutionStrategy > 0, but we do this anyway
   solution_direct(this->get_solution()),
   fe_u_solution_direct(&this->get_velocity_space(), (char*)"Stokes_u_direct",
                        (char*)"Stokes_u_direct", solution_direct.block(0), 
                        solution_direct.length(0), 2),
   fe_p_solution_direct(&this->get_pressure_space(), (char*)"Stokes_p_old", 
                        (char*)"Stokes_p_old", solution_direct.block(2), 
                        solution_direct.length(2)),
   interface(in), icond(c), mat_aux(nullptr), etaToBd(nullptr), 
   map_sol2eta(nullptr), map_sol2eta_rhs(), periodic_dofs(), eta_hom(nullptr),
   global_DOF_interface(), normals(), assemble_on_input(false), 
   assemble_on_return(true)
{
  //OutPut("StokesProblem constructor\n");
  
  // the matrix blocks in this->NSE2D::matrix all have the same structure (where
  // possible), since we want to change those structures individually, we copy
  // the structure often enough so that each matrix block has its own structure
  const TStructure & sq = 
    this->get_matrix().get_A_block(0)->GetStructure();
  this->get_matrix().get_A_block(1)->SetStructure(std::make_shared<TStructure>(sq));
  this->get_matrix().get_A_block(2)->SetStructure(std::make_shared<TStructure>(sq));
  this->get_matrix().get_A_block(3)->SetStructure(std::make_shared<TStructure>(sq));
  const TStructure & rs1 = this->get_matrix().get_BT_block(0)->GetStructure();
  this->get_matrix().get_BT_block(1)->SetStructure(std::make_shared<TStructure>(rs1));
  const TStructure & rs2 = this->get_matrix().get_B_block(0)->GetStructure();
  this->get_matrix().get_B_block(1)->SetStructure(std::make_shared<TStructure>(rs2));
  
  // copy values from right hand side. 
  DrhsNSE = this->get_rhs(); // From now on DrhsNSE is not changed any more
  this->get_rhs().reset();
  
  
  // values for coupled problem or fixed point formulation
  if(TDatabase::ParamDB->StoDa_problemType == 1
     && TDatabase::ParamDB->StoDa_updatingStrategy == 4)
  {
    // D-RR
    assemble_on_input = true;
    assemble_on_return = false;
  }
  
  if(TDatabase::ParamDB->StoDa_solutionStrategy == 3)
  {
    // if Stecklov Poincare is solved
    bool stokes_first = TDatabase::ParamDB->StoDa_StokesFirst == 1;
    switch(typeOf_bci)
    {
      case Neumann:
        assemble_on_input = stokes_first ? true : false;
        assemble_on_return = stokes_first ? true : false;
        break;
      case Dirichlet:
        assemble_on_input = true;
        assemble_on_return = true;
        // note that here it should be stokes_first == false, because in the
        // other case we don't need a Dirichlet Stokes problem
        if(TDatabase::ParamDB->EXAMPLE == 0) // where else could this be set?
          TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
        break;
      case Robin:
        ErrThrow("Robin interface conditions for Stecklov-Poincare not ",
                 "possible yet");
        break;
      case weakRobin:
        ErrThrow("weak Robin interface conditions for Stecklov-Poincare not ",
                 "possible yet");
        break;
      case DirichletSTAB:
        ErrThrow("DirichletSTAB interface conditions for Stecklov-Poincare ",
                 "not possible yet");
        break;
      default:
        ErrThrow("unsupported boundary conditions on interface");
        break;
    }
  }
}

/** ************************************************************************ */
StokesProblem::~StokesProblem()
{
  //OutPut("StokesProblem destructor\n");
}

/** ************************************************************************ */
void StokesProblem::mat_NSE_iIntegrals(BlockMatrixNSE2D* m)
{
  if(m == nullptr)
    m = &this->get_matrix();
  Output::print<1>(" Assemble interface integrals into matrix (NSE)");
  
  const int ID = this->get_velocity_space().GetCollection()->GetCell(0)
                 ->GetReference_ID();
  const int vActiveBound = this->get_velocity_space().GetActiveBound();
  
  
  bool velocity_pressure_coupling = true;
  if(typeOf_bci == Neumann || typeOf_bci == Robin
     || (typeOf_bci == Dirichlet && TDatabase::ParamDB->StoDa_weakGamma <= 0.0))
    velocity_pressure_coupling = false;
  bool pressure_pressure_coupling = false;
  if(typeOf_bci == weakRobin)
    pressure_pressure_coupling = true;
  
  // everything needed for the local assembling routine
  local_edge_assembling l;
  
  for(unsigned int j = 0; j < interface.size(); j++)
  {
    TBaseCell *cell = interface[j]->GetNeighbour(0);
    if(cell->GetReference_ID() != ID)
      cell = interface[j]->GetNeighbour(1);
    // index of this edge in this cell (0,1, or 2 for triangles)
    const int iEdge = interface[j]->GetIndexInNeighbor(cell);
    
    FE2D vFEId = this->get_velocity_space().GetFE2D(0, cell);
    FE2D pFEId = this->get_pressure_space().GetFE2D(0, cell);
    TBaseFunct2D *vBaseFunct = TFEDatabase2D::GetBaseFunct2DFromFE2D(vFEId);
    TBaseFunct2D *pBaseFunct = TFEDatabase2D::GetBaseFunct2DFromFE2D(pFEId);
    const int fe_degree_v = vBaseFunct->GetPolynomialDegree();
    const int fe_degree_p = pBaseFunct->GetPolynomialDegree();
    const int degree = fe_degree_v + MAX(fe_degree_v,fe_degree_p);
    
    QuadFormula1D LineQF = TFEDatabase2D::GetQFLineFromDegree(degree);
    int N_LinePoints; // number of quadrature points on this edge
    double *LineWeights, *zeta; // quadrature weights and points
    {
      TQuadFormula1D *qf = TFEDatabase2D::GetQuadFormula1D(LineQF);
      qf->GetFormulaData(N_LinePoints, LineWeights, zeta);
    }
    vBaseFunct->MakeRefElementData(LineQF);
    pBaseFunct->MakeRefElementData(LineQF);
    
    const int N_vBaseFunct = vBaseFunct->GetDimension();
    const int N_pBaseFunct = pBaseFunct->GetDimension();
    
    RefTrans2D RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(vFEId);
    TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
    
    int const * const vDOF = 
      this->get_velocity_space().GetGlobalDOF(cell->GetCellIndex());
    int const * const pDOF = 
      this->get_pressure_space().GetGlobalDOF(cell->GetCellIndex());
    
    double **uref = TFEDatabase2D::GetJointValues2D(vBaseFunct->GetID(), LineQF,
                                                    iEdge);
    double **uxiref = TFEDatabase2D::GetJointDerivatives2D(vBaseFunct->GetID(),
                                                           LineQF, iEdge, D10);
    double **uetaref = TFEDatabase2D::GetJointDerivatives2D(vBaseFunct->GetID(),
                                                            LineQF, iEdge, D01);
    double **pref = TFEDatabase2D::GetJointValues2D(pBaseFunct->GetID(), LineQF,
                                                    iEdge);
    double **pxiref = TFEDatabase2D::GetJointDerivatives2D(pBaseFunct->GetID(),
                                                           LineQF, iEdge, D10);
    double **petaref = TFEDatabase2D::GetJointDerivatives2D(pBaseFunct->GetID(),
                                                            LineQF, iEdge, D01);
    
    double *uorig = new double[N_vBaseFunct];
    double *uxorig = new double[N_vBaseFunct];
    double *uyorig = new double[N_vBaseFunct];
    double *porig = new double[N_pBaseFunct];
    double *pxorig = new double[N_pBaseFunct];
    double *pyorig = new double[N_pBaseFunct];
    
    // compute length of the edge
    const double hE = interface[j]->GetLength();
    // tangential/normal vectors to this boundary (normalized)
    double nx, ny;
    getNormal(cell, interface[j], nx, ny);
    const double tx = -ny, ty = nx;
    
    // these maps correspond to local matrices which will be added to the 
    // respective global matrix blocks (A11, A12, A21, A22, B1, B2, B1T, B2T)
    std::map<int, std::map<int, double> > b_x, b_y, b_xt, b_yt;
    
    l.nt.set(nx, ny, tx, ty);
    l.hE = hE;
    
    // loop of quadrature points
    for(int k = 0; k < N_LinePoints; k++)
    {
      TFEDatabase2D::GetOrigValues(RefTrans, zeta[k], vBaseFunct, iEdge,
                                   uref[k], uxiref[k], uetaref[k], uorig,
                                   uxorig, uyorig);
      TFEDatabase2D::GetOrigValues(RefTrans, zeta[k], pBaseFunct, iEdge,
                                   pref[k], pxiref[k], petaref[k], porig,
                                   pxorig, pyorig);
      
      l.qw = LineWeights[k] * (hE / 2);
      for(int lt = 0; lt < N_vBaseFunct; lt++) // local test index
      {
        const int Test_vDOF = vDOF[lt];
        if(Test_vDOF >= vActiveBound)
          continue;
        const double v = uorig[lt]; // value of velocity test function
        const double vx = uxorig[lt]; // x-derivative of velocity test function
        const double vy = uyorig[lt]; // y-derivative of velocity test function
        
        l.testDOF = Test_vDOF;
        l.at.setTest(v, vx, vy);
        for(int l2 = 0; l2 < N_vBaseFunct; l2++)
        {
          l.at.setAnsatz(uorig[l2], uxorig[l2], uyorig[l2]);
          l.ansatzDOF = vDOF[l2];
          
          localAssembleMat_velocity_velocity(l);
        }
        
        // velocity -- pressure entries
        if(velocity_pressure_coupling)
        {
          for(int l2 = 0; l2 < N_pBaseFunct; l2++)
          {
            l.at.setAnsatz(porig[l2], pxorig[l2], pyorig[l2]);
            l.ansatzDOF = pDOF[l2];
            
            localAssembleMat_velocity_pressure(l);
          }
        }
      } //for(int l=0;l<N_vBaseFunct;l++)
      
      // pressure -- pressure coupling
      if(pressure_pressure_coupling)
      {
        for(int lt = 0; lt < N_pBaseFunct; lt++)
        {
          l.at.setTest(porig[lt], pxorig[lt], pyorig[lt]);
          l.testDOF = pDOF[lt];
          
          for(int l2 = 0; l2 < N_pBaseFunct; l2++)
          {
            l.at.setAnsatz(porig[l2], pxorig[l2], pyorig[l2]);
            l.ansatzDOF = pDOF[l2];
            
            localAssembleMat_pressure_pressure(l);
          }
        }
      }
    } //for(int k=0;k<N_LinePoints;k++)
    // add the computed values of this row to the matrix.  
    
    delete[] uorig;
    delete[] uxorig;
    delete[] uyorig;
    delete[] porig;
    delete[] pxorig;
    delete[] pyorig;
    
  } //for(int j=0; j<N_InnerInterfaceJoints; j++)
  
  m->get_A_block(0)->add(l.m.m[0]);
  m->get_A_block(1)->add(l.m.m[1]);
  m->get_A_block(2)->add(l.m.m[2]);
  m->get_A_block(3)->add(l.m.m[3]);
  m->get_C_block()->add(l.m.m[4]);
  m->get_BT_block(0)->add(l.m.m[5]);
  m->get_BT_block(1)->add(l.m.m[6]);
  m->get_B_block(0)->add(l.m.m[7]);
  m->get_B_block(1)->add(l.m.m[8]);
  
  //reorder_rows(m);
}

/** ************************************************************************ */
void StokesProblem::localAssembleMat_velocity_velocity(local_edge_assembling &l)
{
  const double nu = 1.0 / TDatabase::ParamDB->RE_NR;
  const double nx = l.nt.nx, ny = l.nt.ny;
  const double tx = l.nt.tx, ty = l.nt.ty;
  const double u = l.at.u, ux = l.at.ux, uy = l.at.uy;
  const double v = l.at.v, vx = l.at.vx, vy = l.at.vy;
  const int tDOF = l.testDOF;
  const int aDOF = l.ansatzDOF;
  const double qw = l.qw;
  
  // terms in tangential direction
  switch(getIcond())
  {
    case 0:
    {
      // from Beavers-Joseph-Saffman condition
      const double alpha = TDatabase::ParamDB->StoDa_alpha;
      // normal component is either implemented weakly, or this is not a
      // Dirichlet problem
      l.m.m[0][tDOF][aDOF] += alpha * (v * tx) * (u * tx) * qw;
      l.m.m[1][tDOF][aDOF] += alpha * (v * tx) * (u * ty) * qw;
      l.m.m[2][tDOF][aDOF] += alpha * (v * ty) * (u * tx) * qw;
      l.m.m[3][tDOF][aDOF] += alpha * (v * ty) * (u * ty) * qw;
      break;
    }
    case 1:
    { 
      // interface condition u.t = 0  (always implemented weakly)
      // normal component is either implemented weakly, or this is not a
      // Dirichlet problem
      
      // -(t.T(u,p).n,v.t)  (this is left after integration by parts)
      l.m.m[0][tDOF][aDOF] -= (2 * nx * tx * ux + (ny * tx + nx * ty) * uy) * nu
                              * v * tx * qw;
      l.m.m[1][tDOF][aDOF] -= ((ny * tx + nx * ty) * ux + 2 * ny * ty * uy) * nu
                              * v * tx * qw;
      l.m.m[2][tDOF][aDOF] -= (2 * nx * tx * ux + (ny * tx + nx * ty) * uy) * nu
                              * v * ty * qw;
      l.m.m[3][tDOF][aDOF] -= ((ny * tx + nx * ty) * ux + 2 * ny * ty * uy) * nu
                              * v * ty * qw;
      // enforcing weak boundary condition u.t=0 on interface
      const double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma) / l.hE;
      l.m.m[0][tDOF][aDOF] += gamma * (v * tx) * (u * tx) * qw;
      l.m.m[1][tDOF][aDOF] += gamma * (v * tx) * (u * ty) * qw;
      l.m.m[2][tDOF][aDOF] += gamma * (v * ty) * (u * tx) * qw;
      l.m.m[3][tDOF][aDOF] += gamma * (v * ty) * (u * ty) * qw;
      break;
    }
    default:
      ErrThrow("unsupported interface condition. Choose either the ",
               "Beavers-Joseph-Saffman condition (0) or u.t=0 (1)");
      break;
  }
  
  // terms in normal direction
  switch(typeOf_bci)
  {
    case Neumann:
      break;
    case Robin:
    {
      // (u.n,v.n)
      const double gamma_f = TDatabase::ParamDB->StoDa_gamma_f * qw;
      l.m.m[0][tDOF][aDOF] += gamma_f * (v * nx) * (u * nx);
      l.m.m[1][tDOF][aDOF] += gamma_f * (v * nx) * (u * ny);
      l.m.m[2][tDOF][aDOF] += gamma_f * (v * ny) * (u * nx);
      l.m.m[3][tDOF][aDOF] += gamma_f * (v * ny) * (u * ny);
      break;
    }
    case weakRobin:
    {
      // -(2 nu n.DD(u).n,v.n) 
      const double val = -2 * nu * (nx * ux + ny * uy) * v * qw;
      l.m.m[0][tDOF][aDOF] += val * nx * nx;
      l.m.m[1][tDOF][aDOF] += val * ny * nx;
      l.m.m[2][tDOF][aDOF] += val * nx * ny;
      l.m.m[3][tDOF][aDOF] += val * ny * ny;
      // (gamma_f(n.T(u,p).n) + u_n - eta , v.n - gamma h_E n.T(v,q).n)
      //  *1/(gamma_f+gamma h_E)
      const double gamma_f = TDatabase::ParamDB->StoDa_gamma_f;
      const double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma);
      const double tilde_lambda = 1 / (gamma_f + gamma * l.hE);
      const double ansatzfct = (2 * nu * (nx * ux + ny * uy) + gamma_f * u)
                               * tilde_lambda * qw;
      const double testfct = v - gamma * l.hE * 2 * nu * (nx * vx + ny * vy);
      l.m.m[0][tDOF][aDOF] += ansatzfct * nx * testfct * nx;
      l.m.m[1][tDOF][aDOF] += ansatzfct * ny * testfct * nx;
      l.m.m[2][tDOF][aDOF] += ansatzfct * nx * testfct * ny;
      l.m.m[3][tDOF][aDOF] += ansatzfct * ny * testfct * ny;
      break;
    }
    case DirichletSTAB:
    {
      // (u.n,v.n) *(1/(h gamma_f nu))
      double gamma_f = TDatabase::ParamDB->StoDa_gamma_f
                       / TDatabase::ParamDB->RE_NR;
      gamma_f *= qw / l.hE;
      l.m.m[0][tDOF][aDOF] += gamma_f * (v * nx) * (u * nx);
      l.m.m[1][tDOF][aDOF] += gamma_f * (v * nx) * (u * ny);
      l.m.m[2][tDOF][aDOF] += gamma_f * (v * ny) * (u * nx); // = a_xy
      l.m.m[3][tDOF][aDOF] += gamma_f * (v * ny) * (u * ny);
      
      // -(n.T(u,p).n,v.n)*(gamma_f/(gamma_f+h))
      double val = -2 * nu * (nx * ux + ny * uy) * v;
      val *= TDatabase::ParamDB->StoDa_gamma_f
             / (TDatabase::ParamDB->StoDa_gamma_f + l.hE)
             * qw;
      l.m.m[0][tDOF][aDOF] += val * nx * nx;
      l.m.m[1][tDOF][aDOF] += val * ny * nx;
      l.m.m[2][tDOF][aDOF] += val * nx * ny;
      l.m.m[3][tDOF][aDOF] += val * ny * ny;
      break;
    }
    case Dirichlet:
    {
      if(TDatabase::ParamDB->StoDa_weakGamma <= 0.0)
      {
        // strong normal condition
        // nothing more to be done here
      }
      else
      {
        // weak normal condition
        // -(n.T(u,p).n,v.n)
        const double val = -2 * nu * (nx * ux + ny * uy) * v * qw;
        l.m.m[0][tDOF][aDOF] += val * nx * nx;
        l.m.m[1][tDOF][aDOF] += val * ny * nx;
        l.m.m[2][tDOF][aDOF] += val * nx * ny;
        l.m.m[3][tDOF][aDOF] += val * ny * ny;
        
        // (u.n,v.n) *(1/h)
        double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma);
        
        gamma *= qw / l.hE;
        l.m.m[0][tDOF][aDOF] += gamma * (v * nx) * (u * nx);
        l.m.m[1][tDOF][aDOF] += gamma * (v * nx) * (u * ny);
        l.m.m[2][tDOF][aDOF] += gamma * (v * ny) * (u * nx);
        l.m.m[3][tDOF][aDOF] += gamma * (v * ny) * (u * ny);
      }
      break;
    }
    default:
      ErrThrow("unsupported type of Problem. Choose either ",
               "Neumann_Neumann (0), Robin_Robin (1), weak Robin_Robin(2), ",
               "Dirichlet_DirichletSTAB(3) or Dirichlet_Dirichlet(4)");
      break;
  }
}

/** ************************************************************************ */
void StokesProblem::localAssembleMat_velocity_pressure(local_edge_assembling &l)
{
  const double nu = 1.0 / TDatabase::ParamDB->RE_NR;
  const double nx = l.nt.nx, ny = l.nt.ny;
  //const double tx = l.nt.tx, ty = l.nt.ty;
  const double p = l.at.u; //, px = l.at.ux, py = l.at.uy; // Stokes pressure
  const double v = l.at.v, vx = l.at.vx, vy = l.at.vy; // Stokes velocity
  const int tDOF = l.testDOF;
  const int aDOF = l.ansatzDOF;
  const double qw = l.qw;
  
  switch(typeOf_bci)
  {
    case Neumann:
    case Robin:
      break;
    case weakRobin:
    {
      // (gamma_f(n.T(u,p).n) + u_n , v.n - gamma h_E n.T(v,q).n)
      //  *1/(gamma_f+gamma h_E)
      const double gamma_f = TDatabase::ParamDB->StoDa_gamma_f;
      const double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma);
      
      // Stokes velocity test function, Stokes pressure ansatz function
      const double tilde_lambda = 1 / (gamma_f + gamma * l.hE);
      double ansatzfct = -p * tilde_lambda * qw;
      double testfct = v - gamma * l.hE * 2 * nu * (nx * vx + ny * vy);
      l.m.m[5][tDOF][aDOF] += ansatzfct * testfct * nx;
      l.m.m[6][tDOF][aDOF] += ansatzfct * testfct * ny;
      // Stokes pressure test function, Stokes velocity ansatz function
      ansatzfct = (2 * nu * (nx * vx + ny * vy) + gamma_f * v) * qw;
      testfct = tilde_lambda * gamma * l.hE * p;
      l.m.m[7][aDOF][tDOF] += ansatzfct * testfct * nx;
      l.m.m[8][aDOF][tDOF] += ansatzfct * testfct * ny;
      
      break;
    }
    case DirichletSTAB:
    {
      // -(n.T(u,p).n+eta,v.n)*(gamma_f/(gamma_f+h))
      double val = p * v * qw;
      val *= TDatabase::ParamDB->StoDa_gamma_f
             / (TDatabase::ParamDB->StoDa_gamma_f + l.hE);
      l.m.m[5][tDOF][aDOF] += val * nx;
      l.m.m[6][tDOF][aDOF] += val * ny;
      break;
    }
    case Dirichlet:
    {
      if(TDatabase::ParamDB->StoDa_weakGamma <= 0.0)
      {
        // strong normal condition
        // nothing more to be done here
      }
      else
      {
        // weak normal condition
        // -(n.T(u,p).n,v.n)
        const double val = p * v * qw;
        l.m.m[5][tDOF][aDOF] += val * nx;
        l.m.m[6][tDOF][aDOF] += val * ny;
      }
      break;
    }
    default:
      ErrThrow("unsupported type of Problem. Choose either ",
               "Neumann_Neumann (0), Robin_Robin (1), weak Robin_Robin(2), ",
               "Dirichlet_DirichletSTAB(3) or Dirichlet_Dirichlet(4)");
  }
}

/** ************************************************************************ */
void StokesProblem::localAssembleMat_pressure_pressure(local_edge_assembling &l)
{
  //const double nx = l.nt.nx, ny = l.nt.ny;
  //const double tx = l.nt.tx, ty = l.nt.ty;
  const double p = l.at.u; //, px = l.at.ux, py = l.at.uy;
  const double q = l.at.v; //, qx = l.at.vx, qy = l.at.vy;
  const int tDOF = l.testDOF;
  const int aDOF = l.ansatzDOF;
  const double qw = l.qw;
  switch(typeOf_bci)
  {
    case Neumann:
    case Robin:
    case DirichletSTAB:
    case Dirichlet:
      // nothing needs to be done, no pressure-pressure coupling
      // we should never be here, if 'pressure_pressure_coupling' is 
      // set correctly in 'mat_NSE_iIntegrals'
      break;
    case weakRobin:
    {
      // (gamma_f(n.T(u,p).n) + u_n , v.n - gamma h_E n.T(v,q).n)
      //  *1/(gamma_f+gamma h_E)
      const double gamma_f = TDatabase::ParamDB->StoDa_gamma_f;
      const double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma);
      
      const double tilde_lambda = 1 / (gamma_f + gamma * l.hE);
      l.m.m[4][tDOF][aDOF] += -p * gamma * l.hE * q * tilde_lambda * qw;
      break;
    }
    default:
      ErrThrow("unsupported type of Problem. Choose either ",
               "Neumann_Neumann (0), Robin_Robin (1), weak Robin_Robin(2), ",
               "Dirichlet_DirichletSTAB(3) or Dirichlet_Dirichlet(4)");
      break;
  }
}

/** ************************************************************************ */
void StokesProblem::localAssembleEtaToBd_velocity(local_edge_assembling &l)
{
  const double nu = 1.0 / TDatabase::ParamDB->RE_NR;
  const int N_SveloDOF = this->get_velocity_space().GetN_DegreesOfFreedom();
  const double nx = l.nt.nx, ny = l.nt.ny;
  //const double tx = l.nt.tx, ty = l.nt.ty;
  const double xi = l.at.u;
  const double v = l.at.v, vx = l.at.vx, vy = l.at.vy;
  //const bool nt = (POW(normal(l.testDOF).nx,2) >= POW(normal(l.testDOF).ny,2));
  //const int tDOF = l.testDOF + (nt ? 0 : N_SveloDOF);
  const int tDOF = l.testDOF;
  const int aDOF = l.ansatzDOF;
  
  switch(typeOf_bci)
  {
    case Neumann:
    case Robin:
    {
      // (eta,v.n)
      // this will be put into the row normal direction, so 'v' here is already
      // multiplied by the normal
      l.m.m[0][tDOF][aDOF] += xi * v * l.qw;  
      //l.m.m[0][tDOF][aDOF] += xi * v * nx * l.qw;
      //l.m.m[0][tDOF + N_SveloDOF][aDOF] += xi * v * ny * l.qw;
      //l.m.m[0][tDOF + 2*N_SveloDOF][aDOF] += 0.0; // pressure component
      if(kink(tDOF%N_SveloDOF))
      {
        // we need an entry in the 'tangential' direction, too. Maybe just
        // xi*v*l.qw as well
        
        // this is a generalized tangential, defined as a orthogonal vector
        // to the generalized normal
        double tx = normal(tDOF%N_SveloDOF).tx;
        double ty = normal(tDOF%N_SveloDOF).ty;
        // tmp is not zero, because (tx,ty) is not the tangential to (nx,ny)
        double tmp = nx*tx + ny*ty;
        l.m.m[0][(tDOF+N_SveloDOF)%(2*N_SveloDOF)][aDOF] += xi * v * l.qw * tmp;
      }
      break;
    }
    case weakRobin:
    {
      ErrThrow("localAssembleEtaToBd_velocity: check with nt");
      // (eta , v.n - gamma h_E n.T(v,q).n) *1/(gamma_f+gamma h_E)
      const double gamma_f = TDatabase::ParamDB->StoDa_gamma_f;
      const double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma);
      const double tilde_lambda = 1 / (gamma_f + gamma * l.hE);
      double val2 = xi * (v - gamma * l.hE * 2 * nu * (vx * nx + vy * ny));
      val2 *= tilde_lambda;
      l.m.m[0][tDOF][aDOF] += val2 * nx * l.qw;
      l.m.m[0][tDOF + N_SveloDOF][aDOF] += val2 * ny * l.qw;
      break;
    }
    case DirichletSTAB:
    {
      ErrThrow("Dirichlet_DirichletSTAB does not work iteratively yet");
      break;
    }
    case Dirichlet:
    {
      ErrThrow("localAssembleEtaToBd_velocity: check with nt");
      double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma);
      // (eta,v.n) * (1/h)
      l.m.m[0][tDOF][aDOF] += xi * v * nx * gamma * l.qw / l.hE;
      l.m.m[0][tDOF + N_SveloDOF][aDOF] += xi * v * ny * gamma * l.qw / l.hE;
      break;
    }
    default:
      ErrThrow("unsupported type of Problem. Choose either ",
               "Neumann_Neumann (0), Robin_Robin (1), weak Robin_Robin(2), ",
               "Dirichlet_DirichletSTAB(3) or Dirichlet_Dirichlet(4)");
      break;
  }
}

/** ************************************************************************ */
void StokesProblem::localAssembleEtaToBd_pressure(local_edge_assembling &l)
{
  //const int N_SveloDOF = this->get_velocity_space().GetN_DegreesOfFreedom();
  //const int tDOF = l.testDOF;
  //const int aDOF = l.ansatzDOF;
  //l.m.m[0][tDOF+2*N_SveloDOF][aDOF] += 0.0;
}

/** ************************************************************************ */
void StokesProblem::reorder_rows(BlockMatrixNSE2D* m)
{
  if(m == NULL)
    m = &this->get_matrix();
  
  const int ID = 
    this->get_velocity_space().GetCollection()->GetCell(0)->GetReference_ID();
  const int vActiveBound = this->get_velocity_space().GetActiveBound();
  
  std::map<int, std::map<int, double> > new_rows_a11;
  std::map<int, std::map<int, double> > new_rows_a12;
  std::map<int, std::map<int, double> > new_rows_a21;
  std::map<int, std::map<int, double> > new_rows_a22;
  std::map<int, std::map<int, double> > new_rows_b1t;
  std::map<int, std::map<int, double> > new_rows_b2t;
  
  int n_v_dof = this->get_velocity_space().GetN_DegreesOfFreedom();
  const bool strong_normal = TDatabase::ParamDB->StoDa_weakGamma <= 0.0
                             && typeOf_bci == Dirichlet;
  
  for(unsigned int j = 0; j < interface.size(); j++)
  {
    TBaseCell *cell = interface[j]->GetNeighbour(0);
    if(cell->GetReference_ID() != ID)
      cell = interface[j]->GetNeighbour(1);
    // index of this edge in this cell (0,1, or 2 for triangles)
    
    const int cell_index = cell->GetCellIndex();
    int const * const vDOF = get_velocity_space().GetGlobalDOF(cell_index);
    FE2D vFEId = this->get_velocity_space().GetFE2D(0, cell);
    TBaseFunct2D *vBaseFunct = TFEDatabase2D::GetBaseFunct2DFromFE2D(vFEId);
    const int N_vBaseFunct = vBaseFunct->GetDimension();
    
    for(int lvdof = 0; lvdof < N_vBaseFunct; lvdof++) // local velocity dof
    {
      const int this_vdof = vDOF[lvdof];
      if(global_DOF_interface[this_vdof][0].second == -1)
        continue; // not an interface dof
      if(this_vdof >= vActiveBound)
        continue; // not an active dof (dirichlet boundary)
        
      if(new_rows_a11.count(this_vdof) == 1
         || new_rows_a21.count(this_vdof) == 1)
        continue; // this dof has been handled in neighboring interface edge
      
      const double nx = normal(this_vdof).nx;
      const double ny = normal(this_vdof).ny;
      const double tx = normal(this_vdof).tx;
      const double ty = normal(this_vdof).ty;
      const bool nt = nx*nx - ny*ny > -1e-12; // (nx * nx >= ny * ny);
      
      if(strong_normal && getIcond() == 1) // u.t=0 and u.n=eta
      {
        new_rows_a11[this_vdof][this_vdof] = nt ? nx : tx;
        new_rows_a12[this_vdof][this_vdof] = nt ? ny : ty;
        new_rows_b1t[this_vdof]; // set entire row to zero
        new_rows_a21[this_vdof][this_vdof] = nt ? tx : nx;
        new_rows_a22[this_vdof][this_vdof] = nt ? ty : ny;
        new_rows_b2t[this_vdof]; // set entire row to zero
      }
      else if(TDatabase::ParamDB->StoDa_weakGamma <= 0.0 && getIcond() == 1)
      {
        // set tangential component
        if(nt)
        {
          new_rows_a21[this_vdof][this_vdof] = tx;
          new_rows_a22[this_vdof][this_vdof] = ty;
          new_rows_b2t[this_vdof]; // set entire row to zero
        }
        else
        {
          new_rows_a11[this_vdof][this_vdof] = tx;
          new_rows_a12[this_vdof][this_vdof] = ty;
          new_rows_b1t[this_vdof]; // set entire row to zero
        }
      }
      else if(!strong_normal || !nt)
      {
        // first component ( if nt==true, this is the normal component)
        copy_row(this_vdof, m->get_A_block(0), new_rows_a11, nt ? nx : tx);
        copy_row(this_vdof, m->get_A_block(2), new_rows_a11, nt ? ny : ty);
        copy_row(this_vdof, m->get_A_block(1), new_rows_a12, nt ? nx : tx);
        copy_row(this_vdof, m->get_A_block(3), new_rows_a12, nt ? ny : ty);
        copy_row(this_vdof, m->get_BT_block(0), new_rows_b1t, nt ? nx : tx);
        copy_row(this_vdof, m->get_BT_block(1), new_rows_b1t, nt ? ny : ty);
      }
      else
      {
        new_rows_a11[this_vdof][this_vdof] = nx;
        new_rows_a12[this_vdof][this_vdof] = ny;
        new_rows_b1t[this_vdof]; // set entire row to zero
        // this is done at the end of Assemble_Neumann_map_to_interface
        //(*rhs)[row] = 0; 
      }
      if(!strong_normal || nt)
      {
        // second component ( if nt==true, this is the tangential component)
        copy_row(this_vdof, m->get_A_block(0), new_rows_a21, nt ? tx : nx);
        copy_row(this_vdof, m->get_A_block(2), new_rows_a21, nt ? ty : ny);
        copy_row(this_vdof, m->get_A_block(1), new_rows_a22, nt ? tx : nx);
        copy_row(this_vdof, m->get_A_block(3), new_rows_a22, nt ? ty : ny);
        copy_row(this_vdof, m->get_BT_block(0), new_rows_b2t, nt ? tx : nx);
        copy_row(this_vdof, m->get_BT_block(1), new_rows_b2t, nt ? ty : ny);
      }
      else
      {
        new_rows_a21[this_vdof][this_vdof] = nx;
        new_rows_a22[this_vdof][this_vdof] = ny;
        new_rows_b2t[this_vdof]; // set entire row to zero
        // this is done at the end of Assemble_Neumann_map_to_interface
        //(*rhs[0])[row + n_dof] = 0;
      }
      double new_rhs1 = this->get_rhs().at(this_vdof) * (nt ? nx : tx)
                     + this->get_rhs().at(this_vdof + n_v_dof) * (nt ? ny : ty);
      double new_rhs2 = this->get_rhs().at(this_vdof) * (nt ? tx : nx)
                     + this->get_rhs().at(this_vdof + n_v_dof) * (nt ? ty : ny);
      this->get_rhs().at(this_vdof) = new_rhs1;
      this->get_rhs().at(this_vdof + n_v_dof) = new_rhs2;
      
      /*else
      { // kink
        // we have two normals here, one is (nx,ny) from the current (j-th) 
        // interface edge. The other is from the neighboring interface edge 
      
        // enforce two conditions, one for each of the two normal components
        
        // the two normals n1, n2
        n_t_vector n1 = normal(this_vdof, 0);
        n_t_vector n2 = normal(this_vdof, 1);
        n_t_vector n = normal(this_vdof); // take modified normal
        if(!strong_normal)
        {
          copy_row(this_vdof, m->squareBlock(0), new_rows_a11, n1.nx);
          copy_row(this_vdof, m->squareBlock(2), new_rows_a11, n1.ny);
          copy_row(this_vdof, m->squareBlock(1), new_rows_a12, n1.nx);
          copy_row(this_vdof, m->squareBlock(3), new_rows_a12, n1.ny);
          copy_row(this_vdof, m->rectBlock(0), new_rows_b1t, n1.nx);
          copy_row(this_vdof, m->rectBlock(1), new_rows_b1t, n1.ny);
          
          copy_row(this_vdof, m->squareBlock(0), new_rows_a21, n2.nx);
          copy_row(this_vdof, m->squareBlock(2), new_rows_a21, n2.ny);
          copy_row(this_vdof, m->squareBlock(1), new_rows_a22, n2.nx);
          copy_row(this_vdof, m->squareBlock(3), new_rows_a22, n2.ny);
          copy_row(this_vdof, m->rectBlock(0), new_rows_b2t, n2.nx);
          copy_row(this_vdof, m->rectBlock(1), new_rows_b2t, n2.ny);
        }
        else
        { // strong normal
          new_rows_a11[this_vdof][this_vdof] = n1.nx;
          new_rows_a12[this_vdof][this_vdof] = n1.ny;
          new_rows_b1t[this_vdof]; // set entire row to zero
          new_rows_a21[this_vdof][this_vdof] = n2.nx;
          new_rows_a22[this_vdof][this_vdof] = n2.ny;
          new_rows_b2t[this_vdof]; // set entire row to zero
        }
        
        double new_rhs1 = rhs[0]->at(this_vdof) * n1.nx
                          + rhs[0]->at(this_vdof + n_v_dof) * n1.ny;
        double new_rhs2 = rhs[0]->at(this_vdof) * n2.nx
                          + rhs[0]->at(this_vdof + n_v_dof) * n2.ny;
        rhs[0]->at(this_vdof) = new_rhs1;
        rhs[0]->at(this_vdof + n_v_dof) = new_rhs2;
      }*/
    }
  }
  m->get_A_block(0)->changeRows(new_rows_a11);
  m->get_A_block(1)->changeRows(new_rows_a12);
  m->get_A_block(2)->changeRows(new_rows_a21);
  m->get_A_block(3)->changeRows(new_rows_a22);
  m->get_BT_block(0)->changeRows(new_rows_b1t);
  m->get_BT_block(1)->changeRows(new_rows_b2t);
}

/** ************************************************************************ */
void StokesProblem::reorder_rows(std::shared_ptr<TMatrix> m, bool only_normal)
{
  const int ID = 
    this->get_velocity_space().GetCollection()->GetCell(0)->GetReference_ID();
  const int vActiveBound = this->get_velocity_space().GetActiveBound();
    
  std::map<int, std::map<int, double> > new_rows;
  
  int n_vdof = this->get_velocity_space().GetN_DegreesOfFreedom();
  if(this->get_rhs().length() - m->GetN_Rows() != 0)
    ErrThrow("wrong matrix. Currently only the Stokes coupling matrix can be ",
             "reordered ", m->GetN_Rows(), " ", this->get_rhs().length());
  
  const bool strong_normal = TDatabase::ParamDB->StoDa_weakGamma <= 0.0
                             && typeOf_bci == Dirichlet;
  only_normal = only_normal || strong_normal;
  
  if(strong_normal)
    ErrThrow("Dirichlet problem with strong normal in big system");
  
  for(unsigned int j = 0; j < interface.size(); j++)
  {
    TBaseCell *cell = interface[j]->GetNeighbour(0);
    if(cell->GetReference_ID() != ID)
      cell = interface[j]->GetNeighbour(1);
    // index of this edge in this cell (0,1, or 2 for triangles)
    
    const int cell_index = cell->GetCellIndex();
    int const * const vDOF = get_velocity_space().GetGlobalDOF(cell_index);
    FE2D vFEId = this->get_velocity_space().GetFE2D(0, cell);
    TBaseFunct2D *vBaseFunct = TFEDatabase2D::GetBaseFunct2DFromFE2D(vFEId);
    const int N_vBaseFunct = vBaseFunct->GetDimension();
    
    for(int lvdof = 0; lvdof < N_vBaseFunct; lvdof++) // local velocity dof
    {
      const int this_vdof = vDOF[lvdof];
      if(global_DOF_interface[this_vdof][0].second == -1)
        continue; // not an interface dof
      if(this_vdof >= vActiveBound)
        continue; // not an active dof (dirichlet boundary)
        
      if(new_rows.count(this_vdof) == 1
         || new_rows.count(this_vdof + n_vdof) == 1)
        continue; // this dof has been handled in neighboring interface edge
      
      const double nx = normal(this_vdof).nx;
      const double ny = normal(this_vdof).ny;
      const double tx = normal(this_vdof).tx;
      const double ty = normal(this_vdof).ty;
      const bool nt = nx*nx - ny*ny > -1e-12; // (nx * nx >= ny * ny);
      
      if(!only_normal || !nt)
      {
        // first component ( if nt==true, this is the normal component)
        copy_row(this_vdof,        m.get(), new_rows, nt ? nx : tx, this_vdof);
        copy_row(this_vdof+n_vdof, m.get(), new_rows, nt ? ny : ty, this_vdof);
      }
      else
      {
        new_rows[this_vdof]; // set entire row to zero
        //new_rows[this_vdof + n_vdof]; // set entire row to zero
      }
      if(!only_normal || nt)
      {
        // second component ( if nt==true, this is the tangential component)
        copy_row(this_vdof, m.get(), new_rows, nt ? tx : nx, 
                 this_vdof + n_vdof);
        copy_row(this_vdof + n_vdof, m.get(), new_rows, nt ? ty : ny,
                 this_vdof + n_vdof);
      }
      else
      {
        //new_rows[this_vdof]; // set entire row to zero
        new_rows[this_vdof + n_vdof]; // set entire row to zero
      }
      
      /*else
      { // kink
        // we have two normals here, one is (nx,ny) from the current (j-th) 
        // interface edge. The other is from the neighboring interface edge 
        
        // enforce two conditions, one for each of the two normal components
        
        // the two normals n1, n2
        n_t_vector n1 = normal(this_vdof, 0);
        n_t_vector n2 = normal(this_vdof, 1);
        copy_row(this_vdof,          m, new_rows, n1.nx, this_vdof);
        copy_row(this_vdof + n_vdof, m, new_rows, n1.ny, this_vdof);
        copy_row(this_vdof,          m, new_rows, n2.nx, this_vdof + n_vdof);
        copy_row(this_vdof + n_vdof, m, new_rows, n2.ny, this_vdof + n_vdof);
      }*/
    }
  }
  m->changeRows(new_rows);
}

/** ************************************************************************ */
void StokesProblem::prepare_map_solution_to_interface(
    const InterfaceFunction &eta)
{
  // the map which is used to map a computed solution to an interface function
  // will be affine linear (even linear for Neumann problems).
  
  /** create the structure of the underlying linear map (matrix[0]) */
  create_structure_of_map_solution_to_interface(eta.length());
  
  /** assemble the linear part of this map, and the affine part (vector) */
  // If this is a Neumann problem, Dirichlet values are returned. In this case
  // we have to assemble a mass matrix on the interface. If however this is a
  // Dirichlet problem, we return Neumann data, which means we have to assemble
  // integrals in the domain (not on the interface)
  // local matrices store all changes made to the linear part of this map
  // (matrix), both from integrals over the interface and over the domain
  local_matrices m;
  
  const int pt = TDatabase::ParamDB->StoDa_problemType;
  const int solution_strategy = TDatabase::ParamDB->StoDa_solutionStrategy;
  const int updating_strategy = TDatabase::ParamDB->StoDa_updatingStrategy;
  const bool C_RR = typeOf_bci == Robin && pt == 0 && updating_strategy == 3
                    && (solution_strategy == 1 || solution_strategy == -1);
  const bool D_RR = typeOf_bci == Robin && pt == 1 && updating_strategy == 4;
  
  if(typeOf_bci == Neumann || C_RR || D_RR)
  {
    // assemble a mass matrix on the interface (return Dirichlet dof)
    // if D-RR is true, we assemble (p-K grad(p).n, q), i.e. strong version of
    // Robin interface data
    assemble_Dirichlet_map_to_interface(m, D_RR);
  }
  else if(typeOf_bci == Dirichlet) // return Neumann data
  {
    map_sol2eta_rhs.copy_structure((const BlockVector&)eta);
    // assemble an integral in (parts of) the domain
    assemble_Neumann_map_to_interface(m);
    // assemble (u.n,v.n) on the interface
    assemble_Dirichlet_map_to_interface(m, D_RR); // D_RR is false here
  }
  else if(typeOf_bci == Robin || typeOf_bci == weakRobin)
  {
    map_sol2eta_rhs.copy_structure((const BlockVector&)eta);
    // D_RR should be false here
    // return -gamma_p * flux + pressure
    assemble_Dirichlet_map_to_interface(m, D_RR,
                                        TDatabase::ParamDB->StoDa_gamma_p); 
    assemble_Neumann_map_to_interface(m, -1.0);
  }
  else if(typeOf_bci == DirichletSTAB)
  {
    ErrThrow("DirichletSTAB not yet supported");
  }
  else
  {
    ErrThrow("unknown type of boundary condition on the interface");
  }
  map_sol2eta->add(m.m[0]);
}

/** ************************************************************************ */
void StokesProblem::find_global_DOF_interface(const InterfaceFunction &eta)
{
  if(global_DOF_interface.size() != 0)
    return; // no need to call this a second time
  
  // either do (1) or do not (-1) invert the direction
  // To each interface edge there is a direction, i.e. a vector pointing from
  // the first point to the second. This direction has to coincide with the 
  // direction of this edge in the neighboring Stokes cell. In there, the 
  // direction is counterclockwise. 
  // Each edge direction can only coincide with either the neighboring Stokes
  // or the neighboring Darcy cell.
  // if invertDirection is -1 then the reference transformations F_S for the 
  // Stokes dofs and F_I the interface dofs ((each) from [-1,1] to the edge 
  // [\vec a,\vec b]) have a different sign in front of their derivative: E.g.,
  // F_S(-1) = \vec a = F_I(1)  and  F_S(1) = \vec b = F_I(-1).
  // Therefore if invertDirection is -1 we simply use F_I(-z) instead of F_I(z)
  // to evaluate the interface function at the same quadrature point as in the 
  // Stokes part, F_S(z) = F_I(-z). This relies on the fact that the quadrature
  // points and weights are symmetric in [-1,1]
  int invertDirection = 1;
  
  global_DOF_interface.resize(
    this->get_velocity_space().GetN_DegreesOfFreedom());
  normals.resize(this->get_velocity_space().GetN_DegreesOfFreedom());
  // Id of the Stokes space
  const int S_id =
      this->get_velocity_space().GetCollection()->GetCell(0)->GetReference_ID();
  
  // everything needed for local assembling routines
  
  for(unsigned int iEdge = 0; iEdge < interface.size(); iEdge++)
  {
    const TInnerInterfaceJoint * thisEdge = interface[iEdge];
    TBaseCell *s_cell = thisEdge->GetNeighbour(0);
    if(s_cell->GetReference_ID() != S_id)
      s_cell = thisEdge->GetNeighbour(1);
    
    FE2D uFEId = this->get_velocity_space().GetFE2D(0, s_cell);
    TFE2D *uFE = TFEDatabase2D::GetFE2D(uFEId);
    const int N_uBf = uFE->GetN_DOF();
    int *uDOF = this->get_velocity_space().GetGlobalDOF(s_cell->GetCellIndex());
    // Basis functions for Stokes velocity
    TBaseFunct2D * uBf = uFE->GetBaseFunct2D();
    //RefTrans2D refTransS = uFE->GetRefTransID();
    
    //FE2D pFEId = pressure_space[0]->GetFE2D(0, s_cell);
    //TFE2D *pFE = TFEDatabase2D::GetFE2D(pFEId);
    //const int N_pBf = pFE->GetN_DOF();
    //int *pDOF = pressure_space[0]->GetGlobalDOF(s_cell->GetCellIndex());
    // Basis functions for Stokes velocity
    //TBaseFunct2D * pBf = pFE->GetBaseFunct2D();
    
    // number of local degrees of freedom on the interface 
    const int N_iBaseFunct = abs(eta.getSpaceType())+1;
    if(N_iBaseFunct != 3)
      ErrThrow("only quadratic interface functions supported so far");
    
    // index of this edge in the adjacent cell in Stokes subdomain
    const int eI = thisEdge->GetIndexInNeighbor(s_cell); // edge Index
    // normal (pointing out of the Stokes subdomain) and tangential 
    double nx, ny, tx, ty;
    getNormal(s_cell, thisEdge, nx, ny);
    n_t_vector this_normal(nx, ny);
    // check if edge is directed the same way as in the cell 's_cell'
    thisEdge->GetTangent(tx, ty);
    if(tx == -ny && ty == nx)
      invertDirection = 1;
    else
      //if(tx==ny && ty==-nx)
      invertDirection = -1;
    
    // find Stokes DOFs which correspond to the three interface DOFs on this 
    // edge, we do this only on the reference cell, which is enough since the 
    // transformation would not change the situation
    double points_on_interface[3] = {-1.0, 0.0, 1.0}; //points on reference edge
    if(invertDirection == -1)
    {
      points_on_interface[0] = 1.0;
      points_on_interface[2] = -1.0;
    }
    double *values[N_iBaseFunct];
    for(int i = 0; i < N_iBaseFunct; i++)
      values[i] = new double[N_uBf];
    uBf->GetValues(N_iBaseFunct, points_on_interface, eI, values);
    // now check which of the Stokes velocity basis functions are one at a point
    // on the interface and zero elsewhere
    for(int i = 0; i < N_iBaseFunct; i++)
    {
      const int iDOF = eta.getDOF(iEdge, i); // dof on interface
      for(int j = 0; j < N_uBf; j++)
      {
        if(values[i][j] == 1.0)
        {
          // the i-th interface dof corresponds to the j-th Stokes velocity dof
          global_DOF_interface[uDOF[j]].push_back(std::make_pair(s_cell, iDOF));
          if(normals[uDOF[j]].empty())
          {
            normals[uDOF[j]].push_back(this_normal);
          } 
          else if( normal(uDOF[j], 0) != this_normal )
             //&& (abs(normals[uDOF[j]].nx - nx) > 1e-4 
             //    || abs(normals[uDOF[j]].ny - ny) > 1e-4))
          { // we have a kink here
            //OutPut("old normal " << normal(uDOF[j], 0).nx << " " << 
            //       normal(uDOF[j], 0).ny << 
            //       " new " << this_normal.nx << " " << this_normal.ny << 
            //       " iDOF " << iDOF << endl);
            //const double old_nx = normals[uDOF[j]].nx;
            //const double old_ny = normals[uDOF[j]].ny;
            //const double one_over_det = 1 / (old_nx * ny - old_ny * nx);
            //normals[uDOF[j]].nx = (ny - old_ny) * one_over_det;
            //normals[uDOF[j]].ny = (old_nx - nx) * one_over_det;
            normals[uDOF[j]].push_back(this_normal);
          }
          // else // nothing more to do
          continue;
        }
      }
    }
    for(int i = 0; i < N_iBaseFunct; i++)
      delete[] values[i];
  }
  
  // find cells not sharing an edge with the interface but only a vertex and add
  // those to 'global_DOF_interface'
  TCollection * s_coll = this->get_velocity_space().GetCollection();
  const int n_cells = s_coll->GetN_Cells();
  for(int icell = 0; icell < n_cells; icell++)
  {
    TBaseCell *s_cell = s_coll->GetCell(icell);
    FE2D vFEId = this->get_velocity_space().GetFE2D(0, s_cell);
    TFE2D *vFE = TFEDatabase2D::GetFE2D(vFEId);
    const int N_vBf = vFE->GetN_DOF(); // number of Darcy pressure basis functions
    int *sDOF = this->get_velocity_space().GetGlobalDOF(s_cell->GetCellIndex());
    
    for(int j = 0; j < N_vBf; j++)
    {
      int s_dof = sDOF[j]; // this Stokes velocity dof
      // this vector of pairs
      std::vector<std::pair<TBaseCell*, int>>& tvop = 
        global_DOF_interface[s_dof];
      // check if this dof is on the interface and if so, enter this loop,
      if(tvop.empty())
        continue;
      // check if it is associated to this cell already
      bool already_included = false;
      for(unsigned int n = 0; n < tvop.size(); n++)
      {
        if(tvop[n].first == s_cell)
        {
          already_included = true;
          continue;
        }
      }
      if(!already_included)
      {
        int iDOF = tvop[0].second; // exists since tvop.size() is not zero
        // iDOF could be different in this cell for discontinuous interface
        // functions. That is only the case for D-RR, and then this is not used 
        // anyway. becuase integrals over the domain are not necessary for D-RR
        tvop.push_back(std::make_pair(s_cell, iDOF));
      }
    }
  }
  
  // fill entries in global_DOF_interface which are not yet filled with defaults
  for(int i = 0; i < this->get_velocity_space().GetN_DegreesOfFreedom(); i++)
  {
    if(global_DOF_interface[i].empty())
      global_DOF_interface[i].push_back(std::pair<TBaseCell*, int>(NULL, -1));
  }
}

/** ************************************************************************ */
void StokesProblem::create_structure_of_map_solution_to_interface(
    unsigned int l)
{
  // degrees of freedom for Stokes velocity and pressure 
  const int N_veloDOF = this->get_velocity_space().GetN_DegreesOfFreedom();
  
  // everything needed to call the constructor of TStructure
  const int n_rows = l;
  const int n_cols = this->get_rhs().length(); // = 2*N_veloDOF+N_presDOF;
  int N_entries = 0;
  int *cols, *rows = new int[n_rows + 1];
  
  // for each row (i.e, interface dof) we have one vector containing all 
  // (Stokes) dofs which couple with this dof. For all dofs which do not belong
  // to a cell adjacent to the interface, the corresponding columns will be of  
  // length 0
  std::vector<std::vector<int> > coupledDOF(n_rows);
  
  // array containing the number of basis functions for all finite elements
  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  
  // D-RR or regular Robin-Robin method 
  bool include_pressure = (typeOf_bci == Robin
                           && TDatabase::ParamDB->StoDa_problemType == 1)
                          || typeOf_bci == Dirichlet;
  
  // the collection of all Stokes cells
  const TCollection* s_coll = this->get_velocity_space().GetCollection();
  // number of Stokes cells
  const unsigned int n_cells = (unsigned int) s_coll->GetN_Cells();
  
  // loop over all cells
  for(unsigned int icell = 0; icell < n_cells; icell++)
  {
    TBaseCell *s_cell = s_coll->GetCell(icell);
    
    // velocity
    FE2D uFEId = this->get_velocity_space().GetFE2D(0, s_cell);
    const int N_uBaseFunct = N_BaseFunct[uFEId];
    int *uStokesDOF = get_velocity_space().GetGlobalDOF(s_cell->GetCellIndex());
    int *uDOF = this->get_velocity_space().GetGlobalDOF(s_cell->GetCellIndex());
    // pressure
    FE2D pFEId = this->get_pressure_space().GetFE2D(0, s_cell);
    const int N_pBaseFunct = N_BaseFunct[pFEId];
    int *pStokesDOF = get_pressure_space().GetGlobalDOF(s_cell->GetCellIndex());
    
    // loop over all velocity degrees of freedom of this Stokes cell
    for(int row = 0; row < N_uBaseFunct; row++)
    {
      int iDOF = global_DOF_interface[uDOF[row]][0].second;
      if(iDOF == -1)
        continue; // not an interface dof
      int l = 0;
      while(global_DOF_interface[uDOF[row]][l].first != s_cell)
        l++;
      // this dof belongs to more than one interface dof (this happens for
      // discontinuous interface functions), take the correct one.
      iDOF = global_DOF_interface[uDOF[row]][l].second;
      
      // loop over all velocity degrees of freedom of this Sokes cell
      for(int jS = 0; jS < N_uBaseFunct; jS++)
      {
        const int stokes_vDOF = uStokesDOF[jS];
        // first Stokes velocity component
        coupledDOF[iDOF].push_back(stokes_vDOF);
        // second Stokes velocity component
        coupledDOF[iDOF].push_back(stokes_vDOF + N_veloDOF);
      }
      if(include_pressure)
      for(int jS = 0; jS < N_pBaseFunct; jS++)
      {
        const int stokes_pDOF = pStokesDOF[jS];
        // Stokes pressure component
        coupledDOF[iDOF].push_back(stokes_pDOF + 2 * N_veloDOF);
      }
    }
  }
  
  // fill the array 'rows' and compute 'N_entries'
  rows[0] = 0;
  for(int i = 0; i < n_rows; i++)
  {
    // sort the interface dofs for the i-th interface dof
    std::sort(coupledDOF[i].begin(), coupledDOF[i].end());
    // remove all duplicates 
    std::vector<int>::iterator it = std::unique(coupledDOF[i].begin(),
                                                coupledDOF[i].end());
    // resize to the size without duplicates
    coupledDOF[i].resize(std::distance(coupledDOF[i].begin(), it));
    // number of Stokes dofs which couple with this interface dof  
    N_entries += coupledDOF[i].size();
    rows[i + 1] = N_entries;
  }
  
  cols = new int[N_entries];
  // fill the array 'cols'
  for(int i = 0; i < n_rows; i++)
  {
    for(int j = 0; j < rows[i + 1] - rows[i]; j++)
    {
      cols[rows[i] + j] = coupledDOF[i].at(j);
    }
  }
  // generate sparse matrix
  std::shared_ptr<TStructure> structure(new TStructure(n_rows, n_cols, 
                                                       N_entries, cols, rows));
  map_sol2eta.reset(new TMatrix(structure)); // empty matrix
}

/** ************************************************************************ */
void StokesProblem::assemble_Dirichlet_map_to_interface(local_matrices & m,
                                                        bool D_RR, double a)
{
  // number of velocity dofs
  const int N_SveloDOF = this->get_velocity_space().GetN_DegreesOfFreedom();
  
  if(D_RR)
  {
    // D-RR is somewhat difficult, so we separate it from the rest
    this->assemble_Dirichlet_map_to_interface_D_RR(m, a);
    return;
  }
  
  if(!assemble_on_return)
  {
    // instead of assembling terms on the interface, we only return the nodal
    // values and the assembling is done in the object which recieves the data
    // on the interface from this object
    for(int idof = 0; idof < N_SveloDOF; idof++)
    {
      int interface_dof = global_DOF_interface[idof][0].second;
      if(interface_dof != -1)
      {
        // the normal is not exactly a unit normal if this dof is at a kink
        m.m[0][interface_dof][idof             ] += normal(idof).nx;
        m.m[0][interface_dof][idof + N_SveloDOF] += normal(idof).ny;
      }
    }
    return;
  }
  // else: assemble terms 
  
  if(TDatabase::ParamDB->StoDa_interfaceType == 1 && typeOf_bci == Dirichlet)
    return; // no tangential mass term on interface needed
  
  // Id of the Stokes space
  const int S_id =
      this->get_velocity_space().GetCollection()->GetCell(0)->GetReference_ID();
  
  for(unsigned int iEdge = 0; iEdge < interface.size(); iEdge++)
  {
    const TInnerInterfaceJoint * thisEdge = interface[iEdge];
    TBaseCell *s_cell = thisEdge->GetNeighbour(0);
    if(s_cell->GetReference_ID() != S_id)
      s_cell = thisEdge->GetNeighbour(1);
    
    FE2D uFEId = this->get_velocity_space().GetFE2D(0, s_cell);
    TFE2D *uFE = TFEDatabase2D::GetFE2D(uFEId);
    const int N_uBf = uFE->GetN_DOF();
    int *uDOF = this->get_velocity_space().GetGlobalDOF(s_cell->GetCellIndex());
    // Basis functions for Stokes velocity
    TBaseFunct2D * uBf = uFE->GetBaseFunct2D();
    RefTrans2D refTransS = uFE->GetRefTransID();
    
    FE2D pFEId = this->get_pressure_space().GetFE2D(0, s_cell);
    TFE2D *pFE = TFEDatabase2D::GetFE2D(pFEId);
    const int N_pBf = pFE->GetN_DOF();
    int *pDOF = this->get_pressure_space().GetGlobalDOF(s_cell->GetCellIndex());
    // Basis functions for Stokes velocity
    TBaseFunct2D * pBf = pFE->GetBaseFunct2D();
    
    // set the coordinate transformation from the reference cell to this cell
    TFEDatabase2D::SetCellForRefTrans(s_cell, refTransS);
    
    // get the quadrature formula
    QuadFormula1D QFId; // quadrature formula id
    int N_LinePoints; // number of qudrature points
    double *LineWeights, *zeta; // quadrature weights and points (on [-1,1])
    {
      // polynomial degree of finite element, needed for choosing an 
      // appropriate quadrature formula, 2 is the polynomial degree of functions 
      // on the interface
      const int fe_degree = uBf->GetPolynomialDegree() * 2;
      QFId = TFEDatabase2D::GetQFLineFromDegree(fe_degree);
      TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(QFId);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
      // qf1 no longer needed, only local scope
    }
    
    // make sure all functions & derivatives are available for this quadrature
    uBf->MakeRefElementData(QFId);
    pBf->MakeRefElementData(QFId);
    
    // compute length of the edge
    const double hE = thisEdge->GetLength();
    // index of this edge in the adjacent cell in Stokes subdomain
    const int eI = thisEdge->GetIndexInNeighbor(s_cell); // edge Index
    // normal (pointing out of the Stokes subdomain) and tangential 
    double nx, ny, tx, ty;
    getNormal(s_cell, thisEdge, nx, ny);
    thisEdge->GetTangent(tx, ty);
    
    // values of all functions and derivatives on reference element 
    double **uref = TFEDatabase2D::GetJointValues2D(uBf->GetID(), QFId, eI);
    double **pref = TFEDatabase2D::GetJointValues2D(pBf->GetID(), QFId, eI);
    
    // values and derivatives of all functions at one quadrature point
    double *uorig = new double[N_uBf];
    double *porig = new double[N_pBf];
    
    for(int k = 0; k < N_LinePoints; k++)
    {
      // quadrature weight and determinant of tranformation
      const double qw = LineWeights[k] * (hE / 2);
      
      TFEDatabase2D::GetOrigValues(refTransS, zeta[k], uBf, eI, uref[k], NULL,
                                   NULL, uorig, NULL, NULL);
      TFEDatabase2D::GetOrigValues(refTransS, zeta[k], pBf, eI, pref[k], NULL,
                                   NULL, porig, NULL, NULL);
      
      for(int row = 0; row < N_uBf; row++)
      {
        int testDOF = global_DOF_interface[uDOF[row]][0].second;
        if(testDOF == -1)
          continue; // not an interface dof
        if(global_DOF_interface[uDOF[row]][0].first != s_cell)
          // this dof belongs to two interface dofs (this happens for
          // discontinuous interface functions), take the correct one.
          testDOF = global_DOF_interface[uDOF[row]][1].second;
        
        const double v = uorig[row];
        //const double vx = uxorig[row];
        //const double vy = uyorig[row];
        
        //nx = normal(uDOF[row]).nx;
        //ny = normal(uDOF[row]).ny;
        
        // velocity
        for(int col = 0; col < N_uBf; col++)
        {
          const int ansatzDOF = uDOF[col];
          const double u = uorig[col];
          // (u.t,v.t) 
          if(TDatabase::ParamDB->StoDa_interfaceType == 0 
             && (typeOf_bci == Dirichlet || typeOf_bci == Robin))
          {
            if(kink(uDOF[row]))
            {
              // kink, we set n to neither one of the two possible normals, 
              // n1, n2, but to a (uniqe) linear combination such that n.n1 = 1 
              // and n.n2 = 1.
              
              // this 'n' is not unit normal vector. The following is an 
              // appropriate test function vector: v *(nx, ny)
              nx = normal(uDOF[row]).nx;
              ny = normal(uDOF[row]).ny;
              // ignore 'a' here because it is gamma_p for Robin problems
              double test = v * (nx * tx + ny * ty);
              if(typeOf_bci == Robin)
                test *= -1.0;
              m.m[0][testDOF][ansatzDOF] += test * u * tx * qw;
              m.m[0][testDOF][ansatzDOF + N_SveloDOF] += test * u * ty * qw;
            }
            if(typeOf_bci == Dirichlet)
              continue; // further assembling only for Robin problems
          }
          m.m[0][testDOF][ansatzDOF             ] += a * v * u * qw * nx;
          m.m[0][testDOF][ansatzDOF + N_SveloDOF] += a * v * u * qw * ny;
        }
        
        // pressure
        if(D_RR)
        {
          for(int col = 0; col < N_pBf; col++)
          {
            const int ansatzDOF = pDOF[col];
            m.m[0][testDOF][ansatzDOF + 2 * N_SveloDOF] += a * porig[col]
                                                           * uorig[row] * qw;
          }
        }
      }
    }
    delete[] uorig;
    delete[] porig;
  }
}

/** ************************************************************************ */
void StokesProblem::assemble_Dirichlet_map_to_interface_D_RR(local_matrices & m,
                                                             double a)
{
  // in case of assemble_on_return == true, we have to use a quadrature rule
  // and evaluate the basis functions at the quadrature points
  // if assemble_on_return == false we have to evaluate the basis functions
  // at the three points on each edge.
  // so either way we have to evaluate all basis functions at points on the 
  // interface.
  // Note that if assemble_on_return == false the interface space should be
  // discontinuous!
  
  // number of velocity dofs
  const int N_SveloDOF = this->get_velocity_space().GetN_DegreesOfFreedom();
  const double nu = 1.0 / TDatabase::ParamDB->RE_NR; // viscosity
                    
  // Id of the Stokes space
  const int S_id =
      this->get_velocity_space().GetCollection()->GetCell(0)->GetReference_ID();
  
  for(unsigned int iEdge = 0; iEdge < interface.size(); iEdge++)
  {
    const TInnerInterfaceJoint * thisEdge = interface[iEdge];
    TBaseCell *s_cell = thisEdge->GetNeighbour(0);
    if(s_cell->GetReference_ID() != S_id)
      s_cell = thisEdge->GetNeighbour(1);
    
    FE2D uFEId = this->get_velocity_space().GetFE2D(0, s_cell);
    TFE2D *uFE = TFEDatabase2D::GetFE2D(uFEId);
    const int N_uBf = uFE->GetN_DOF();
    int *uDOF = this->get_velocity_space().GetGlobalDOF(s_cell->GetCellIndex());
    // Basis functions for Stokes velocity
    TBaseFunct2D * uBf = uFE->GetBaseFunct2D();
    RefTrans2D refTransS = uFE->GetRefTransID();
    
    FE2D pFEId = this->get_pressure_space().GetFE2D(0, s_cell);
    TFE2D *pFE = TFEDatabase2D::GetFE2D(pFEId);
    const int N_pBf = pFE->GetN_DOF();
    int *pDOF = this->get_pressure_space().GetGlobalDOF(s_cell->GetCellIndex());
    // Basis functions for Stokes velocity
    TBaseFunct2D * pBf = pFE->GetBaseFunct2D();
    
    // set the coordinate transformation from the reference cell to this cell
    TFEDatabase2D::SetCellForRefTrans(s_cell, refTransS);
    
    // get the quadrature formula
    QuadFormula1D QFId; // quadrature formula id
    int N_LinePoints; // number of qudrature points
    double *LineWeights, *zeta; // quadrature weights and points (on [-1,1])
    if(assemble_on_return)
    {
      // do a quadrature
      
      // polynomial degree of finite element, needed for choosing an 
      // appropriate quadrature formula, 2 is the polynomial degree of functions 
      // on the interface
      const int fe_degree = uBf->GetPolynomialDegree() * 2;
      QFId = TFEDatabase2D::GetQFLineFromDegree(fe_degree);
      TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(QFId);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
      // qf1 no longer needed, only local scope
    }
    else
    {
      // point evaluations at three nodes on this edge
      N_LinePoints = 3;
      LineWeights = new double[3];
      LineWeights[0] = LineWeights[1] = LineWeights[2] = 1.;
      zeta = new double[3];
      zeta[0] = -1.;
      zeta[1] = 0.;
      zeta[2] = 1.;
    }
    
    // make sure all functions & derivatives are available for this quadrature
    uBf->MakeRefElementData(QFId);
    pBf->MakeRefElementData(QFId);
    
    // compute length of the edge
    const double hE = thisEdge->GetLength();
    // index of this edge in the adjacent cell in Stokes subdomain
    const int eI = thisEdge->GetIndexInNeighbor(s_cell); // edge Index
    // normal (pointing out of the Stokes subdomain) and tangential 
    double nx, ny, tx, ty;
    getNormal(s_cell, thisEdge, nx, ny);
    thisEdge->GetTangent(tx, ty);
    
    // values of all functions and derivatives on reference element
    double **uref, **uxiref, **uetaref, **pref;
    if(assemble_on_return)
    {
      uref = TFEDatabase2D::GetJointValues2D(uBf->GetID(), QFId, eI);
      uxiref = TFEDatabase2D::GetJointDerivatives2D(uBf->GetID(), QFId, eI,
                                                    D10);
      uetaref = TFEDatabase2D::GetJointDerivatives2D(uBf->GetID(), QFId, eI,
                                                     D01);
      pref = TFEDatabase2D::GetJointValues2D(pBf->GetID(), QFId, eI);
    }
    else
    {
      // we could define a quadrature rule with only the three points, but
      // the functions 'TFEDatabase2D::GetJointDerivatives2D' only know 
      // predefined quadrature formulas, we have to do this on our own here
      // disadvantage: these results are not stored in the database and have to
      // be recomputed on every interface edge
      uref = new double*[N_LinePoints]; // == 3
      uxiref = new double*[N_LinePoints]; // == 3
      uetaref = new double*[N_LinePoints]; // == 3
      pref = new double*[N_LinePoints]; // == 3
      for(int i = 0; i < N_LinePoints; i++)
      {
        uref[i] = new double[N_uBf];
        uxiref[i] = new double[N_uBf];
        uetaref[i] = new double[N_uBf];
        pref[i] = new double[N_pBf];
      }
      uBf->GetDerivatives(D00, N_LinePoints, zeta, eI, uref);
      uBf->GetDerivatives(D10, N_LinePoints, zeta, eI, uxiref);
      uBf->GetDerivatives(D01, N_LinePoints, zeta, eI, uetaref);
      pBf->GetDerivatives(D00, N_LinePoints, zeta, eI, pref);
    }
    // values and derivatives of all functions at one quadrature point
    double *uorig = new double[N_uBf];
    double *uxorig = new double[N_uBf];
    double *uyorig = new double[N_uBf];
    double *porig = new double[N_pBf];
    
    for(int k = 0; k < N_LinePoints; k++)
    {
      // quadrature weight and determinant of tranformation
      const double qw = LineWeights[k] * (assemble_on_return ? (hE / 2) : 1.0);
      
      TFEDatabase2D::GetOrigValues(refTransS, zeta[k], uBf, eI, uref[k],
                                   uxiref[k], uetaref[k], uorig, uxorig,
                                   uyorig);
      TFEDatabase2D::GetOrigValues(refTransS, zeta[k], pBf, eI, pref[k], NULL,
                                   NULL, porig, NULL, NULL);
      
      for(int row = 0; row < N_uBf; row++)
      {
        int testDOF = global_DOF_interface[uDOF[row]][0].second;
        if(testDOF == -1)
          continue; // not an interface dof
        if(global_DOF_interface[uDOF[row]][0].first != s_cell)
          // this dof belongs to two interface dofs (this happens for
          // discontinuous interface functions), take the correct one.
          testDOF = global_DOF_interface[uDOF[row]][1].second;
        
        const double v = uorig[row];
        //const double vx = uxorig[row];
        //const double vy = uyorig[row];
        
        //nx = normal(uDOF[row]).nx;
        //ny = normal(uDOF[row]).ny;
        
        // velocity
        for(int col = 0; col < N_uBf; col++)
        {
          const int ansatzDOF = uDOF[col];
          const double u = uorig[col];
          const double ux = uxorig[col];
          const double uy = uyorig[col];
          
          double b = u * TDatabase::ParamDB->StoDa_gamma_p
                     - 2 * nu * (ux * nx + uy * ny);
          m.m[0][testDOF][ansatzDOF] += a * b * v * qw * nx;
          m.m[0][testDOF][ansatzDOF + N_SveloDOF] += a * b * v * qw * ny;
        }// end for loop over ansatz velocity basis functions
        
        // pressure
        for(int col = 0; col < N_pBf; col++)
        {
          const int ansatzDOF = pDOF[col];
          m.m[0][testDOF][ansatzDOF + 2 * N_SveloDOF] += a * porig[col]
                                                         * uorig[row] * qw;
        } // end for loop over ansatz pressure basis functions
      } // end for loop over test velocity basis functions
    } // end for loop over all quadrature points
    delete[] uorig;
    delete[] uxorig;
    delete[] uyorig;
    delete[] porig;
    if(!assemble_on_return)
    {
      // delete arrays which have been created on the heap using new ...[];
      delete[] LineWeights;
      delete[] zeta;
      for(int i = 0; i < N_LinePoints; i++) // i < 3
      {
        delete[] uref[i];
        delete[] uxiref[i];
        delete[] uetaref[i];
        delete[] pref[i];
      }
      delete[] uref;
      delete[] uxiref;
      delete[] uetaref;
      delete[] pref;
    } // end if !assemble_on_return
  } // end for loop over all interface edges
}

/** ************************************************************************ */
void StokesProblem::assemble_Neumann_map_to_interface(local_matrices &m,
                                                      double a)
{
  if(assemble_on_return == false)
    // when returning Neumann data we have to assemble integrals in the domain.
    ErrThrow("Trying to return Neumann data without assembling integrals");
  
  const double nu = TDatabase::ParamDB->RE_NR;
  
  // the collection of all Darcy cells
  TCollection* s_coll = this->get_velocity_space().GetCollection();
  // number of Darcy cells
  const unsigned int n_cells = (unsigned int) s_coll->GetN_Cells();
  const int n_active = this->get_velocity_space().GetActiveBound();
  const int n_v_dof = this->get_velocity_space().GetN_DegreesOfFreedom();
  
  // loop over all cells, skip cells which are not at the interface, i.e., cells
  // which have no dof on the interface. That means cells which only share a 
  // vertex with the interface are considered, even if they have no common edge
  // with the interface. Then assemble the necessary terms
  for(unsigned int icell = 0; icell < n_cells; icell++)
  {
    TBaseCell *s_cell = s_coll->GetCell(icell);
    
    FE2D uFEId = this->get_velocity_space().GetFE2D(0, s_cell);
    TFE2D *uFE = TFEDatabase2D::GetFE2D(uFEId);
    const int N_uBf = uFE->GetN_DOF();
    int *uDOF = this->get_velocity_space().GetGlobalDOF(s_cell->GetCellIndex());
    // Basis functions for Darcy velocity
    TBaseFunct2D * uBf = uFE->GetBaseFunct2D();
    
    // check if any dofs are on the interface
    bool on_interface = false;
    for(int row = 0; row < N_uBf; row++)
    {
      if(global_DOF_interface[uDOF[row]][0].second != -1)
      {
        on_interface = true;
        break;
      }
    }
    if(!on_interface)
      continue;
    
    FE2D pFEId = this->get_pressure_space().GetFE2D(0, s_cell);
    TFE2D *pFE = TFEDatabase2D::GetFE2D(pFEId);
    const int N_pBf = pFE->GetN_DOF(); // # Stokes pressure basis functions
    int *pDOF = this->get_pressure_space().GetGlobalDOF(s_cell->GetCellIndex());
    // Basis functions for Stokes pressure
    TBaseFunct2D * pBf = pFE->GetBaseFunct2D();
    
    RefTrans2D refTrans = pFE->GetRefTransID();
    TFEDatabase2D::SetCellForRefTrans(s_cell, refTrans);
    
    const int fe_degree = 2 * uBf->GetPolynomialDegree(); // 2 -> bilinear form
    QuadFormula2D QF_2D_Id = TFEDatabase2D::GetQFFromDegree(
        fe_degree, TFEDatabase2D::GetRefElementFromFE2D(uFEId));
    int N_qPoints; // number of qudrature points
    double *qWs, *x, *y; // quadrature weights and points on reference cell
    {
      TQuadFormula2D *qf2 = TFEDatabase2D::GetQuadFormula2D(QF_2D_Id);
      qf2->GetFormulaData(N_qPoints, qWs, x, y);
      // qf2 no longer needed, only local scope
    }
    
    // make sure all functions & derivatives are available for this quadrature
    uBf->MakeRefElementData(QF_2D_Id);
    pBf->MakeRefElementData(QF_2D_Id);
    
    // get values and derivatives on reference cell
    double **uref = TFEDatabase2D::GetRefElementValues(uBf->GetID(), QF_2D_Id,
                                                       D00);
    double **uxiref = TFEDatabase2D::GetRefElementValues(uBf->GetID(), QF_2D_Id,
                                                         D10);
    double **uetaref = TFEDatabase2D::GetRefElementValues(uBf->GetID(),
                                                          QF_2D_Id, D01);
    double **pref = TFEDatabase2D::GetRefElementValues(pBf->GetID(), QF_2D_Id,
                                                       D00);
    
    double *uorig = new double[N_uBf];
    double *uxorig = new double[N_uBf];
    double *uyorig = new double[N_uBf];
    double *porig = new double[N_pBf];
    
    double *X, *Y, *det; // Points in this cell, determinant of the transformation
    X = new double[N_qPoints];
    Y = new double[N_qPoints];
    det = new double[N_qPoints];
    TFEDatabase2D::GetOrigFromRef(refTrans, N_qPoints, x, y, X, Y, det);
    
    // assemble
    for(int k = 0; k < N_qPoints; k++)// loop over all quadrature points
    {
      // transformation to this grid cell
      // we use NULL in the 'TGridCell'-argument because it is not needed. 
      TFEDatabase2D::GetOrigValues(refTrans, x[k], y[k], uBf, s_coll, NULL,
                                   uref[k], uxiref[k], uetaref[k], uorig,
                                   uxorig, uyorig);
      TFEDatabase2D::GetOrigValues(refTrans, x[k], y[k], pBf, s_coll, NULL,
                                   pref[k], NULL, NULL, porig, NULL, NULL);
      
      double qw = qWs[k] * det[k];
      
      for(int row = 0; row < N_uBf; row++)
      {
        // local interface dof which is exactly this Stokes velocity dof (on the 
        // interface)
        int testDOF = global_DOF_interface[uDOF[row]][0].second;
        if(testDOF == -1)
          continue; // not an interface dof
        if(global_DOF_interface[uDOF[row]][0].first != s_cell)
          // this dof belongs to two interface dofs (this happens for
          // discontinuous interface functions), take the correct one.
          testDOF = global_DOF_interface[uDOF[row]][1].second;
        
        // normal for this dof
        const double nx = normal(uDOF[row]).nx; 
        const double ny = normal(uDOF[row]).ny;
        // the test function is uorig[row]*(nx, ny)
        //const double v  = uorig[row];
        const double vx = uxorig[row];
        const double vy = uyorig[row];
        
        for(int col = 0; col < N_uBf; col++)
        {
          const int ansatzDOF = uDOF[col];
          //const double u = uorig[col];
          const double ux = uxorig[col];
          const double uy = uyorig[col];
          // a_f(u,v), where the test function v = uorig[row]*(nx, ny)
          m.m[0][testDOF][ansatzDOF] += a * nu * (2 * ux * vx + uy * vy) * nx
                                        * qw;
          m.m[0][testDOF][ansatzDOF] += a * nu * uy * vx * ny * qw;
          m.m[0][testDOF][ansatzDOF + n_v_dof] += a * nu * ux * vy * nx * qw;
          m.m[0][testDOF][ansatzDOF + n_v_dof] += a * nu
                                                  * (2 * uy * vy + ux * vx) * ny
                                                  * qw;
        }
        for(int col = 0; col < N_pBf; col++)
        {
          const int ansatzDOF = pDOF[col] + 2 * n_v_dof;
          // b_f(v,p)
          m.m[0][testDOF][ansatzDOF] -= a * porig[col] * vx * nx * qw;
          m.m[0][testDOF][ansatzDOF] -= a * porig[col] * vy * ny * qw;
        }
        
        // take care of 'map_sol2eta_rhs'
        if(uDOF[row] > n_active) // Dirichlet (nonactive) dof
        {
          // we need to assemble (f,v) here, which is not present in the 
          // right hand side (it is only present for active dofs)
          
          // the call to the linear Coefficient function should be done before
          // the loop over the quadrature points. However we only want to call
          // that function when it is necessary. So here we do the entire
          // quadrature (ugly)
          if(k == 0)
          {
            // note that, if you take the solution of this Dirichlet Stokes 
            // problem and map it to eta, using 'map_solution_to_interface', 
            // then take that eta as Neumann boundary data of a Stokes problem,
            // the non-active dofs are negelected anyway. That means you can 
            // do anything in this entry in eta. This quadrature is therefore
            // not necessary in that case. However it might be important, if 
            // eta is passed to a Darcy problem as interface data.
            double **coeffs = new double*[N_qPoints];
            for(int i = 0; i < N_qPoints; i++)
              coeffs[i] = new double[4];
            
            CoeffFct2D* linCoeffs = get_example().get_coeffs();
            linCoeffs(N_qPoints, X, Y, NULL, coeffs);
            double fvx = 0.0, fvy = 0.0; // assemble (f,v) into this double
            for(int i = 0; i < N_qPoints; i++)
            {
              // note that pref and porig are the same, otherwise we could not
              // access porig at other quadrature points than the 'k'-th
              fvx += (coeffs[i][1] * uref[i][row]) * nx * qWs[i] * det[i];
              fvy += (coeffs[i][2] * uref[i][row]) * ny * qWs[i] * det[i];
            }
            map_sol2eta_rhs(testDOF) += a *fvx;
            map_sol2eta_rhs(testDOF) += a *fvy;
            
            // delete created arrays
            for(int i = 0; i < N_qPoints; i++)
              delete[] coeffs[i];
            delete[] coeffs;
          }
        }
        else
        {
          if(k == 0)
          {
            // only do this once (not for every quadrature point)
            const bool nt = nx*nx - ny*ny > -1e-12; // nx*nx >= ny*ny;
            //map_sol2eta_rhs(testDOF) = a * this->get_rhs()[uDOF[row]] * nx;
            //map_sol2eta_rhs(testDOF) += 
            //                         a*this->get_rhs()[uDOF[row]+n_v_dof]*ny;
            double val = a * this->get_rhs()[uDOF[row] + (nt ? 0 : n_v_dof)];
            map_sol2eta_rhs(testDOF) = val;
          }
        }
      }
    }
    delete [] det;
    delete [] X;
    delete [] Y;
    delete [] uorig;
    delete [] uxorig;
    delete [] uyorig;
    delete [] porig;
  }
  if (TDatabase::ParamDB->StoDa_weakGamma <= 0.0)
  {
    // remove entries in rhs, so that only eta will be there
    for(int i = 0; i < n_active; i++)
    {
      int testDOF = global_DOF_interface[i][0].second;
      if(testDOF == -1)
        continue; // not an interface dof
      const bool nt = POW(normal(i).nx, 2) - POW(normal(i).tx, 2) > -1e-12;
      //const bool nt = (POW(normal(i).nx, 2) >= POW(normal(i).tx, 2));
      
      if(typeOf_bci == Dirichlet)
      {
        // reset value in normal component (eta will be added here later)
        this->get_rhs()[i + (nt ? 0 : n_v_dof)] = 0;
        DrhsNSE[i + (nt ? 0 : n_v_dof)] = 0;
      }
      if(getIcond() == 1)
      {
        // reset value in tangential component (u.t = 0)
        this->get_rhs()[i + (nt ? n_v_dof : 0)] = 0;
        DrhsNSE[i + (nt ? n_v_dof : 0)] = 0;
      }
      /*
       { // kink
       this->get_rhs()[i] = 0;
       DrhsNSE[i] = 0;
       this->get_rhs()[i + n_v_dof] = 0;
       [i + n_v_dof] = 0;
       }*/
    }
  }
}

/** ************************************************************************ */
void StokesProblem::map_solution_to_interface(InterfaceFunction &eta, double a,
                                              bool old)
{
  if(map_sol2eta == NULL)
  {
    ErrMsg("map_sol2eta has not been assembled yet in a " << typeOf_bci << 
           " (Navier-) Stokes problem. It is done now.");
    prepare_map_solution_to_interface(eta);
  }
  
  const double* const solution = old ? getSolOld() 
                                     : get_solution().get_entries();
  //const double* const solution = solution_direct->get_entries();
  
  // multiply solution with linear part of the map
  map_sol2eta->multiply(solution, eta.get_entries(), a);
  
  const int pt = TDatabase::ParamDB->StoDa_problemType;
  const int solution_strategy = TDatabase::ParamDB->StoDa_solutionStrategy;
  const int updating_strategy = TDatabase::ParamDB->StoDa_updatingStrategy;
  const bool C_RR = typeOf_bci == Robin && pt == 0 && updating_strategy == 3
                    && (solution_strategy == 1 || solution_strategy == -1);
  const bool D_RR = typeOf_bci == Robin && pt == 1 && updating_strategy == 4;
  if(typeOf_bci == Neumann || C_RR || D_RR)
  {} //there is nothing more to do
  else if(typeOf_bci == Dirichlet)
  {
    // substract right hand side, -1 means full vector
    eta.BlockVector::add(map_sol2eta_rhs.get_entries(), -1, -a);
  }
  else if(typeOf_bci == Robin || typeOf_bci == weakRobin)
  {
    // substract right hand side, -1 means full vector
    eta.BlockVector::add(map_sol2eta_rhs.get_entries(), -1, -a);
  }
  else if(typeOf_bci == DirichletSTAB)
  {
    ErrThrow("DirichletSTAB not yet supported");
  }
  else
  {
    ErrThrow("unknown type of boundary condition on the interface");
  }
}

/** ************************************************************************ */
void StokesProblem::update(InterfaceFunction& eta_p, InterfaceFunction* eta_f)
{
  // update eta_p
  const int solution_strategy = TDatabase::ParamDB->StoDa_solutionStrategy;
  
  if(solution_strategy == 1 || solution_strategy == -1)
  {
    const int updatingProcedure = TDatabase::ParamDB->StoDa_updatingStrategy;
    const double theta = TDatabase::ParamDB->StoDa_theta_p; // damping
    // regular Neumann-Neumann, Robin-Robin scheme, or C-RR or D-RR
    switch(updatingProcedure)
    {
      case 1:
      {
        // damping with old interface function
        eta_p.scale(1 - theta);
        this->map_solution_to_interface(eta_p, theta);
        break;
      }
      case 2: // damping with old Stokes solution
      {
        eta_p.reset();
        this->map_solution_to_interface(eta_p, theta);
        this->map_solution_to_interface(eta_p, 1.0 - theta, true);
        break;
      }
      case 3: // C-RR
      {
        if(eta_f == NULL)
          ErrThrow("For the C-RR method, you have to specify eta_f");
        
        double gamma_sum = TDatabase::ParamDB->StoDa_gamma_f
                           + TDatabase::ParamDB->StoDa_gamma_p;
        eta_p.scale(1 - theta);
        this->map_solution_to_interface(eta_p, theta * gamma_sum);
        // substract eta_f
        eta_p.add(eta_f, -theta);
        break;
      }
      case 4: // D-RR
        // get Dirichlet data and normal component of normal stress
        eta_p.scale(1 - theta);
        this->map_solution_to_interface(eta_p, theta);
        break;
      default:
        ErrThrow("unknown updating strategy");
        break;
    }
  }
  else if(solution_strategy == 2 || solution_strategy == -2)
  {
    // solving fixed point equation on interface 
    // turn off damping here, it is done in 'solve_fixed_point(eta)'
    eta_p.reset();
    this->map_solution_to_interface(eta_p, 1.0);
  }
  else if(solution_strategy == 3 || solution_strategy == -3)
  {
    // solving a Stecklov-Poincare equation
    // this is done without this update function
    ErrThrow("Stecklov-Poincare");
    this->map_solution_to_interface(eta_p, 1.0);
  }
}


/** ************************************************************************ */
void StokesProblem::create_etaToBd(const InterfaceFunction &eta)
{
  // degrees of freedom for Stokes velocity and pressure 
  const int N_veloDOF = this->get_velocity_space().GetN_DegreesOfFreedom();
  const int N_iDOF = eta.length(); // number of degrees of freedom on interface
      
  // everything needed to call the constructor of TStructure
  const int n_rows = this->get_rhs().length(); // = 2*N_veloDOF+N_presDOF;
  const int n_cols = N_iDOF;
  int N_entries = 0;
  int *cols, *rows = new int[n_rows + 1];
  
  // for each row (i.e, Stokes dof) we have one vector containing all 
  // (interface) dofs which couple with this dof. For all dofs which do not 
  // belong to a cell adjacent to the interface, the corresponding vector will 
  // be of length 0
  std::vector<std::vector<int> > coupledDOF(n_rows);
  
  // array containing the number of basis functions for all finite elements
  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  
  // Id of the Stokes space
  const int S_id =
     this->get_velocity_space().GetCollection()->GetCell(0)->GetReference_ID();
  
  // loop over the interface edges
  for(unsigned int iEdge = 0; iEdge < interface.size(); iEdge++)
  {
    const TInnerInterfaceJoint * thisEdge = interface[iEdge];
    TBaseCell *s_cell = thisEdge->GetNeighbour(0);
    if(s_cell->GetReference_ID() != S_id)
      s_cell = thisEdge->GetNeighbour(1);
    
    // velocity
    FE2D uFEId = this->get_velocity_space().GetFE2D(0, s_cell);
    const int N_uBaseFunct = N_BaseFunct[uFEId];
    int *uStokesDOF = 
       this->get_velocity_space().GetGlobalDOF(s_cell->GetCellIndex());
    // pressure
    FE2D pFEId = this->get_pressure_space().GetFE2D(0, s_cell);
    const int N_pBaseFunct = N_BaseFunct[pFEId];
    int *pStokesDOF = 
      this->get_pressure_space().GetGlobalDOF(s_cell->GetCellIndex());
    // number of basis functions on each interface edge
    const int N_iBaseFunct = 3; //abs(eta.spaceType)+1;
    
    // loop over all degrees of freedom of this interface function on this edge
    for(int j = 0; j < N_iBaseFunct; j++)
    {
      const int iDOF = eta.getDOF(iEdge, j);
      // loop over all velocity degrees of freedom of this Sokes cell
      for(int jS = 0; jS < N_uBaseFunct; jS++)
      {
        const int ansatz_vDOF = uStokesDOF[jS];
        // first Stokes velocity component
        coupledDOF[ansatz_vDOF].push_back(iDOF);
        // second Stokes velocity component
        coupledDOF[ansatz_vDOF + N_veloDOF].push_back(iDOF);
      }
      for(int jS = 0; jS < N_pBaseFunct; jS++)
      {
        const int ansatz_pDOF = pStokesDOF[jS];
        // Stokes pressure component
        coupledDOF[ansatz_pDOF + 2 * N_veloDOF].push_back(iDOF);
      }
    }
  }
  
  // fill the array 'rows' and compute 'N_entries'
  rows[0] = 0;
  for(int i = 0; i < n_rows; i++)
  {
    // sort the interface dofs for the i-th Stokes dof
    std::sort(coupledDOF[i].begin(), coupledDOF[i].end());
    // remove all duplicates 
    std::vector<int>::iterator it = std::unique(coupledDOF[i].begin(),
                                                coupledDOF[i].end());
    // resize to the size without duplicates
    coupledDOF[i].resize(std::distance(coupledDOF[i].begin(), it));
    // number of interface dofs which couple with this Stokes dof  
    N_entries += coupledDOF[i].size();
    rows[i + 1] = N_entries;
  }
  
  cols = new int[N_entries];
  // fill the array 'cols'
  for(int i = 0; i < n_rows; i++)
  {
    for(int j = 0; j < rows[i + 1] - rows[i]; j++)
    {
      cols[rows[i] + j] = coupledDOF[i].at(j);
    }
  }
  // generate sparse matrix
  std::shared_ptr<TStructure> structure(new TStructure(n_rows, n_cols, 
                                                       N_entries, cols, rows));
  etaToBd.reset(new TMatrix(structure)); // empty matrix
      
  Assemble_etaToBd(eta);
  if(TDatabase::ParamDB->StoDa_periodicBoundary)
    makePeriodicBoundary(etaToBd);
}

/** ************************************************************************ */
void StokesProblem::Assemble_etaToBd(const InterfaceFunction &eta)
{
  // number of Stokes velocity degrees of freedom for each component
  const int N_Sactive = this->get_velocity_space().GetActiveBound();
  const int N_SveloDOF = this->get_velocity_space().GetN_DegreesOfFreedom();
  
  if(!assemble_on_input)
  {
    for(int row = 0; row < N_SveloDOF; row++)
    {
      int i_dof = global_DOF_interface[row][0].second;
      if(i_dof != -1 && row < N_Sactive)
      {
        const double nx = normal(row).nx;
        const double tx = normal(row).tx;
        // set entries such that larger numbers are on the diagonal
        const bool nt = nx*nx - tx*tx > -1e-12; // (nx * nx >= tx * tx);
        for(unsigned int i = 0; i < global_DOF_interface[row].size(); i++)
        {
          i_dof = global_DOF_interface[row][i].second;
          (*etaToBd)(row + (nt ? 0 : N_SveloDOF), i_dof) = normal(row).norm();
          /*else // kink
          {
            (*etaToBd)(row,              i_dof) = 1.0;
            (*etaToBd)(row + N_SveloDOF, i_dof) = 1.0;
          }*/
          //SIGN((nt ? nx : -tx));
          //normals[row].nx + normals[row].ny;
        }
      }
    }
  }
  else
  {
    if(typeOf_bci == Dirichlet && TDatabase::ParamDB->StoDa_weakGamma <= 0.0)
    {
      // do essentially the same as in the 'assemble_on_input==false' case
      for(int row = 0; row < N_SveloDOF; row++)
      {
        if(row >= N_Sactive)
          continue;
        int interface_dof = global_DOF_interface[row][0].second;
        if(interface_dof != -1)
        {
          const double nx = normal(row).nx;
          const double tx = normal(row).tx;
          // set entries such that larger numbers are on the diagonal
          const bool nt = nx*nx - tx*tx > -1e-12; // (nx * nx >= tx * tx);
          // in case of a kink we only enforce two conditions in the two normal 
          // directions. No tangential conditions are imposed in this dof
          if(getIcond() == 1)
          {
            // u.t = 0 and u.n = eta
            (*etaToBd)(row + (nt ? 0 : N_SveloDOF), interface_dof) = 1.0;
            // set rhs to zero here, so that only eta will be there in normal
            // direction, and zero in tangential direction
            this->get_rhs()[row             ] = 0.0;
            this->get_rhs()[row + N_SveloDOF] = 0.0;
          }
          else
          {
            // Beavers-Joseph-Saffman
            (*etaToBd)(row + (nt ? 0 : N_SveloDOF), interface_dof) = 1.0;
            // set normal component of rhs to zero here, so that only eta 
            // will be there
            this->get_rhs()[row + (nt ? 0 : N_SveloDOF)] = 0.0;
          }
          /*else // kink
          {
            (*etaToBd)(row             , interface_dof) = 1.0;
            (*etaToBd)(row + N_SveloDOF, interface_dof) = 1.0;
            // set rhs to zero here, so that only eta will be there
            this->get_rhs()[row             ] = 0.0;
            this->get_rhs()[row + N_SveloDOF] = 0.0;
          }*/
        }
      }
      return; // no more assembling 
    }
    
    // assemble (eta, v.n)_Gamma here
    
    // Id of the Stokes space
    const int S_id = this->get_velocity_space().GetCollection()->GetCell(0)
        ->GetReference_ID();
    
    // everything needed for local assembling
    local_edge_assembling l;
    
    for(unsigned int iEdge = 0; iEdge < interface.size(); iEdge++)
    {
      const TInnerInterfaceJoint * thisEdge = interface[iEdge];
      TBaseCell *s_cell = thisEdge->GetNeighbour(0);
      if(s_cell->GetReference_ID() != S_id)
        s_cell = thisEdge->GetNeighbour(1);
      
      FE2D uFEId = this->get_velocity_space().GetFE2D(0, s_cell);
      TFE2D *uFE = TFEDatabase2D::GetFE2D(uFEId);
      const int N_uBf = uFE->GetN_DOF();
      int *uDOF = 
        this->get_velocity_space().GetGlobalDOF(s_cell->GetCellIndex());
      // Basis functions for Stokes velocity
      TBaseFunct2D * uBf = uFE->GetBaseFunct2D();
      //RefTrans2D refTransS = uFE->GetRefTransID();
      
      FE2D pFEId = this->get_pressure_space().GetFE2D(0, s_cell);
      TFE2D *pFE = TFEDatabase2D::GetFE2D(pFEId);
      const int N_pBf = pFE->GetN_DOF();
      int *pDOF = 
         this->get_pressure_space().GetGlobalDOF(s_cell->GetCellIndex());
      // Basis functions for Stokes pressure
      TBaseFunct2D * pBf = pFE->GetBaseFunct2D();
      
      RefTrans2D refTrans = uFE->GetRefTransID();
      TFEDatabase2D::SetCellForRefTrans(s_cell, refTrans);
      
      // get the quadrature formula
      QuadFormula1D QFId; // quadrature formula id
      int N_LinePoints; // number of qudrature points
      double *LineWeights, *zeta; // quadrature weights and points (on [-1,1])
      {
        // polynomial degree of finite element, needed for choosing an 
        // appropriate quadrature formula
        const int fe_degree = uBf->GetPolynomialDegree()
                              * abs(eta.getSpaceType());
        QFId = TFEDatabase2D::GetQFLineFromDegree(fe_degree);
        TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(QFId);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
        // qf1 no longer needed, only local scope
      }
      
      // make sure all functions & derivatives are available for this quadrature
      uBf->MakeRefElementData(QFId);
      pBf->MakeRefElementData(QFId);
      
      // compute length of the edge
      const double hE = thisEdge->GetLength();
      // index of this edge in the adjacent cell in Stokes subdomain
      const int eI = thisEdge->GetIndexInNeighbor(s_cell); // edge Index
      // normal (pointing out of the Stokes subdomain) and tangential 
      double nx, ny;
      getNormal(s_cell, thisEdge, nx, ny);
      
      bool nt = nx*nx - ny*ny > -1e-12; //(nx*nx >= ny*ny);
      // values of all functions and derivatives on reference element 
      double **uref = TFEDatabase2D::GetJointValues2D(uBf->GetID(), QFId, eI);
      double **uxiref = TFEDatabase2D::GetJointDerivatives2D(uBf->GetID(), QFId,
                                                             eI, D10);
      double **uetaref = TFEDatabase2D::GetJointDerivatives2D(uBf->GetID(),
                                                              QFId, eI, D01);
      double **pref = TFEDatabase2D::GetJointValues2D(pBf->GetID(), QFId, eI);
      
      // values and derivatives of all functions at one quadrature point
      double *uorig = new double[N_uBf];
      double *uxorig = new double[N_uBf];
      double *uyorig = new double[N_uBf];
      double *porig = new double[N_pBf];
      
      l.nt.set(nx, ny); // tangent vector not needed
      l.hE = hE;
      
      for(int k = 0; k < N_LinePoints; k++)
      {
        // quadrature weight and determinant of tranformation
        l.qw = LineWeights[k] * (hE / 2);
        
        TFEDatabase2D::GetOrigValues(refTrans, zeta[k], uBf, eI, uref[k],
                                     uxiref[k], uetaref[k], uorig, uxorig,
                                     uyorig);
        TFEDatabase2D::GetOrigValues(refTrans, zeta[k], pBf, eI, pref[k], NULL,
                                     NULL, porig, NULL, NULL);
        
        for(int row = 0; row < N_uBf; row++)
        {
          if(uDOF[row] >= N_Sactive)
            continue;
          if(kink(uDOF[row]))
            nt = normal(uDOF[row]).nx * normal(uDOF[row]).nx
                 - normal(uDOF[row]).ny * normal(uDOF[row]).ny
                 > -1e-12;
          else
            nt = nx*nx - ny*ny > -1e-12; //(nx*nx >= ny*ny);
          l.testDOF = uDOF[row] + (nt ? 0 : N_SveloDOF);
          l.at.setTest(uorig[row], uxorig[row], uyorig[row]);
          //l.nt = normal(uDOF[row]); // exists only for dofs on the interface
          
          for(int col = 0; col < N_uBf; col++)
          {
            l.ansatzDOF = global_DOF_interface[uDOF[col]][0].second;
            if(l.ansatzDOF == -1)
              continue; // not an interface dof
            if(global_DOF_interface[uDOF[col]][0].first != s_cell)
            {
              // this dof belongs to two interface dofs (this happens for
              // discontinuous interface functions), take the correct one.
              l.ansatzDOF = global_DOF_interface[uDOF[col]].at(1).second;
              if(global_DOF_interface[uDOF[col]][1].first != s_cell)
                ErrThrow("wrong stokes cell ",
                         global_DOF_interface[uDOF[col]].size());
            }
            // derivatives not needed
            l.at.setAnsatz(uorig[col], uxorig[col], uyorig[col]);
            
            localAssembleEtaToBd_velocity(l);
          }
        }
        
        for(int row = 0; row < N_pBf; row++)
        {
          l.testDOF = pDOF[row];
          l.at.setTest(porig[row], 1e10, 1e10); // derivatives not needed
          for(int col = 0; col < N_uBf; col++)
          {
            l.ansatzDOF = global_DOF_interface[uDOF[col]][0].second;
            if(l.ansatzDOF == -1)
              continue; // not an interface dof
            if(global_DOF_interface[uDOF[col]][0].first != s_cell)
              // this dof belongs to two interface dofs (this happens for
              // discontinuous interface functions), take the correct one.
              l.ansatzDOF = global_DOF_interface[uDOF[col]][1].second;
            
            // derivatives not needed
            l.at.setAnsatz(uorig[col], uxorig[col], uyorig[col]);
                           
            localAssembleEtaToBd_pressure(l);
          }
        }
      }
      delete[] uorig;
      delete[] uxorig;
      delete[] uyorig;
      delete[] porig;
    }
    etaToBd->add(l.m.m[0]);
  }
  if(typeOf_bci == Neumann) // we want to impose -nTn = eta rather than nTn=eta
    *etaToBd *= -1.0;
}

/** ************************************************************************ */
void StokesProblem::addEta(const InterfaceFunction& eta, double factor)
{
  if(etaToBd == NULL)
  {
    // etaToBd should have been done already
    ErrMsg("etaToBd has not been assembled in a " << typeOf_bci << 
           " (Navier-) Stokes problem. It will be done now.\n");
    create_etaToBd(eta);
  }
  double const * const entries = etaToBd->GetEntries();
  int const * const rows = etaToBd->GetRowPtr();
  int const * const cols = etaToBd->GetKCol();
  const int n_rows = etaToBd->GetN_Rows();
  
  for(int i = 0; i < n_rows; i++)
  {
    double value = 0;
    for(int j = rows[i]; j < rows[i + 1]; j++)
      value += entries[j] * eta(cols[j]);
    this->get_rhs()[i] += factor * value;
  }
}

/** ************************************************************************ */
void StokesProblem::copy_rhs()
{
  DrhsNSE = this->get_rhs();
}

/** ************************************************************************ */
void StokesProblem::solve(const InterfaceFunction &eta)
{
  // copy right hand side without interface integrals
  this->get_rhs() = DrhsNSE;
  // copy old solution for damping later
  this->solution_old = this->get_solution();
  
  //this->NSE2D::get_rhs().print("NSE_rhs_before");
  addEta(eta);
  //this->NSE2D::get_rhs().print("NSE_rhs_after");
  
  const int Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  // check for max number of nonlinear iterations.
  // if it is 0, solve the linear stokes problem.
  if(Max_It == 0)
  {
    /* =========================   Solve Stokes   ========================== */
    Output::print<1>("Solving a ", typeOf_bci, " Stokes problem");
    // the projection of the pressure will be done after the solve, 
    TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
    
    //this->NSE2D::solve();
    UmfpackSolverWrapper usw(this->get_matrix().get_combined_matrix().get());
    usw.solve(this->get_rhs().get_entries(), 
              this->get_solution().get_entries());
    
    if(typeOf_bci == Dirichlet)
    {
      // possibly project pressure. This is necessary if all outer boundary
      // conditions in the Stokes subdomain are of Dirichlet type
      switch(TDatabase::ParamDB->EXAMPLE)
      {
        case 0:
          this->get_pressure().project_into_L20();
          break;
        case 1:
        {
          double new_mean = 2 / TDatabase::ParamDB->RE_NR
                            + 1 / (3 * TDatabase::ParamDB->SIGMA_PERM);
          this->get_pressure().project_into_L20(new_mean);
          break;
        }
        default: // no projection for other examples
          break;
      }
    }
  }
  else
  {
    /* ======================   Solve Navier-Stokes   ====================== */
    Output::print<1>("Solving a ", typeOf_bci, " Navier-Stokes problem");
    //print_assemble_input_output();
    
    ErrThrow("Navier-Stokes does not work currently");
    
    this->NSE2D::normOfResidual();
    // nonlinear iteration
    if(!TDatabase::ParamDB->StoDa_periodicBoundary)
    {
      for(int k = 0;; k++)
      {
        if(TDatabase::ParamDB->SC_VERBOSE > 0)
          OutPut("  NONLINEAR ITERATION " << k << endl);
        
        this->assemble_nonlinear_term();
        //this->add_linear_matrix();
        if(this->stopIt(k))
          break;
        this->NSE2D::solve();
      }
    }
    else
    {
      // unfortunately for periodic boundary conditions (the riverbed example)
      // everything is much more difficult (ugly even)
      TFEFunction2D *feFunctions[2] = { this->get_velocity_component(0), 
                                       this->get_velocity_component(1) };
      for(int k = 0;; k++)
      {
        if(TDatabase::ParamDB->SC_VERBOSE > 0)
          OutPut("NONLINEAR ITERATION " << k << endl);
        
        BoundValueFunct2D *BV[3] = {example.get_bd(0), example.get_bd(1),
                                    example.get_bd(2)};
        BlockMatrixNSE2D itMat(this->get_velocity_space(), 
                               this->get_pressure_space(), BV);
        
        { // assemble nonlinear term into 'itMat'
          if(TDatabase::ParamDB->SC_VERBOSE > 1)
            OutPut(" assemble nonlinear term " << endl);
          LocalAssembling2D la(LocalAssembling2D_type::NSE2D_Galerkin_Nonlinear,
                               feFunctions, get_example().get_coeffs());
          itMat.Assemble(la, this->get_rhs());
          // remove ones on diagonal for Dirichlet rows
          //itMat.SetNonactiveDiagonal(0, 0.0);
          //itMat.SetNonactiveDiagonal(3, 0.0);
        }
        if(TDatabase::ParamDB->StoDa_periodicBoundary)
        {
          this->makePeriodicBoundary(itMat.block(0), false, true);
          this->makePeriodicBoundary(itMat.block(1), false, true);
          this->makePeriodicBoundary(itMat.block(2), false, true); // BT1
          this->makePeriodicBoundary(itMat.block(3), false, true);
          this->makePeriodicBoundary(itMat.block(4), false, true);
          this->makePeriodicBoundary(itMat.block(5), false, true); // B2T
        }
        
        if(this->stopIt(k))
          // && (k > 0))
          break;
        
        typedef preconditioner <BlockMatrix, BlockVector> block_prec;
        typedef solver <BlockMatrix, BlockVector, block_prec> block_solver;
        
        // preconditioner object
        block_prec M(&itMat);
        // solver object
        block_solver s(&itMat, &this->get_solution(), &this->get_rhs(), &M);
        
        // solve
        s.solve();
      }
    }
  }
  if(TDatabase::ParamDB->StoDa_periodicBoundary)
    checkPeriodicDOFs();
}

/** ************************************************************************ */
void StokesProblem::findPeriodicDOFs()
{
  // threshold for two doubles to be equal
  const double eps = 1e-8;
  
  // left and right boundary as periodic boundary
  const double x_left = 0.0;
  const double x_right = 2.0;
  if(TDatabase::ParamDB->EXAMPLE != 2)
    ErrThrow("only the riverbed example is supported for periodic boundaries");
  
  const TFESpace2D * fespace = &this->get_velocity_space();
  const TCollection* coll = fespace->GetCollection();
  
  const int n_dofs = fespace->GetN_DegreesOfFreedom();
  
  // first loop over cells
  for(int i1_cell = 0, N_Cells = coll->GetN_Cells(); i1_cell < N_Cells;
      i1_cell++)
  {
    TBaseCell *cell1 = coll->GetCell(i1_cell);
    // check if face on boundary
    for(int j1 = 0, N_edges = cell1->GetN_Edges(); j1 < N_edges; j1++)
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
      
      // check if vertex x = x_left
      bool left = false;
      if(fabs(x11 - x_left) > eps || fabs(x12 - x_left) > eps)
      { // one vertex does not lie on the left boundary
        if(fabs(x11 - x_right) > eps || fabs(x12 - x_right) > eps)
          continue; // one vertex does not lie on the right boundary
      }
      else
        left = true;
      
      const FE2D FEid1 = fespace->GetFE2D(i1_cell, cell1);
      const TFE2D * FE1 = TFEDatabase2D::GetFE2D(FEid1);
      
      // global indices of all degrees of freedom in this cell 
      const int *globalDOF1 = fespace->GetGlobalDOF(i1_cell);
      // local degrees of freedom which correspond to this edge
      const int* localDOF1 = FE1->GetFEDesc2D()->GetJointDOF(j1);
      const int N_localDOF1 = FE1->GetFEDesc2D()->GetN_JointDOF();
      
      // midpoint of edge
      const double y1 = (y11 + y12) / 2;
      
      // find cell which should be coupled to this cell
      // inner loop over the cells
      for(int i2_cell = i1_cell + 1; i2_cell < N_Cells; i2_cell++)
      {
        TBaseCell *cell2 = coll->GetCell(i2_cell);
        // check if face on boundary
        for(int j2 = 0, N_edges2 = cell2->GetN_Edges(); j2 < N_edges2; j2++)
        {
          const TJoint *joint2 = cell2->GetJoint(j2);
          //not on boundary
          if(joint2->GetType() != BoundaryEdge)
            continue;
          double x21, x22, y21, y22;
          cell2->GetVertex(j2)->GetCoords(x21, y21);
          cell2->GetVertex((j2 + 1) % n_Vert)->GetCoords(x22, y22);
          
          if(fabs(x21 - (left ? x_right : x_left)) > eps
             || fabs(x22 - (left ? x_right : x_left)) > eps)
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
          
          if(TDatabase::ParamDB->SC_VERBOSE > 2)
            OutPut(
                " creating a vertical periodic boundary at y=(" << y21 << "," << y22 << ")\n");
          for(int edge_dof = 0; edge_dof < N_localDOF1; edge_dof++)
          {
            // due to counterclockwise numbering in each cell we have to go 
            // through one edge the opposite way:
            const int dof1 = globalDOF1[localDOF1[N_localDOF1 - 1 - edge_dof]];
            const int dof2 = globalDOF2[localDOF2[edge_dof]];
            
            if(left)
            {
              periodic_dofs[dof1] = dof2;
              periodic_dofs[dof1 + n_dofs] = dof2 + n_dofs; // second component
              //OutPut(" dofs " << dof1 << "\t" << dof2 << endl);
              //OutPut("      " << dof1+n_dofs << "\t" << dof2+n_dofs << endl)
            }
            else
            {
              periodic_dofs[dof2] = dof1;
              periodic_dofs[dof2 + n_dofs] = dof1 + n_dofs; // second component
              //OutPut(" dofs " << dof2 << "\t" << dof1 << endl);
              //OutPut("      " << dof2+n_dofs << "\t" << dof1+n_dofs << endl)
            }
          }
        }
      }
    }
  }
  OutPut(
      "There are " << periodic_dofs.size()/2 << " periodic Stokes degrees of freedom\n");
}

/** ************************************************************************ */
void StokesProblem::makePeriodicBoundary(std::shared_ptr<TMatrix> mat,
                                         bool stokesMat, bool p )
{
  if(periodic_dofs.empty())
    ErrThrow("called StokesProblem::makePeriodicBoundary with map ",
             "'periodic_dofs' not yet set");
  if(!mat)
  {
    makePeriodicBoundary(this->get_matrix().block(0), true, true);  // A11
    makePeriodicBoundary(this->get_matrix().block(1), false, true);
    makePeriodicBoundary(this->get_matrix().block(2), false, true); // B!t
    makePeriodicBoundary(this->get_matrix().block(3), false, true);
    makePeriodicBoundary(this->get_matrix().block(4), true, true);  // A22
    makePeriodicBoundary(this->get_matrix().block(5), false, true); // B2T
    return;
  }
  
  // assume all matrices A.. have the same structure
  int const * const rowPtr = mat->GetRowPtr();
  int const * const colPtr = mat->GetKCol();
  double const * const entries = mat->GetEntries();
  const int n_rows = mat->GetN_Rows();
  // this map will be passed to "mat->changeRows(new_rows)" 
  std::map<int, std::map<int, double> > new_rows;
  
  for(std::map<int, int>::iterator ii = periodic_dofs.begin();
      ii != periodic_dofs.end() && ii->first < n_rows; ++ii)
  {
    //OutPut( ii->first << ": " << ii->second << endl);
    // the row with number "ii->second" is added to the row number "ii->first". 
    // Then row number "ii->second" is replaced by a row with two entries, 1 
    // on the diagonal and -1 on the "ii->first"-th entry. Also the right hand 
    // side is changed to be zero
    
    // loop over all entries in row "ii->first" of A-matrices
    for(int i = rowPtr[ii->first]; i < rowPtr[1 + ii->first]; i++)
    {
      if(entries[i] != 0.0 || p)
        (new_rows[ii->first])[colPtr[i]] += entries[i];
      
    } // up to here this row would simply be copied.
    // loop over all entries in row "ii->second" of A-matrices
    for(int i = rowPtr[ii->second]; i < rowPtr[1 + ii->second]; i++)
    {
      if(entries[i] != 0.0 || p)
        (new_rows[ii->first])[colPtr[i]] += entries[i];
    }
    if(stokesMat)
    {
      (new_rows[ii->second])[ii->second] = 1.0; // diagonal
      (new_rows[ii->second])[ii->first] = -1.0; // coupling
      // here two entries for the right hand side are handled.
      // One would be enough, but which one depends on which matrix this is 
      // (A11 or A22).
      this->get_rhs()[ii->second] = 0;
      this->get_rhs()[ii->second + this->get_rhs().length(0)] = 0;
    }
    else if(!p || n_rows != mat->GetN_Columns())
      (new_rows[ii->second]); // set entire row to zero
    else
    {
      (new_rows[ii->second])[ii->second] = 0.0; // diagonal
      (new_rows[ii->second])[ii->first] = 0.0; // coupling
    }
  }
  mat->changeRows(new_rows);
}

/** ************************************************************************ */
void StokesProblem::checkPeriodicDOFs()
{
  std::map<int, int>::iterator it;
  double const * const vals = this->get_solution().get_entries();
  for(it = periodic_dofs.begin(); it != periodic_dofs.end(); ++it)
  {
    double diff = vals[it->first] - vals[it->second];
    if(fabs(diff) > 1e-6)
      Error("WRONG Stokes here " << it->first << "\t" << it->second << "\t" << 
            setprecision(10) << diff << endl);
  }
}

/** ************************************************************************ */
n_t_vector& StokesProblem::normal(unsigned int dof, unsigned int i)
{
  if(i != 0 && i != 1 && i != 10)
    ErrThrow("'i' must be 0, 1, or 10 (which means let me choose)");
  if(dof >= normals.size())
    ErrThrow("trying to get normal of dof ", dof, " which does not exist");
  // change to the correct i, in case i==10 (default)
  if(i == 10)
    i = kink(dof) ? 2 : 0;
  if( i == 2 && normals[dof].size() == 2) // we have a kink, if i==2
  {
    // create auxilary normal at kink, then normals[dof].size()==3
    const double det = normal(dof, 0).nx * normal(dof, 1).ny
                         - normal(dof, 0).ny * normal(dof, 1).nx;
    const double nx = (normal(dof, 1).ny - normal(dof, 0).ny) / det;
    const double ny = (normal(dof, 0).nx - normal(dof, 1).nx) / det;
    normals[dof].push_back(n_t_vector(nx, ny));
  }
  if(i >= normals[dof].size())
    ErrThrow("there are only ", normals[dof].size(), " normal to this dof ",
             dof);
  return normals[dof][i];
}

/** ************************************************************************ */
void StokesProblem::print_normals()
{
  OutPut("printing normals ...\n");
  
  const int ID = 
    this->get_velocity_space().GetCollection()->GetCell(0)->GetReference_ID();
  const int vActiveBound = this->get_velocity_space().GetActiveBound();
  
  int n_vdof = this->get_velocity_space().GetN_DegreesOfFreedom();
  
  OutPut("number of Stokes velocity dofs " << n_vdof << endl);
  
  for(unsigned int j = 0; j < interface.size(); j++)
  {
    TBaseCell *cell = interface[j]->GetNeighbour(0);
    if(cell->GetReference_ID() != ID)
      cell = interface[j]->GetNeighbour(1);
    // index of this edge in this cell (0,1, or 2 for triangles)
    
    const int cell_index = cell->GetCellIndex();
    int const * const vDOF = 
      this->get_velocity_space().GetGlobalDOF(cell_index);
    FE2D vFEId = this->get_velocity_space().GetFE2D(0, cell);
    TBaseFunct2D *vBaseFunct = TFEDatabase2D::GetBaseFunct2DFromFE2D(vFEId);
    const int N_vBaseFunct = vBaseFunct->GetDimension();
    
    for(int lvdof = 0; lvdof < N_vBaseFunct; lvdof++) // local velocity dof
    {
      const int this_vdof = vDOF[lvdof];
      if(global_DOF_interface[this_vdof][0].second == -1)
      {
        // not an interface dof
        OutPut("interface " << j << " cell " << cell_index << " dof " << 
               this_vdof << endl);
      }
      else
      {
        OutPut("interface " << j << " cell " << cell_index << " dof " << 
               this_vdof << " ");        
        if(kink(this_vdof))
        {
          OutPut("kink!\n" << "\t" << normal(this_vdof, 0) << "\n\t" << 
                 normal(this_vdof, 1) << "\n\t" << normal(this_vdof));
        }
        else
        {
          OutPut(normal(this_vdof));
        }
        OutPut("\t" << global_DOF_interface[this_vdof].size() << " " << 
               global_DOF_interface[this_vdof][0].second);
        
        if(this_vdof >= vActiveBound)
        {
          OutPut(" Dirichlet " << endl);
        }
        else
          OutPut(endl);
      }
    }
  }
}
