#include <algorithm>
#include <DarcyPrimal.h>
#include <FEDatabase2D.h>

/** ************************************************************************ */
DarcyPrimal::DarcyPrimal(const TDomain& domain, const Example_CD2D& ex, 
                         InterfaceCondition t,
                         std::vector<const TInnerInterfaceJoint* >& in,
                         int reference_id) :
    CD2D(domain, ex, reference_id),
    typeOf_bci(t), IrhsDARCY(this->get_rhs()), 
    // actually the old solution is only needed if 
    // TDatabase::ParamDB->StoDa_solutionStrategy != 0, but we do this anyway
    solution_old(this->get_solution()), 
    fe_solution_old(&this->get_space(), (char*)"Darcy_p_old", 
                    (char*)"Darcy_p_old", solution_old.get_entries(),
                    solution_old.length()),
    // actually the direct solution is only needed if 
    // TDatabase::ParamDB->StoDa_solutionStrategy > 0, but we do this anyway
    solution_direct(this->get_solution()),
    fe_solution_direct(&this->get_space(), (char*) "Darcy_p_direct", 
                       (char*) "Darcy_p_direct", solution_direct.get_entries(),
                       solution_direct.length()),
    interface(in), stokesToDarcy(nullptr),
    mat_aux(nullptr), etaToBd(nullptr), map_sol2eta(nullptr), 
    map_sol2eta_rhs(), periodic_dofs(), eta_hom(nullptr), 
    global_DOF_interface(), assemble_on_input(false), assemble_on_return(true)
{
  
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
        assemble_on_input = stokes_first ? false : true;
        assemble_on_return = stokes_first ? false : true;
        break;
      case Dirichlet:
        assemble_on_input = true;
        assemble_on_return = true;
        // note that here it should be stokes_first == true, because in the
        // other case we don't need a Dirichlet Darcy problem
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
DarcyPrimal::~DarcyPrimal()
{
  //Output::print<1>("DarcyPrimal destructor");
  delete eta_hom;
}

/** ************************************************************************ */
void DarcyPrimal::copy_rhs()
{
  IrhsDARCY = this->get_rhs();
}

/** ************************************************************************ */
void DarcyPrimal::mat_DARCY_iIntegrals(BlockMatrix * m)
{
  if(m == NULL)
    m = &this->get_matrix();
  Output::print<1>(" Assemble interface integrals into matrix (Darcy)");
  
  // velocity-velocity
  TSquareMatrix2D* Amat = this->get_matrix().get_matrix();
  
  const TFESpace2D & pSpace = this->get_space();
  if(typeOf_bci == Dirichlet && TDatabase::ParamDB->StoDa_weakGamma <= 0.0)
  {
    // set Dirichlet dofs exactly, instead of using weak Dirichlet conditions)
    // this map will be used to call Amat->changeRows(new_rows, true) later on
    std::map<int, std::map<int, double> > new_rows;
    const int n_dof = pSpace.GetN_DegreesOfFreedom();
    for(int row = 0; row < n_dof; row++)
    {
      if(global_DOF_interface[row][0].second == -1)
        continue; // not an interface dof
      new_rows[row][row] = 1.0; // set diagonal entry only, all others are zero
    }
    Amat->changeRows(new_rows); // this is quite expensive
    return; // no more assembling here
  }
  
  const int ActiveBound = pSpace.GetActiveBound();
  const int ID = pSpace.GetCollection()->GetCell(0)->GetReference_ID();
  local_edge_assembling l;
  
  for(unsigned int j = 0; j < interface.size(); j++)
  {
    TBaseCell *cell = interface[j]->GetNeighbour(0);
    if(cell->GetReference_ID() != ID)
      cell = interface[j]->GetNeighbour(1);
    FE2D pFEId = pSpace.GetFE2D(0, cell);
    TFE2D *pFE = TFEDatabase2D::GetFE2D(pFEId);
    const int N_pBf = pFE->GetN_DOF(); // number of Darcy pressure basis functs
    int *pDOF = pSpace.GetGlobalDOF(cell->GetCellIndex());
    TBaseFunct2D *pBf = pFE->GetBaseFunct2D();
    
    RefTrans2D refTrans = pFE->GetRefTransID();
    TFEDatabase2D::SetCellForRefTrans(cell, refTrans);
    
    int fe_degree = 2 * pBf->GetPolynomialDegree();
    
    QuadFormula1D QFId = TFEDatabase2D::GetQFLineFromDegree(fe_degree);
    int N_LinePoints; // number of qudrature points
    double *LineWeights, *zeta; // quadrature weights and points (on [-1,1])
    {
      TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(QFId);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
    }
    
    pBf->MakeRefElementData(QFId);
    
    // index of this edge in this cell (0,1, or 2 for triangles)
    const int eI = interface[j]->GetIndexInNeighbor(cell);
    
    // compute length of the edge
    const double hE = interface[j]->GetLength();
    // tangential/normal vectors to this boundary (normalized)
    double nx, ny;
    getNormal(cell, interface[j], nx, ny);
    
    double **pref = TFEDatabase2D::GetJointValues2D(pBf->GetID(), QFId, eI);
    double **pxiref = TFEDatabase2D::GetJointDerivatives2D(pBf->GetID(), QFId,
                                                           eI, D10);
    double **petaref = TFEDatabase2D::GetJointDerivatives2D(pBf->GetID(), QFId,
                                                            eI, D01);
    
    double *porig = new double[N_pBf];
    double *pxorig = new double[N_pBf];
    double *pyorig = new double[N_pBf];
    
    l.nt.set(nx, ny, 0, 0);
    l.hE = hE;
    
    // loop of quadrature points
    for(int k = 0; k < N_LinePoints; k++)
    {
      TFEDatabase2D::GetOrigValues(refTrans, zeta[k], pBf, eI, pref[k],
                                   pxiref[k], petaref[k], porig, pxorig,
                                   pyorig);
      // quadrature weight and determinant of transformation
      l.qw = LineWeights[k] * (hE / 2);
      for(int lt = 0; lt < N_pBf; lt++)
      {
        l.testDOF = pDOF[lt];
        if(l.testDOF >= ActiveBound)
          continue; // do not change non-active rows
        l.at.setTest(porig[lt], pxorig[lt], pyorig[lt]);
        
        for(int l2 = 0; l2 < N_pBf; l2++)
        {
          l.ansatzDOF = pDOF[l2];
          l.at.setAnsatz(porig[l2], pxorig[l2], pyorig[l2]);
          
          localAssembleMat_pressure_pressure(l);
        } //for(int l2=0; l2<N_pBaseFunct; l2++)
      } //for(int l=0;l<N_vBaseFunct;l++)
    } //for(int k=0;k<N_LinePoints;k++)
    delete[] porig;
    delete[] pxorig;
    delete[] pyorig;
    
  } //for(int j=0; j<N_InnerInterfaceJoints; j++)
  
  Amat->add(l.m.m[0]);
}

/** ************************************************************************ */
void DarcyPrimal::localAssembleMat_pressure_pressure(local_edge_assembling &l)
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  const double nx = l.nt.nx, ny = l.nt.ny; // =n_p = -n_f
  const double p = l.at.u, px = l.at.ux, py = l.at.uy;
  const double q = l.at.v, qx = l.at.vx, qy = l.at.vy;
  const int tDOF = l.testDOF;
  const int aDOF = l.ansatzDOF;
  const double qw = l.qw;
  switch(typeOf_bci)
  {
    case Neumann:
    {
      // no further assembling needed
      break;
    }
    case Robin:
    {
      // (phi,psi) (1/gamma_p)
      l.m.m[0][tDOF][aDOF] += p * q * qw / (TDatabase::ParamDB->StoDa_gamma_p);
      break;
    }
    case weakRobin:
    {
      // -(K (grad phi).n,psi)
      l.m.m[0][tDOF][aDOF] += -K * (px * nx + py * ny) * q * qw;
      // (gamma_p(K (grad phi).n) + phi, 
      //   psi - gamma h_E K (grad psi).n)  *1/(gamma_f+gamma h_E)
      const double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma);
      const double gamma_p = TDatabase::ParamDB->StoDa_gamma_p;
      const double tilde_lambda = 1 / (gamma_p + gamma * l.hE);
      // note normal is pointing out of Darcy subdomain, into Stokes 
      double val2 = gamma_p * K * (px * nx + py * ny) + p; // ansatz fuction
      val2 *= q + gamma * l.hE * K * (qx * nx + qy * ny); // test function
      l.m.m[0][tDOF][aDOF] += val2 * tilde_lambda * qw;
      
      break;
    }
    case DirichletSTAB:
    {
      // (phi,psi)* (gamma_p K /h)
      l.m.m[0][tDOF][aDOF] += p * q * TDatabase::ParamDB->StoDa_gamma_p * K * qw
                              / l.hE;
      // ((K (grad phi).n)-eta,psi)*(gamma_p/(gamma_p+h))
      double val2 = -K * (px * nx + py * ny) * q * qw
                    * TDatabase::ParamDB->StoDa_gamma_p
                    / (TDatabase::ParamDB->StoDa_gamma_p + l.hE);
      l.m.m[0][tDOF][aDOF] += val2;
      break;
    }
    case Dirichlet:
    {
      const double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma);
      // (phi,psi)* (gamma/h)
      l.m.m[0][tDOF][aDOF] += p * q * gamma * qw / l.hE;
      // -((K (grad phi).n),psi)
      l.m.m[0][tDOF][aDOF] -= K * (px * nx + py * ny) * q * qw;
      l.m.m[0][tDOF][aDOF] -= K * (qx * nx + qy * ny) * p * qw; // symmetric
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
void DarcyPrimal::localAssembleEtaToBd_pressure(local_edge_assembling &l)
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  const double nx = l.nt.nx, ny = l.nt.ny;
  const double xi = l.at.u; //, px = l.at.ux, py = l.at.uy;
  const double q = l.at.v, qx = l.at.vx, qy = l.at.vy;
  const int tDOF = l.testDOF;
  const int aDOF = l.ansatzDOF;
  const double qw = l.qw;
  
  switch(typeOf_bci)
  {
    case Neumann:
    {
      // (eta,psi)
      l.m.m[0][tDOF][aDOF] += xi * q * qw;
      break;
    }
    case Robin:
    {
      // (eta,psi) / gamma_p
      l.m.m[0][tDOF][aDOF] += xi * q * qw / TDatabase::ParamDB->StoDa_gamma_p;
      break;
    }
    case weakRobin:
    {
      // (eta, psi - gamma h_E K (grad psi).n) * 1/(gamma_p+gamma h_E)
      const double gamma_p = TDatabase::ParamDB->StoDa_gamma_p;
      const double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma);
      l.m.m[0][tDOF][aDOF] += xi * (q - gamma * l.hE * K * (qx * nx + qy * ny))
                                 * qw / (gamma_p + gamma * l.hE);
      break;
    }
    case DirichletSTAB:
    {
      ErrThrow("Problem type DirichletSTAB does not work yet");
      break;
    }
    case Dirichlet:
    {
      const double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma);
      // (eta,psi) * (gamma/h)
      l.m.m[0][tDOF][aDOF] += xi * q * gamma * qw / l.hE;
      l.m.m[0][tDOF][aDOF] -= xi * K * (qx * nx + qy * ny) * qw; // symmetric
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
void DarcyPrimal::prepare_map_solution_to_interface(
    const InterfaceFunction &eta)
{
  // the map which is used to map a computed solution to an interface function
  // will be affine linear (even linear for Neumann problems).

  /** create the structure of the underlying linear map (matrix) */
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
  else if(typeOf_bci == Dirichlet)
    // || typeOf_bci == Robin || typeOf_bci == weakRobin)
  {
    map_sol2eta_rhs.copy_structure((const BlockVector&)eta);
    // assemble an integral in (parts of) the domain
    assemble_Neumann_map_to_interface(m);
    
  }
  else if(typeOf_bci == Robin || typeOf_bci == weakRobin)
  {
    map_sol2eta_rhs.copy_structure((const BlockVector&)eta);
    // D_RR should be false here
    // return -gamma_f * flux - pressure 
    assemble_Dirichlet_map_to_interface(m, D_RR, -1.0);
    // no minus here due to different normal
    assemble_Neumann_map_to_interface(m, TDatabase::ParamDB->StoDa_gamma_f);
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
void DarcyPrimal::find_global_DOF_interface(const InterfaceFunction &eta)
{
  if(this->global_DOF_interface.size() != 0)
    return; // no need to call this a second time
  
  const TFESpace2D& pSpace = this->get_space();
  
  // either do (1) or do not (-1) invert the direction
  // To each interface edge there is a direction, i.e. a vector pointing from
  // the first point to the second. This direction has to coincide with the 
  // direction of this edge in the neighboring cell. In there the direction is 
  // counterclockwise. 
  // Each edge direction can only coincide with either the neighboring Stokes
  // or the neighboring Darcy cell.
  // if invertDirection is -1 then the reference transformations F_D for the 
  // Darcy dofs and F_I the interface dofs ((each) from [-1,1] to the edge 
  // [\vec a,\vec b]) have a different sign in front of their derivative: E.g.,
  // F_D(-1) = \vec a = F_I(1)  and  F_D(1) = \vec b = F_I(-1).
  // Therefore if invertDirection is -1 we simply use F_I(-z) instead of F_I(z)
  // to evaluate the interface function at the same quadrature point as in the 
  // Darcy part, F_D(z) = F_I(-z). This relies on the fact that the quadrature
  // points and weights are symmetric in [-1,1]
  int invertDirection = 1;
  
  this->global_DOF_interface.resize(pSpace.GetN_DegreesOfFreedom());
  //,vector<pair<TBaseCell *, int> >(1,pair<TBaseCell *, int>(NULL, -1)));
  
  // Id of the Darcy space
  const int D_id = pSpace.GetCollection()->GetCell(0)->GetReference_ID();
  
  for(unsigned int iEdge = 0; iEdge < interface.size(); iEdge++)
  {
    const TInnerInterfaceJoint * thisEdge = interface[iEdge];
    TBaseCell *d_cell = thisEdge->GetNeighbour(0);
    if(d_cell->GetReference_ID() != D_id)
      d_cell = thisEdge->GetNeighbour(1);
    
    FE2D pFEId = pSpace.GetFE2D(0, d_cell);
    TFE2D *pFE = TFEDatabase2D::GetFE2D(pFEId);
    // number of Darcy pressure basis functions
    const int N_pBf = pFE->GetN_DOF(); 
    int *pDOF = pSpace.GetGlobalDOF(d_cell->GetCellIndex());
    // Basis functions for Stokes velocity
    TBaseFunct2D * pBf = pFE->GetBaseFunct2D();
    
    // number of local degrees of freedom on the interface 
    const int N_iBaseFunct = 3; //abs(eta.spaceType)+1;
    
    // index of this edge in the adjacent cell in Darcy subdomain
    const int eI = thisEdge->GetIndexInNeighbor(d_cell); // edge Index
        
    // normal (pointing out of the Darcy subdomain) and tangential 
    double nx, ny, tx, ty;
    getNormal(d_cell, thisEdge, nx, ny);
    thisEdge->GetTangent(tx, ty);
    // check if edge is directed the same way as in the cell 'd_cell'
    if(tx == -ny && ty == nx)
      invertDirection = 1;
    else
      //if(tx==ny && ty==-nx)
      invertDirection = -1;
    
    double points_on_interface[3] = {-1.0, 0.0, 1.0}; //points on reference edge
    if(invertDirection == -1)
    {
      points_on_interface[0] = 1.0;
      points_on_interface[2] = -1.0;
    }
    double *values[N_iBaseFunct];
    for(int i = 0; i < N_iBaseFunct; i++)
      values[i] = new double[N_pBf];
    pBf->GetValues(N_iBaseFunct, points_on_interface, eI, values);
    
    // now check which of the Darcy pressure basis functions are one at a point
    // on the interface and zero elsewhere
    for(int i = 0; i < N_iBaseFunct; i++)
    {
      int iDOF = eta.getDOF(iEdge, i);
      for(int j = 0; j < N_pBf; j++)
      {
        if(values[i][j] == 1.0)
        {
          // the i-th interface dof corresponds to the j-th Darcy pressure dof
          this->global_DOF_interface[pDOF[j]].push_back(std::make_pair(d_cell,
                                                                       iDOF));
          
          //Output::print<1>("add ", pDOF[j], " ", d_cell, " ", iDOF);
          continue;
        }
      }
    }
    for(int i = 0; i < N_iBaseFunct; i++)
      delete[] values[i];
  }
  
  // find cells not sharing an edge with the interface but only a vertex and add
  // those to 'global_DOF_interface'
  TCollection * d_coll = pSpace.GetCollection();
  const int n_cells = d_coll->GetN_Cells();
  for(int icell = 0; icell < n_cells; icell++)
  {
    TBaseCell *d_cell = d_coll->GetCell(icell);
    FE2D pFEId = pSpace.GetFE2D(0, d_cell);
    TFE2D *pFE = TFEDatabase2D::GetFE2D(pFEId);
    // number of Darcy pressure basis functions
    const int N_pBf = pFE->GetN_DOF(); 
    int *pDOF = pSpace.GetGlobalDOF(d_cell->GetCellIndex());
    
    for(int j = 0; j < N_pBf; j++)
    {
      int d_dof = pDOF[j];
      // this vector of pairs
      std::vector<std::pair<TBaseCell*,int>>& tvop =global_DOF_interface[d_dof];
      // check if this dof is on the interface and if so, enter this loop,
      if(tvop.empty())
        continue;
      // check if it is associated to this cell already
      bool already_included = false;
      for(unsigned int n = 0; n < tvop.size(); n++)
      {
        if(tvop[n].first == d_cell)
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
        // anyway, because no integrals over the domain are necessary for D-RR
        tvop.push_back(std::make_pair(d_cell, iDOF));
      }
    }
  }
  
  // fill entries in global_DOF_interface which are not yet filled with defaults
  for(int i = 0; i < pSpace.GetN_DegreesOfFreedom(); i++)
  {
    if(global_DOF_interface[i].empty())
      global_DOF_interface[i].push_back(std::pair<TBaseCell *, int>(NULL, -1));
  }
}

/** ************************************************************************ */
void DarcyPrimal::create_structure_of_map_solution_to_interface(unsigned int l)
{
  const TFESpace2D& pSpace = this->get_space();
  
  // everything needed to call the constructor of TStructure
  const int n_rows = l;
  const int n_cols = IrhsDARCY.length(); // N_presDOF;
  int N_entries = 0;
  int *cols, *rows = new int[n_rows + 1];
  
  // for each row (i.e, interface dof) we have one vector containing all 
  // (Darcy) dofs which couple with this dof. For all dofs which do not 
  // belong to a cell adjacent to the interface, the corresponding column will 
  // be completely zero
  std::vector<std::vector<int> > coupledDOF(n_rows);
  
  // array containing the number of basis functions for all finite elements
  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  
  // the collection of all Darcy cells
  const TCollection* d_coll = pSpace.GetCollection();
  // number of Darcy cells
  const unsigned int n_cells = (unsigned int) d_coll->GetN_Cells();
  
  // loop over the interface edges
  for(unsigned int icell = 0; icell < n_cells; icell++)
  {
    TBaseCell *d_cell = d_coll->GetCell(icell);
    
    // pressure
    FE2D pFEId = pSpace.GetFE2D(0, d_cell);
    const int N_pBaseFunct = N_BaseFunct[pFEId];
    int *pDOF = pSpace.GetGlobalDOF(d_cell->GetCellIndex());
    
    // loop over all degrees of freedom of this Darcy cell
    for(int row = 0; row < N_pBaseFunct; row++)
    {
      int iDOF = global_DOF_interface[pDOF[row]][0].second;
      if(iDOF == -1)
        continue; // not an interface dof
      int l = 0;
      while(global_DOF_interface[pDOF[row]][l].first != d_cell)
        l++;
      // this dof belongs to more than one interface dof (this happens for
      // discontinuous interface functions), take the correct one.
      iDOF = global_DOF_interface[pDOF[row]][l].second;
      
      // loop over all pressure degrees of freedom of this Darcy cell
      for(int col = 0; col < N_pBaseFunct; col++)
      {
        const int ansatz_pDOF = pDOF[col];
        // pressure component
        coupledDOF[iDOF].push_back(ansatz_pDOF);
      }
    }
  }
  
  // fill the array 'rows' and compute 'N_entries'
  rows[0] = 0;
  for(int i = 0; i < n_rows; i++)
  {
    // sort the Darcy dofs for the i-th interface dof
    std::sort(coupledDOF[i].begin(), coupledDOF[i].end());
    // remove all duplicates 
    const std::vector<int>::iterator it = unique(coupledDOF[i].begin(),
                                                 coupledDOF[i].end());
    // resize to the size without duplicates
    coupledDOF[i].resize(distance(coupledDOF[i].begin(), it));
    // number of Darcy dofs which couple with this interface dof  
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
void DarcyPrimal::assemble_Dirichlet_map_to_interface(local_matrices & m,
                                                       bool D_RR, double a)
{
  if(D_RR)
  {
    // D-RR is somewhat difficult, so we separate it from the rest
    this->assemble_Dirichlet_map_to_interface_D_RR(m, a);
    return;
  }
  
  const TFESpace2D & pSpace = this->get_space();
  
  if(!assemble_on_return)
  {
    // instead of assembling terms on the interface, we only return the nodal
    // values and the assembling is done in the object which recieves the data
    // on the interface from this object
    if(D_RR)
      ErrMsg("assemble_on_return is false and D-RR true. This is not correct");
    const int n_dofs = pSpace.GetN_DegreesOfFreedom();
    for(int idof = 0; idof < n_dofs; idof++)
    {
      int interface_dof = global_DOF_interface[idof][0].second;
      if(interface_dof != -1)
      {
        m.m[0][interface_dof][idof] += 1.0;
      }
    }
    return;
  }
  // assemble terms 
  
  // Id of the Darcy space
  const int D_id = pSpace.GetCollection()->GetCell(0)->GetReference_ID();
  
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  
  for(unsigned int iEdge = 0; iEdge < interface.size(); iEdge++)
  {
    const TInnerInterfaceJoint * thisEdge = interface[iEdge];
    TBaseCell *d_cell = thisEdge->GetNeighbour(0);
    if(d_cell->GetReference_ID() != D_id)
      d_cell = thisEdge->GetNeighbour(1);
    
    FE2D pFEId = pSpace.GetFE2D(0, d_cell);
    TFE2D *pFE = TFEDatabase2D::GetFE2D(pFEId);
    // number of Darcy pressure basis functions
    const int N_pBf = pFE->GetN_DOF();
    int *pDOF = pSpace.GetGlobalDOF(d_cell->GetCellIndex());
    // Basis functions for Stokes velocity
    TBaseFunct2D * pBf = pFE->GetBaseFunct2D();
    
    RefTrans2D refTrans = pFE->GetRefTransID();
    TFEDatabase2D::SetCellForRefTrans(d_cell, refTrans);
    
    // polynomial degree of finite element, needed for choosing an 
    // appropriate quadrature formula, 2 is the polynomial degree of functions 
    // on the interface
    const int fe_degree = 2 * pBf->GetPolynomialDegree();
    // get the quadrature formula
    // quadrature formula id:
    QuadFormula1D QFId = TFEDatabase2D::GetQFLineFromDegree(fe_degree);
    int N_LinePoints; // number of qudrature points
    double *LineWeights, *zeta; // quadrature weights and points (on [-1,1])
    {
      TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(QFId);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
      // qf1 no longer needed, only local scope
    }
    
    // make sure all functions & derivatives are available for this quadrature
    pBf->MakeRefElementData(QFId);
    
    // compute length of the edge
    const double hE = thisEdge->GetLength();
    
    // index of this edge in the adjacent cell in Darcy subdomain
    const int eI = thisEdge->GetIndexInNeighbor(d_cell); // edge Index
    // normal (pointing out of the Darcy subdomain) and tangential 
    double nx, ny;
    if(D_RR)
      getNormal(d_cell, thisEdge, nx, ny);
    
    double **pref = TFEDatabase2D::GetJointValues2D(pBf->GetID(), QFId, eI);
    double **pxiref =
        D_RR ?
            TFEDatabase2D::GetJointDerivatives2D(pBf->GetID(), QFId, eI, D10) :
            NULL;
    double **petaref =
        D_RR ?
            TFEDatabase2D::GetJointDerivatives2D(pBf->GetID(), QFId, eI, D01) :
            NULL;
    
    // values and derivatives of all functions at one quadrature point
    double *porig = new double[N_pBf];
    double *pxorig = D_RR ? new double[N_pBf] : NULL;
    double *pyorig = D_RR ? new double[N_pBf] : NULL;
    
    for(int k = 0; k < N_LinePoints; k++)
    {
      if(D_RR)
        TFEDatabase2D::GetOrigValues(refTrans, zeta[k], pBf, eI, pref[k],
                                     pxiref[k], petaref[k], porig, pxorig,
                                     pyorig);
      else
        TFEDatabase2D::GetOrigValues(refTrans, zeta[k], pBf, eI, pref[k], NULL,
                                     NULL, porig, NULL, NULL);
      
      // quadrature weight and determinant of tranformation
      const double qw = LineWeights[k] * (hE / 2);
      
      for(int row = 0; row < N_pBf; row++)
      {
        int testDOF = global_DOF_interface[pDOF[row]][0].second;
        if(testDOF == -1)
          continue; // not an interface dof
        if(global_DOF_interface[pDOF[row]][0].first != d_cell)
          // this dof belongs to two interface dofs (this happens for
          // discontinuous interface functions), take the correct one.
          testDOF = global_DOF_interface[pDOF[row]][1].second;
        const double q = porig[row]; // test function
        
        for(int col = 0; col < N_pBf; col++)
        {
          const int ansatzDOF = pDOF[col];
          const double p = porig[col]; // ansatz function
          
          if(!D_RR)
            m.m[0][testDOF][ansatzDOF] += a * q * p * qw;
          else
          {
            m.m[0][testDOF][ansatzDOF] -= a * q * p * qw;
            m.m[0][testDOF][ansatzDOF] +=
                a * K * TDatabase::ParamDB->StoDa_gamma_f * q
                * (pxorig[col] * nx + pyorig[col] * ny) * qw;
          }
        }
      }
    }
    delete[] porig;
    delete[] pxorig;
    delete[] pyorig;
  }
}

/** ************************************************************************ */
void DarcyPrimal::assemble_Dirichlet_map_to_interface_D_RR(local_matrices & m,
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
  
  // note that this is very similar to 
  // 'DarcyPrimal::assemble_Dirichlet_map_to_interface'
  const TFESpace2D & pSpace = this->get_space();
  
  const int D_id = pSpace.GetCollection()->GetCell(0)->GetReference_ID();
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  
  for(unsigned int iEdge = 0; iEdge < interface.size(); iEdge++)
  {
    const TInnerInterfaceJoint * thisEdge = interface[iEdge];
    TBaseCell *d_cell = thisEdge->GetNeighbour(0);
    if(d_cell->GetReference_ID() != D_id)
      d_cell = thisEdge->GetNeighbour(1);
    
    FE2D pFEId = pSpace.GetFE2D(0, d_cell);
    TFE2D *pFE = TFEDatabase2D::GetFE2D(pFEId);
    // number of Darcy pressure basis functions
    const int N_pBf = pFE->GetN_DOF();
    int *pDOF = pSpace.GetGlobalDOF(d_cell->GetCellIndex());
    // Basis functions for Stokes velocity
    TBaseFunct2D * pBf = pFE->GetBaseFunct2D();
    
    RefTrans2D refTrans = pFE->GetRefTransID();
    TFEDatabase2D::SetCellForRefTrans(d_cell, refTrans);
    
    
    // get the quadrature formula
    // quadrature formula id:
    QuadFormula1D QFId;
    int N_LinePoints; // number of qudrature points
    double *LineWeights, *zeta; // quadrature weights and points (on [-1,1])
    if(assemble_on_return)
    {
      // do a quadrature
      
      // polynomial degree of finite element, needed for choosing an 
      // appropriate quadrature formula, 2 is the polynomial degree of 
      // functions on the interface
      const int fe_degree = 2 * pBf->GetPolynomialDegree();
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
    pBf->MakeRefElementData(QFId);
    
    // compute length of the edge
    const double hE = thisEdge->GetLength();
    
    // index of this edge in the adjacent cell in Darcy subdomain
    const int eI = thisEdge->GetIndexInNeighbor(d_cell); // edge Index
    // normal (pointing out of the Darcy subdomain) and tangential 
    double nx, ny;
    getNormal(d_cell, thisEdge, nx, ny);
    
    // values of all functions and derivatives on reference element
    double **pref, **pxiref, **petaref;
    if(assemble_on_return)
    {
      pref = TFEDatabase2D::GetJointValues2D(pBf->GetID(), QFId, eI);
      pxiref = TFEDatabase2D::GetJointDerivatives2D(pBf->GetID(), QFId, eI,
                                                    D10);
      petaref = TFEDatabase2D::GetJointDerivatives2D(pBf->GetID(), QFId, eI,
                                                     D01);
    }
    else
    {
      // we could define a quadrature rule with only the three points, but
      // the functions 'TFEDatabase2D::GetJointDerivatives2D' only know 
      // predefined quadrature formulas, we have to do this on our own here
      // disadvantage: these results are not stored in the database and have to
      // be recomputed on every interface edge
      pref = new double*[N_LinePoints]; // == 3
      pxiref = new double*[N_LinePoints]; // == 3
      petaref = new double*[N_LinePoints]; // == 3
      for(int i = 0; i < N_LinePoints; i++)
      {
        pref[i] = new double[N_pBf];
        pxiref[i] = new double[N_pBf];
        petaref[i] = new double[N_pBf];
      }
      pBf->GetDerivatives(D00, N_LinePoints, zeta, eI, pref);
      pBf->GetDerivatives(D10, N_LinePoints, zeta, eI, pxiref);
      pBf->GetDerivatives(D01, N_LinePoints, zeta, eI, petaref);
    }
    // values and derivatives of all functions at one quadrature point
    //double *uorig = (mixed) ? new double[N_uBf * BaseVectDim] : NULL;
    double *porig = new double[N_pBf];
    double *pxorig = new double[N_pBf];
    double *pyorig = new double[N_pBf];
    
    for(int k = 0; k < N_LinePoints; k++)
    {
      // transform values and derivatives for this quadrature point
      TFEDatabase2D::GetOrigValues(refTrans, zeta[k], pBf, eI, pref[k],
                                   pxiref[k], petaref[k], porig, pxorig,
                                   pyorig);
      
      // quadrature weight and determinant of tranformation
      const double qw = LineWeights[k] * (assemble_on_return ? (hE / 2) : 1.0);
      
      for(int row = 0; row < N_pBf; row++)
      {
        int testDOF = global_DOF_interface[pDOF[row]][0].second;
        if(testDOF == -1)
          continue; // not an interface dof
        if(global_DOF_interface[pDOF[row]][0].first != d_cell)
          // this dof belongs to two interface dofs (this happens for
          // discontinuous interface functions), take the correct one.
          testDOF = global_DOF_interface[pDOF[row]][1].second;
        const double q = porig[row]; // test function
        
        for(int col = 0; col < N_pBf; col++)
        {
          const int ansatzDOF = pDOF[col];
          const double p = porig[col]; // ansatz function
          
          m.m[0][testDOF][ansatzDOF] -= a * q * p * qw;
          m.m[0][testDOF][ansatzDOF] += a * K
                                        * TDatabase::ParamDB->StoDa_gamma_f * q
                                        * (pxorig[col] * nx + pyorig[col] * ny)
                                        * qw;
        }// end for loop over ansatz pressure basis functions
      } // end for loop over test pressure basis functions
    } // end for loop over all quadrature points
    //delete[] uorig;
    delete[] porig;
    delete[] pxorig;
    delete[] pyorig;
    if(!assemble_on_return)
    {
      // delete arrays which have been created on the heap using new ...[];
      delete [] LineWeights;
      delete [] zeta;
      for(int i = 0; i < N_LinePoints; i++) // i < 3
      {
        delete [] pref[i];
        delete [] pxiref[i];
        delete [] petaref[i];
      }
      delete [] pref;
      delete [] pxiref;
      delete [] petaref;
    } // end if !assemble_on_return
  }// end for loop over all interface edges
}

/** ************************************************************************ */
void DarcyPrimal::assemble_Neumann_map_to_interface(local_matrices &m,
                                                    double a)
{
  if(assemble_on_return == false)
  {
    // when returning Neumann data we have to assemble integrals in the domain.
    ErrThrow("Trying to return Neumann data without assembling integrals");
  }
  
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  
  const TFESpace2D & pSpace = this->get_space();
  const int n_active = pSpace.GetN_ActiveDegrees();
  
  // the collection of all Darcy cells
  TCollection* d_coll = pSpace.GetCollection();
  // number of Darcy cells
  const unsigned int n_cells = (unsigned int) d_coll->GetN_Cells();
  
  // loop over all cells, skip cells which are not at the interface, i.e., cells
  // which have no dof on the interface. That means cells which only share a 
  // vertex with the interface are considered, even if they have no common edge
  // with the interface. Then assemble the necessary terms
  for(unsigned int icell = 0; icell < n_cells; icell++)
  {
    TBaseCell *d_cell = d_coll->GetCell(icell);
    
    FE2D pFEId = pSpace.GetFE2D(0, d_cell);
    TFE2D *pFE = TFEDatabase2D::GetFE2D(pFEId);
    // number of Darcy pressure basis functions
    const int N_pBf = pFE->GetN_DOF(); 
    int *pDOF = pSpace.GetGlobalDOF(d_cell->GetCellIndex());
    // Basis functions for Stokes velocity
    TBaseFunct2D * pBf = pFE->GetBaseFunct2D();
    
    RefTrans2D refTrans = pFE->GetRefTransID();
    TFEDatabase2D::SetCellForRefTrans(d_cell, refTrans);
    
    const int fe_degree = 2 * pBf->GetPolynomialDegree(); // 2 -> bilinear form
    QuadFormula2D QF_2D_Id = TFEDatabase2D::GetQFFromDegree(
        fe_degree, TFEDatabase2D::GetRefElementFromFE2D(pFEId));
    int N_qPoints; // number of qudrature points
    double *qWeights, *x, *y; // quadrature weights and points on reference cell
    {
      TQuadFormula2D *qf2 = TFEDatabase2D::GetQuadFormula2D(QF_2D_Id);
      qf2->GetFormulaData(N_qPoints, qWeights, x, y);
      // qf2 no longer needed, only local scope
    }
    
    // make sure all functions & derivatives are available for this quadrature
    pBf->MakeRefElementData(QF_2D_Id);
    
    // get values and derivatives on reference cell
    double **pref = TFEDatabase2D::GetRefElementValues(pBf->GetID(), QF_2D_Id,
                                                       D00);
    double **pxiref = TFEDatabase2D::GetRefElementValues(pBf->GetID(), QF_2D_Id,
                                                         D10);
    double **petaref = TFEDatabase2D::GetRefElementValues(pBf->GetID(),
                                                          QF_2D_Id, D01);
    double *porig = new double[N_pBf];
    double *pxorig = new double[N_pBf];
    double *pyorig = new double[N_pBf];
    
    // Points in this cell, determinant of the transformation
    double *X, *Y, *det;
    X = new double[N_qPoints];
    Y = new double[N_qPoints];
    det = new double[N_qPoints];
    TFEDatabase2D::GetOrigFromRef(refTrans, N_qPoints, x, y, X, Y, det);
    
    // assemble
    for(int k = 0; k < N_qPoints; k++) // loop over all quadrature points
    {
      // transformation to this grid cell
      // we use NULL in the 'TGridCell'-argument because it is not needed. 
      TFEDatabase2D::GetOrigValues(refTrans, x[k], y[k], pBf, d_coll, NULL,
                                   pref[k], pxiref[k], petaref[k], porig,
                                   pxorig, pyorig);
      
      double qw = qWeights[k] * det[k];
      
      for(int row = 0; row < N_pBf; row++)
      {
        int testDOF = global_DOF_interface[pDOF[row]][0].second;
        if(testDOF == -1)
          continue; // not an interface dof
        if(global_DOF_interface[pDOF[row]][0].first != d_cell)
          // this dof belongs to two interface dofs (this happens for
          // discontinuous interface functions), take the correct one.
          testDOF = global_DOF_interface[pDOF[row]][1].second;
        
        // local Darcy pressure dof which is exactly this interface dof (on the 
        // interface)
        for(int col = 0; col < N_pBf; col++)
        {
          int ansatzDOF = pDOF[col];
          // minus here due to different normal in Darcy subdomain
          m.m[0][testDOF][ansatzDOF] -= a * K
                                        * (pxorig[row] * pxorig[col]
                                           + pyorig[row] * pyorig[col])
                                        * qw;
        }
        
        // take care of 'map_sol2eta_rhs'
        if(pDOF[row] > n_active)
        {
          // we need to assemble (f,q) here, which is not present in the 
          // right hand side (it is only present for active dofs)
          
          // the call to the linear Coefficient function should be done before
          // the loop over the quadrature points. However we only want to call
          // that function when it is necessary. So here we do the entire
          // quadrature (ugly)
          if(k == 0)
          {
            // note that, if you take the solution of this Dirichlet Darcy 
            // problem and map it to eta, using 'map_solution_to_interface', 
            // then take that eta as Neumann boundary data of a Darcy problem,
            // the non-active dofs are negelected anyway. That means you can 
            // do anything in this entry in eta. This quadrature is therefore
            // not necessary in that case. However it might be important, if 
            // eta is passed to a Stokes problem as interface data (where the 
            // corresponding dof can be active).
            double **coeffs = new double*[N_qPoints];
            for(int i = 0; i < N_qPoints; i++)
              coeffs[i] = new double[5];
            
            CoeffFct2D* linCoeffs = this->get_example().get_coeffs();
            linCoeffs(N_qPoints, X, Y, NULL, coeffs);
            double fq = 0.0; // assemble (f,q) into this double
            for(int i = 0; i < N_qPoints; i++)
            {
              // note that pref and porig are the same, otherwise we could not
              // access porig at other quadrature points than the 'k'-th
              fq += coeffs[i][4] * pref[i][row] * qWeights[i] * det[i];
            }
            // minus here due to different normal in Darcy subdomain
            map_sol2eta_rhs(testDOF) -= a * fq;
            
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
            // minus here due to different normal in Darcy subdomain
            map_sol2eta_rhs(testDOF) = -a * this->getRhs()[pDOF[row]];
          }
        }
        
      }
    } // loop over all quadrature points
    delete[] X;
    delete[] Y;
    delete[] det;
    delete[] porig;
    delete[] pxorig;
    delete[] pyorig;
  } // loop over all cells
  if(typeOf_bci == Dirichlet)
  {
    // remove entries in rhs, so that only eta will be there
    for(int i = 0; i < n_active; i++)
    {
      int testDOF = global_DOF_interface[i][0].second;
      if(testDOF == -1)
        continue; // not an interface dof
      this->getRhs()[i] = 0;
      IrhsDARCY[i] = 0;
    }
  }
}

/** ************************************************************************ */
void DarcyPrimal::map_solution_to_interface(InterfaceFunction &eta, double a,
                                             bool old)
{
  if(map_sol2eta == NULL)
  {
    ErrMsg("map_sol2eta has not been assembled yet in a " << typeOf_bci << 
           " Darcy problem. It is done now.");
    prepare_map_solution_to_interface(eta);
  }
  
  const double* const solution = old ? this->getSolOld() : this->getSol();
  //const double* const solution = solution_direct->get_entries();

  /*this->get_sol()->print("Darcy_sol");
  eta.PrintVals("eta_before");*/
  
  // multiply solution with linear part of the map
  map_sol2eta->multiply(solution, eta.get_entries(), a);
  
  /*eta.PrintVals("eta_after");
  InterfaceFunction eta2(eta);
  eta2.add(this->getP(), 0, -1);
  eta2.add(this->getP(), 1, TDatabase::ParamDB->StoDa_gamma_f
           * TDatabase::ParamDB->SIGMA_PERM);
  eta -= eta;
  eta2.PrintVals("eta2");*/
  
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
void DarcyPrimal::update(InterfaceFunction& eta_f, InterfaceFunction* eta_p)
{
  // update eta_f
  const int solution_strategy = TDatabase::ParamDB->StoDa_solutionStrategy;
  
  if(solution_strategy == 1 || solution_strategy == -1)
  {
    const int updatingProcedure = TDatabase::ParamDB->StoDa_updatingStrategy;
    const double theta = TDatabase::ParamDB->StoDa_theta_f; // damping
    
    switch(updatingProcedure)
    {
      case 1: // damping with old interface function
      {
        eta_f.scale(1 - theta);
        this->map_solution_to_interface(eta_f, theta);
        break;
      }
      case 2: // damping with old Darcy solution
      {
        eta_f.reset();
        this->map_solution_to_interface(eta_f, 1.0 - theta, true);
        this->map_solution_to_interface(eta_f, theta);
        break;
      }
      case 3: // C-RR
      {
        if(eta_p == NULL)
          ErrThrow("For the C-RR method, you have to specify eta_p");
        
        double gamma_quotient = TDatabase::ParamDB->StoDa_gamma_f
                                / TDatabase::ParamDB->StoDa_gamma_p;
        eta_f.scale(1 - theta);
        // get Dirichlet data
        this->map_solution_to_interface(eta_f, -theta * (1 + gamma_quotient));
        eta_f.add(eta_p, theta * gamma_quotient);
        break;
      }
      case 4: // D-RR
      {
        eta_f.scale(1 - theta);
        // get Dirichlet data and normal derivatives 
        this->map_solution_to_interface(eta_f, theta);
        break;
      }
      default:
        ErrThrow("unknown updating strategy ", updatingProcedure);
        break;
    }
  }
  else if(solution_strategy == 2 || solution_strategy == -2)
  {
    // solving fixed point equation on interface 
    // turn off damping here, it is done in 'solve_fixed_point(eta)'
    eta_f.reset();
    this->map_solution_to_interface(eta_f, 1.0);
  }
  else if(solution_strategy == 3 || solution_strategy == -3)
  {
    // solving a Stecklov-Poincare equation
    // this is done without this update function
    ErrThrow("Stecklov-Poincare");
    this->map_solution_to_interface(eta_f, 1.0);
  }
}


/** ************************************************************************ */
void DarcyPrimal::create_etaToBd(const InterfaceFunction &eta)
{
  const TFESpace2D & pSpace = this->get_space();
  
  const int N_iDOF = eta.length();
  
  // everything needed to call the constructor of TStructure
  const int n_rows = IrhsDARCY.length();
  const int n_cols = N_iDOF;
  int N_entries = 0;
  int *cols, *rows = new int[n_rows + 1];
  
  // for each row (i.e, Darcy dof) we have one vector containing all 
  // (interface) dofs which couple with this dof. For all dofs which do not 
  // belong to a cell adjacent to the interface, the corresponding vector will 
  // be of length 0
  std::vector<std::vector<int> > coupledDOF(n_rows);
  
  // array containing the number of basis functions for all finite elements
  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  
  // Id of the Stokes space
  const int D_id = pSpace.GetCollection()->GetCell(0)->GetReference_ID();
  
  // loop over the interface edges
  for(unsigned int iEdge = 0; iEdge < interface.size(); iEdge++)
  {
    const TInnerInterfaceJoint * thisEdge = interface[iEdge];
    TBaseCell *d_cell = thisEdge->GetNeighbour(0);
    if(d_cell->GetReference_ID() != D_id)
      d_cell = thisEdge->GetNeighbour(1);
    
    // pressure
    FE2D pFEId = pSpace.GetFE2D(0, d_cell);
    const int N_pBaseFunct = N_BaseFunct[pFEId];
    int *pDOF = pSpace.GetGlobalDOF(d_cell->GetCellIndex());
    // interface
    const int N_iBaseFunct = 3; //abs(eta.spaceType)+1;
    
    // loop over all degrees of freedom of this interface function on this edge
    for(int j = 0; j < N_iBaseFunct; j++)
    {
      const int iDOF = eta.getDOF(iEdge, j);
      // loop over all pressure degrees of freedom of this Darcy cell
      for(int jS = 0; jS < N_pBaseFunct; jS++)
      {
        const int ansatz_pDOF = pDOF[jS];
        // pressure component
        coupledDOF[ansatz_pDOF].push_back(iDOF);
      }
    }
  }
  
  // fill the array 'rows' and compute 'N_entries'
  rows[0] = 0;
  for(int i = 0; i < n_rows; i++)
  {
    // sort the interface dofs for the i-th Darcy dof
    std::sort(coupledDOF[i].begin(), coupledDOF[i].end());
    // remove all duplicates 
    const std::vector<int>::iterator it = unique(coupledDOF[i].begin(),
                                                 coupledDOF[i].end());
    // resize to the size without duplicates
    coupledDOF[i].resize(std::distance(coupledDOF[i].begin(), it));
    // number of interface dofs which couple with this Darcy dof  
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
  std::shared_ptr<TStructure> Cstructure(new TStructure(n_rows, n_cols, 
                                                        N_entries, cols, rows));
  etaToBd.reset(new TMatrix(Cstructure)); // empty matrix
  
  Assemble_etaToBd(eta);
  if(TDatabase::ParamDB->StoDa_periodicBoundary)
    makePeriodicBoundary(etaToBd);
}

/** ************************************************************************ */
void DarcyPrimal::Assemble_etaToBd(const InterfaceFunction &eta)
{
  const TFESpace2D& pSpace =  this->get_space();
  
  if(global_DOF_interface.empty())
    find_global_DOF_interface(eta);
  
  if(!assemble_on_input)
  {
    // number of Darcy pressure degrees of freedom
    //int N_presDOF = pSpace.GetN_DegreesOfFreedom();
    const int N_active = pSpace.GetActiveBound();
    
    const int n_dofs = pSpace.GetN_DegreesOfFreedom();
    double val = 1.0;
    if(typeOf_bci == Robin) // Robin problems
      val /= TDatabase::ParamDB->StoDa_gamma_p;
    
    for(int row = 0; row < n_dofs; row++)
    {
      int interface_dof = global_DOF_interface[row][0].second;
      if(interface_dof != -1 && row < N_active)
      {
        //(*etaToBd)(idof, interface_dof) = val;
        for(unsigned int i = 0; i < global_DOF_interface[row].size(); i++)
        {
          (*etaToBd)(row, global_DOF_interface[row][i].second) = val;
        }
      }
    }
  }
  else
  {
    if(typeOf_bci == Dirichlet && TDatabase::ParamDB->StoDa_weakGamma <= 0.0)
    {
      // do essentially the same as in the 'assemble_on_input==false' case
      const int N_active = pSpace.GetActiveBound();
      // set Dirichlet dofs exactly, instead of using weak Dirichlet conditions)
      const int n_dof = pSpace.GetN_DegreesOfFreedom();
      for(int row = 0; row < n_dof; row++)
      {
        int interface_dof = global_DOF_interface[row][0].second;
        if(interface_dof == -1 || row >= N_active)
          continue; // not an interface dof, or non-active dof
        (*etaToBd)(row, interface_dof) = 1.0; // set right hand side to eta
        // set rhs to zero here, so that only eta will be there
        this->getRhs()[row] = 0.0; 
      }
      return; // no more assembling
    }
    
    // assemble (eta,psi)_Gamma here
    // number of Darcy pressure degrees of freedom
    //int N_presDOF = pSpace.GetN_DegreesOfFreedom();
    const int N_active = pSpace.GetActiveBound();
    
    // Id of the Darcy space
    const int D_id = pSpace.GetCollection()->GetCell(0)->GetReference_ID();
    
    // everything needed for local assembling
    local_edge_assembling l;
    
    for(unsigned int iEdge = 0; iEdge < interface.size(); iEdge++)
    {
      const TInnerInterfaceJoint * thisEdge = interface[iEdge];
      TBaseCell *d_cell = thisEdge->GetNeighbour(0);
      if(d_cell->GetReference_ID() != D_id)
        d_cell = thisEdge->GetNeighbour(1);
      
      FE2D pFEId = pSpace.GetFE2D(0, d_cell);
      TFE2D *pFE = TFEDatabase2D::GetFE2D(pFEId);
      // number of Darcy pressure basis functions
      const int N_pBf = pFE->GetN_DOF(); 
      int *pDOF = pSpace.GetGlobalDOF(d_cell->GetCellIndex());
      // Basis functions for Stokes velocity
      TBaseFunct2D * pBf = pFE->GetBaseFunct2D();
      
      RefTrans2D refTrans = pFE->GetRefTransID();
      TFEDatabase2D::SetCellForRefTrans(d_cell, refTrans);
      
      // polynomial degree of finite element, needed for choosing an 
      // appropriate quadrature formula
      const int fe_degree = abs(eta.getSpaceType())
                            * pBf->GetPolynomialDegree();
      // get the quadrature formula
      // quadrature formula id:
      QuadFormula1D QFId = TFEDatabase2D::GetQFLineFromDegree(fe_degree);
      int N_LinePoints; // number of qudrature points
      double *LineWeights, *zeta; // quadrature weights and points (on [-1,1])
      {
        TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(QFId);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
        // qf1 no longer needed, only local scope
      }
      
      // make sure all functions & derivatives are available for this quadrature
      pBf->MakeRefElementData(QFId);
      
      // compute length of the edge
      const double hE = thisEdge->GetLength();
      
      // index of this edge in the adjacent cell in Darcy subdomain
      const int eI = thisEdge->GetIndexInNeighbor(d_cell); // edge Index
      // normal (pointing out of the Darcy subdomain) and tangential 
      double nx, ny, tx, ty;
      getNormal(d_cell, thisEdge, nx, ny);
      thisEdge->GetTangent(tx, ty);
      
      double **pref = TFEDatabase2D::GetJointValues2D(pBf->GetID(), QFId, eI);
      double **pxiref = TFEDatabase2D::GetJointDerivatives2D(pBf->GetID(), QFId,
                                                             eI, D10);
      double **petaref = TFEDatabase2D::GetJointDerivatives2D(pBf->GetID(),
                                                              QFId, eI, D01);
      
      // values and derivatives of all functions at one quadrature point
      double *porig = new double[N_pBf];
      double *pxorig = new double[N_pBf];
      double *pyorig = new double[N_pBf];
      
      l.nt.set(nx, ny, tx, ty);
      l.hE = hE;
      
      for(int k = 0; k < N_LinePoints; k++)
      {
        TFEDatabase2D::GetOrigValues(refTrans, zeta[k], pBf, eI, pref[k],
                                     pxiref[k], petaref[k], porig, pxorig,
                                     pyorig);
        
        // quadrature weight and determinant of tranformation
        l.qw = LineWeights[k] * (hE / 2);
        
        for(int row = 0; row < N_pBf; row++)
        {
          l.testDOF = pDOF[row];
          if(l.testDOF >= N_active)
            continue; // do not write anything in non-active rows
          l.at.setTest(porig[row], pxorig[row], pyorig[row]);
          
          for(int col = 0; col < N_pBf; col++)
          {
            l.ansatzDOF = global_DOF_interface[pDOF[col]][0].second;
            if(l.ansatzDOF == -1)
              continue; // not an interface dof
            if(global_DOF_interface[pDOF[col]][0].first != d_cell)
              // this dof belongs to two interface dofs (this happens for
              // discontinuous interface functions), take the correct one.
              l.ansatzDOF = global_DOF_interface[pDOF[col]][1].second;
            
            // derivatives not needed
            l.at.setAnsatz(porig[col], pxorig[col], pyorig[col]);
            
            localAssembleEtaToBd_pressure(l);
          }
        }
      }
      delete[] porig;
      delete[] pxorig;
      delete[] pyorig;
    }
    etaToBd->add(l.m.m[0]);
  }
}

/** ************************************************************************ */
void DarcyPrimal::addEta(const InterfaceFunction& eta, double factor)
{
  if(etaToBd == NULL)
  {
    // etaToBd should have been done already
    ErrMsg("etaToBd has not been assembled in a " << typeOf_bci << 
           " Darcy problem. It will be done now.\n");
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
    this->getRhs()[i] += factor * value;
    //Output::print<1>("add ", i, " ", value*factor);
  }
}

/** ************************************************************************ */
void DarcyPrimal::solve(const InterfaceFunction &eta)
{
  Output::print<1>("Solving a ", typeOf_bci, " Darcy problem");
  //print_assemble_input_output();
  
  // copy rhs without interface integrals
  this->get_rhs() = IrhsDARCY;
  
  // copy old solution for damping later
  solution_old = this->get_solution();
                          
  //darcy_primal->get_rhs().print("Darcy_rhs_before");
  addEta(eta);
  //darcy_primal->get_rhs().print("Darcy_rhs_after");
  
  /* ==========================   Solve Darcy   ========================== */
  this->CD2D::solve();
  
  if(TDatabase::ParamDB->StoDa_periodicBoundary)
    checkPeriodicDOFs();
}

/** ************************************************************************ */
void DarcyPrimal::findPeriodicDOFs()
{
  // threshold for two doubles to be equal
  const double eps = 1e-5;
  
  // left and right boundary as periodic boundary
  const double x_left = 0.0;
  const double x_right = 2.0;
  if(TDatabase::ParamDB->EXAMPLE != 2)
  {
    ErrThrow("only the riverbed example is supported for periodic boundaries");
  }
  
  const TFESpace2D * fespace = getP().GetFESpace2D();
  const TCollection* coll = fespace->GetCollection();
  
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
      
      const FE2D FEid1 = fespace->GetFE2D(i1_cell, cell1);
      const TFE2D * FE1 = TFEDatabase2D::GetFE2D(FEid1);
      
      // global indices of all degrees of freedom in this cell 
      const int *globalDOF1 = fespace->GetGlobalDOF(i1_cell);
      // local degrees of freedom which correspond to this edge
      const int* localDOF1 = FE1->GetFEDesc2D()->GetJointDOF(j1);
      const int N_localDOF1 = FE1->GetFEDesc2D()->GetN_JointDOF();
      
      // check if vertex x = x_left
      bool left = false;
      if(fabs(x11 - x_left) > eps || fabs(x12 - x_left) > eps)
      { // one vertex does not lie on the left boundary
        if(fabs(x11 - x_right) > eps || fabs(x12 - x_right) > eps)
          continue; // one vertex does not lie on the right boundary
      }
      else
        left = true;
      
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
            ErrThrow("Error in making periodic boundary. Two different finite ",
                     "elements");
          }
          if(N_localDOF1 != N_localDOF2)
          {
            ErrThrow("Error in making periodic boundary. Different numbers of ",
                     "of dofs on the periodic boundary");
          }
          
          Output::print<3>(" creating a vertical periodic boundary at y=(", y21,
                           ",", y22, ")");
          for(int edge_dof = 0; edge_dof < N_localDOF1; edge_dof++)
          {
            // due to counterclockwise numbering in each cell we have to go 
            // through one edge the opposite way:
            const int dof1 = globalDOF1[localDOF1[N_localDOF1 - 1 - edge_dof]];
            const int dof2 = globalDOF2[localDOF2[edge_dof]];
            
            if(left)
            {
              periodic_dofs[dof1] = dof2;
              //Output::print<1>(" dofs ", dof1, "\t", dof2);
            }
            else
            {
              periodic_dofs[dof2] = dof1;
              //Output::print<1>(" dofs ", dof2, "\t", dof1);
            }
          }
        }
      }
    }
  }
  Output::print<1>("There are ", periodic_dofs.size(), 
                   " periodic Darcy degrees of freedom\n");
}

/** ************************************************************************ */
void DarcyPrimal::makePeriodicBoundary(std::shared_ptr<TMatrix> mat, 
                                        bool darcyMat)
{
  if(periodic_dofs.empty())
  {
    ErrThrow("called DarcyPrimal::makePeriodicBoundary with map ",
             "'periodic_dofs' not yet set");
  }
  
  if(mat == NULL)
  {
    mat = getMat().block(0);
    darcyMat = true;
  }
  const int * const rowPtr = mat->GetRowPtr();
  const int * const colPtr = mat->GetKCol();
  const double * const entries = mat->GetEntries();
  // this map will be passed to "mat->changeRows(new_rows,true)" 
  std::map<int, std::map<int, double> > new_rows;
  
  const int n_rows = mat->GetN_Rows();
  for(int row = 0; row < n_rows; row++)
  {
    const int periodicDof = getPeriodicDOF(row);
    if(periodicDof == -1)
      continue; // not a periodic dof
    // the row with number "row" is added to the row number "periodicDof". 
    // Then row number "row" is replaced by a zero row.
    
    // add the two rows. 
    // loop over all entries in row "row"
    for(int i = rowPtr[row]; i < rowPtr[1 + row]; i++)
    {
      if(entries[i] != 0.0)
        (new_rows[periodicDof])[colPtr[i]] = entries[i];
    }
    // loop over all entries in row "periodicDof"
    for(int i = rowPtr[periodicDof]; i < rowPtr[1 + periodicDof]; i++)
    {
      if(entries[i] != 0.0)
        (new_rows[periodicDof])[colPtr[i]] += entries[i];
    }
    if(darcyMat)
    {
      (new_rows[row])[row] = 1.0;
      (new_rows[row])[periodicDof] = -1.0;
      // rhs set to the prescribed pressure drop
      getRhs()[row] = TDatabase::ParamDB->StoDa_periodicBoundaryPressureDrop;
    }
    else
      new_rows[row]; // empty row
      
  }
  mat->changeRows(new_rows);
}

/** ************************************************************************ */
void DarcyPrimal::checkPeriodicDOFs()
{
  std::map<int, int>::iterator it;
  double const * const vals = this->getP().GetValues();
  for(it = periodic_dofs.begin(); it != periodic_dofs.end(); ++it)
  {
    double diff = vals[it->first] - vals[it->second];
    diff -= TDatabase::ParamDB->StoDa_periodicBoundaryPressureDrop;
    if(fabs(diff) > 1e-6)
      Error("WRONG Darcy here " << it->first << "\t" << it->second << "\t" << 
            setprecision(10) << diff << endl);
  }
}
