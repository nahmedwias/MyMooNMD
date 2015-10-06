#include <InterfaceFunction.h>
#include <FEDatabase2D.h>

#include <sys/stat.h>
#include <sys/types.h>

/** ************************************************************************ */
int main(int argc, char* argv[])
{
  {
  double time = GetTime();
  TDomain Domain;
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  
  /* ==========================================================================
   * read parameter file
   * ========================================================================*/
  Domain.ReadParam(argv[1]);
  OpenFiles();
  
  // ==========================================================================
  // take example from database (input file), see TDatabase::ParamDB->EXAMPLE
  // during the constructor of the example, some values in the database might be
  // changed
  Example_StokesDarcy2D example;
  // write all database parameters to outfile
  Database.WriteParamDB(argv[0]);
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);
  
  // =======================  Parameters ======================================
  check_all_parameters();
  // Neumann_Neumann, Robin_Robin, Dirichlet_Dirichlet, 
  // Neumann_Dirichlet, Dirichlet_Neumann
  const ProblemType pt = (ProblemType) TDatabase::ParamDB->StoDa_problemType;
  const int solution_strategy = TDatabase::ParamDB->StoDa_solutionStrategy;
  const int updating_strategy = TDatabase::ParamDB->StoDa_updatingStrategy;
  // use Raviart-Thomas or Brezzi-Douglas-Marini elements for Darcy part
  bool mixed = false;
  if(pt == Neumann_Dirichlet || pt == Dirichlet_Neumann) 
    mixed = true; // RT for Robin-Robin not yet possible
  // type of interface conditions for tangential Stokes velocity: 
  //  0 : u.t + alpha t.T.n = 0  // Beavers--Joseph--Saffman
  //  1 : u.t = 0                // zero tangential Stokes velocity
  const int icond = TDatabase::ParamDB->StoDa_interfaceType;
  const bool D_RR = (pt == Robin_Robin && updating_strategy == 4);
  const bool C_RR = (pt == Neumann_Neumann && updating_strategy == 3
                    && (solution_strategy == 1 || solution_strategy == -1)); 
  
  
  /* ==========================================================================
   * read boundary parameterization, initialize coarse grid, Refine,
   * ========================================================================*/
  OutPut("Initialize domain\n");
  Domain.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE);
  // vector containing all interface edges (joints)
  std::vector<const TInnerInterfaceJoint *> interface;
  GetInnerInterfaceJoints(interface, Domain); // fill vector 'interface'
  // refine up to user defined coarsest level
  for(int i = 0; i < TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR; i++)
  {
    Domain.RegRefineAll();
    GetInnerInterfaceJoints(interface, Domain); // update interface vector
  }
  /* ==========================================================================
   * get the collections of the two subdomains
   * ========================================================================*/
  OutPut("\nSome info on the subdomains\n");
  {
    // divide the collection
    TCollection *coll_Darcy = Domain.GetCollection(It_Finest, 0, 1);
    TCollection *coll_NSE = Domain.GetCollection(It_Finest, 0, 2);
    
    // write out some information on the grid.
    OutPut(" number of all cells       " << coll_NSE->GetN_Cells() 
                                          + coll_Darcy->GetN_Cells() << endl);
    OutPut(" number of Stokes cells    " << coll_NSE->GetN_Cells() << endl);
    OutPut(" number of Darcy cells     " << coll_Darcy->GetN_Cells() << endl);
    OutPut(" number of interface edges " << interface.size() << endl);
    // hmin is really the minimum of the largest edge of each cell
    // hmin is __not__ the smallest edge, in a uniform grid it equals hmax  
    {
      double hmin, hmax;
      coll_Darcy->GetHminHmax(&hmin, &hmax);
      OutPut(" Darcy  h_min " << hmin << ", hmax " << hmax << endl);
      coll_NSE->GetHminHmax(&hmin, &hmax);
      OutPut(" Stokes h_min " << hmin << ", hmax " << hmax << endl);
    }
    delete coll_Darcy;
    delete coll_NSE;
  }
  if(TDatabase::ParamDB->WRITE_PS)
  {
    // get the total collection on the finest level (both subdomains)
    TCollection *coll_tot = Domain.GetCollection(It_Finest, 0);
    // if you want cell numbers drawn in PS file: uncomment the line 104
    // "#define __GRIDCELL_WITH_NUMBERS__" in src/Geometry/GridCell.C
    std::ostringstream os(std::ios::trunc);
    os << TDatabase::ParamDB->BASENAME << "grid.ps" << ends;
    Domain.PS(os.str().c_str(), coll_tot);
    delete coll_tot;
  }
  OutPut("\n");
  /* =======================   inteface function   ========================= */
  // if D-RR, then -2 (discontinuous) otherwise 2 (continuous)
  int spaceType = D_RR ? -2 : 2;
  
  // no discontinuous interface functions needed, because even for D-RR, the 
  // assembling of the interface terms is done during output of stokes and darcy
  // classes.
  //int spaceType = -2;
  // 'eta' is needed for solution_strategy 2 and 3
  InterfaceFunction eta(interface, spaceType);
  OutPut("number of dofs on the interface: " << eta.length() << endl);
  Out("We use " << (spaceType > 0 ? "" : "dis") << 
      "continuous interface functions\n", 1);
  OutPut(endl);
  
  /* ==========================================================================
   * build matrices , assemble all terms (except for the interface integrals)
   * the matrix does not change during the interface iteration !!
   * ========================       NSE        ============================= */
  std::vector<InterfaceCondition> interface_conditions;
  if(!D_RR && !pt == Robin_Robin)
    interface_conditions.push_back(Neumann);
  if(solution_strategy == 3 || solution_strategy == -3) // Stecklov-Poincare
    interface_conditions.push_back(Dirichlet);
  if(solution_strategy != 3 && solution_strategy != -3 // Stecklov-Poincare
       && (D_RR || pt == Robin_Robin || C_RR))
    interface_conditions.push_back(Robin);
  //interface_conditions.push_back(weakRobin);
  //interface_conditions.push_back(DirichletSTAB);
  
  std::map<InterfaceCondition, StokesProblem*> ns_problems;
  std::shared_ptr<Example_NSE2D> ns_example = example.get_stokes_example();
  for(unsigned int i = 0; i < interface_conditions.size(); i++)
  {
    OutPut("intitializing a " << interface_conditions[i] << 
           " (Navier-) Stokes problem\n");
    StokesProblem * ns = new StokesProblem(Domain, ns_example, 
                                           interface_conditions[i],
                                           interface, icond);
    if(TDatabase::ParamDB->SC_VERBOSE > 0) ns->print_assemble_input_output();
    ns->assemble();
    ns->find_global_DOF_interface(eta);
    ns->mat_NSE_iIntegrals(); // assemble all interface integrals into matrix
    if(TDatabase::ParamDB->StoDa_periodicBoundary) // change this if necessary)
    {
      if(!C_RR) ns->get_matrix().block(8)->add(0, 0, 1e30);
      if(TDatabase::ParamDB->LAPLACETYPE == 1)
      {
        // assemble the additional term (nu (grad u)^T n , v)
        // for Poisseuille flow this is needed to get the correct solution
        //map<int, double> comps;
        //comps.insert(pair<int, double>(3, -1.0));
        //comps.insert(pair<int, double>(5, -1.0));
        //ns->assemble_boundary_term("n nu graduT n v", comps);
        ErrMsg("Stokes_2D::assemble_boundary_term not yet ported from MooNMD");
        exit(1);
      }
      ns->findPeriodicDOFs();
      ns->makePeriodicBoundary();
    }
    ns->reorder_rows();
    ns->copy_rhs(); // copy non-interface-part of rhs
    ns->prepare_input_output_maps(eta);
    if(interface_conditions[i] == Dirichlet
       && (TDatabase::ParamDB->EXAMPLE == 0
           || TDatabase::ParamDB->EXAMPLE == 1))
      ns->get_matrix().block(8)->add(0, 0, 1e30);
    ns_problems[interface_conditions[i]] = ns;
  }
  OutPut(endl);
  
  /* ==========================     Darcy     ============================== */
  // needed for primal Darcy
  TDatabase::ParamDB->ANSATZ_ORDER = TDatabase::ParamDB->VELOCITY_SPACE;
  if(mixed)
    TDatabase::ParamDB->VELOCITY_SPACE = 1002;
  std::shared_ptr<Example_CD2D> d_example = example.get_primal_darcy_example();
  std::map<InterfaceCondition, DarcyPrimal*> d_problems;
  for(unsigned int i = 0; i < interface_conditions.size(); i++)
  {
    OutPut("intitializing a " << (mixed ? "mixed " : "primal ")
           << interface_conditions[i] << " Darcy problem\n");
    DarcyPrimal * d = new DarcyPrimal(Domain, *d_example.get(), 
                                      interface_conditions[i], interface, 1);
    if(TDatabase::ParamDB->SC_VERBOSE > 0) d->print_assemble_input_output();
    d->assemble(); // assemble linear terms
    d->find_global_DOF_interface(eta);
    d->mat_DARCY_iIntegrals(); // assemble all interface integrals into matrix
    if(TDatabase::ParamDB->StoDa_periodicBoundary)
    {
      d->findPeriodicDOFs();
      d->makePeriodicBoundary();
    }
    d->copy_rhs();
    d->prepare_input_output_maps(eta); // assemble maps to and from interface
    d_problems[interface_conditions[i]] = d;
  }
  
  /* ==================   Stokes--Darcy problem class   ==================== */
  // system, we're trying to solve
  StokesDarcy2D sd(ns_problems, d_problems);
  
  Out(std::setfill('*') << std::setw(79) << "*\n" << std::setfill(' '), 0);
  Out("Everything has been set up (withing " << GetTime()-time << " seconds), "
      << "ready to solve\n\n", 0);
  /* ===================   compute direct solution   ======================= */
  time = GetTime();
  if(solution_strategy >= 0)
  {
    sd.solveDirect();
    
    double error = ErrorOnInterface(*(sd.c_stokes()), *(sd.c_darcy()), -1);
    OutPut("  Error on interface " << error << endl)
    WriteVtk_and_measureErrors(*(sd.c_stokes()), *(sd.c_darcy()), NULL, -1);
    Out("Done with direct solution of the coupled system within " << 
        GetTime()-time << " seconds\n", 0);
  }
  
  if(TDatabase::ParamDB->SC_VERBOSE) sd.check_equalities();
  
  if(solution_strategy != 0) // not only direct solve
  {
    Out(std::setfill('*') << std::setw(79) << "*\n" << std::setfill(' '), 0);
    Out("Starting to solve iteratively\n\n", 0);
  }
  
  /* =======================   interface iteration   ======================= */
  time = GetTime();
  if(solution_strategy == 1 || solution_strategy == -1)
  {
    // solve iteratively using blockwise Jacobi or blockwise Gauss--Seidel
    // 'eta_f' and 'eta_p' are needed for solution_strategy 1, they are copies
    // of 'eta' which iteself is not needed here.
    InterfaceFunction eta_f(eta);
    InterfaceFunction eta_p(eta);
    
    int maxit = TDatabase::ParamDB->StoDa_nIterations;
    for(int it = 0; it < maxit; it++)
    {
      OutPut("\nInterface Iteration: " << it << endl);
      OutPut(" Memory: " << setw(8) << GetMemory()/(1048576.0) << " MB\n");
      sd.solve_coupled_system(eta_f, eta_p);
      
      // check stopping criterion, compute errors, make pictures
      if(sd.stopIteration(it))
        break;
    } // end loop for coupling
    Out("Done with classical iterative method within " << GetTime()-time <<
        " seconds\n", 0);
  }
  else if(solution_strategy == 2 || solution_strategy == -2)
  {
    // solve a fixed point formulation
    int maxit = TDatabase::ParamDB->StoDa_nIterations;
    for(int it = 0; it < maxit; it++)
    {
      OutPut("\nInterface Iteration: " << it << endl);
      OutPut(" Memory: " << setw(8) << GetMemory()/(1048576.0) << " MB\n");
      sd.solve_fixed_point(eta);
      
      // check stopping criterion, compute errors, make pictures
      if(sd.stopIteration(it))
        break;
    } // end loop for coupling
    OutPut("norm of final eta " << eta.norm() << endl);
    Out("Done with fixed point iteration within " << GetTime()-time <<
         " seconds\n", 0);
  }
  else if(solution_strategy == 3 || solution_strategy == -3)
  {
    OutPut("\nSolving the Stecklov-Poincare equation\n");
    // solve a Stecklov-Poincare equation using template iterative solvers
    sd.solve_Stecklov_Poincare(eta);
    
    const double error = ErrorOnInterface(*(sd.stokes()), *(sd.darcy()), 0);
    OutPut("  Error on interface " << error << endl)
    WriteVtk_and_measureErrors(*(sd.stokes()), *(sd.darcy()), NULL, 0);
    Out("Done with Steklov-Poincare iteration within " << GetTime()-time <<
        " seconds\n", 0);
  }
  /* ============================   clean up   ============================= */
  for(unsigned int i = 0; i < interface_conditions.size(); i++)
  {
    delete ns_problems.at(interface_conditions[i]);
    delete d_problems.at(interface_conditions[i]);
  }
  }
  OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0) << " MB" << endl);
  OutPut("used time: " << GetTime() << endl);
  OutPut("Program finished normally\n");
  CloseFiles();
  return 0;
}
