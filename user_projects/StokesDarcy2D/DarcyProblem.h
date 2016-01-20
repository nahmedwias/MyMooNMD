#ifndef DARCYPROBLEM_H
#define DARCYPROBLEM_H

class DarcyProblem;

#include <Darcy2D.h>
#include <CD2D.h>
#include <InterfaceFunction.h>
#include <Example_StokesDarcy2D.h>

/** Description of the class DarcyProblem
 * More or less everything related to the Darcy subproblem is done within
 * an instance of this class. Mainly that is
 * - Assembling of interface terms into matrix and right hand side
 * - solve the Darcy subproblem
 */
class DarcyProblem
{
  private:
  
  // inheritance is not useful here because we want to have the primal and the 
  // mixed formulation possible at the same (compile) time. So we include two 
  // member variables for the two formulations
  CD2D* darcy_primal;
  Darcy2D* darcy_mixed;
  std::shared_ptr<Example_CD2D> example_primal;
  std::shared_ptr<Example_Darcy2D> example_mixed;
  
  // The type of boundary condition on the interface (bci) determines which 
  // integrals over the interface should be assembled.
  InterfaceCondition typeOf_bci;
  // use mixed formulation in Darcy subdomain (Raviart-Thomas (RT) or 
  // Brezzi-Douglas-Marini (BDM) elements)
  bool mixed;
  
  // right hand side which only stores the parts of the right hand side which
  // do not originate from interface terms. 
  BlockVector* IrhsDARCY;
  
  // solution (previous iterate)
  BlockVector* solution_old;
  TFEFunction2D* p_Darcy_old;
  TFEFunction2D* u_Darcy_old;
  
  // solution (from big coupled system)
  BlockVector* solution_direct;
  TFEFunction2D* p_Darcy_direct;
  TFEFunction2D* u_Darcy_direct;
  
  // vector of all interface edges
  std::vector<TInnerInterfaceJoint *>& interface;
  
  // coupling matrix, mapping a Stokes solution into the Darcy domain, 
  // e.g. -(u.n,psi)
  // it is assembled only if a direct solution is computed
  std::shared_ptr<TMatrix> stokesToDarcy;
  
  // combined matrix, ready to be passed to solver
  // for RT==true:  ( A BT )   otherwise just A
  //                ( B 0  )
  // for C-RR the iteration matrix is different from the matrix used in the 
  // direct method. 'mat_aux' is the one from the direct method
  mutable std::shared_ptr<TMatrix> mat_aux;
  
  // matrix which maps an interface function eta to boundary conditions
  // In each iteration step we add the product etaToBd*eta to the right hand 
  // side. This corresponds to setting eta as boundary data. (etaToBd depends
  // on the boundary condition, and on the space on the interface)
  std::shared_ptr<TMatrix> etaToBd;
  
  // after solving multiply the solution with this matrix to get the data for 
  // the interface variable. If this is a Neumann problem, then this matrix is
  // a mass matrix. If this is a Dirichlet problem, then this matrix is the 
  // system matrix without the weak Dirichlet conditions and without all rows
  // not corresponsing to the interface. For a Robin problem, this could be a 
  // linear combination of both.
  std::shared_ptr<TMatrix> map_sol2eta;
  BlockVector *map_sol2eta_rhs;
  
  // for the riverbed example with "nonhomogeneous Robin conditions" this map 
  // contains all pairs of velocity dofs which have to be identified.
  std::map<int,int> periodic_dofs;
  
  // The operator represented by this class is affine linear. In order to use 
  // the iterative solvers, we need the linear part of this operator. To 
  // achieve this we save the homogeneous (nonzero) solution. So that the 
  // linear part of this operator can be computed by substracting this 
  // homogeneous solution.
  InterfaceFunction* eta_hom;
  
  
  // some degrees of freedom (dofs) on the interface directly correspond to some 
  // dofs in this Darcy subdomain (that means restricte to the interface they 
  // have the same values). This variable stores this information and is written
  // in 'find_global_DOF_interface(eta)'.
  // global_DOF_interface.size() is the overall number of Darcy dofs. For one
  // such dof 'i' the size of the vector global_DOF_interface[i] is the number 
  // of interface dofs it it eqal to. This is usually one. However if the 
  // interface function eta is discontinuous (as is needed for D-RR), then it is
  // the number of cells this dof is in. In that case the i-th dof is equal to
  // the 'global_DOF_interface[i][0].second'-th dof on the interface restricted 
  // to the cell 'global_DOF_interface[i][0].first'. 
  std::vector<std::vector<std::pair<TBaseCell *, int> > > global_DOF_interface;
  
  
  // When data from one subdomain is passed to the other, usually some integrals
  // have to be evaluated. If Neumann (or Robin) data is returned these are 
  // integrals over (parts of) the domain, if Dirichlet data is returned this is
  // an integral over the interface. While in the first part this has to be done
  // in the originating subdomain, the second case can be assembled in either
  // the originating or the receiving subdomain. 
  // If we solve the coupled problem or the fixed point formulation, we do all 
  // the necessary assemblings in the originating subdomain, such the receiving 
  // subdomain only needs to add some numbers to the rhs-vector. However in the
  // Stecklov-Poincare formulation things get trickier. Here it depends on which
  // of the two formulations is used (depends on StoDa_StokesFirst). The 
  // following to booleans simply save this information, so we don't have to 
  // check the database so often.  
  bool assemble_on_input;
  bool assemble_on_return;
  
  
  // private member methods
  
  /**
   * create 'etaToBd'. 
   */
  void create_etaToBd(const InterfaceFunction &eta);
  
  /** assembling of 'etaToBd'. This is called during the function 
   * create_etaToBd(eta) 
   * */
  void Assemble_etaToBd(const InterfaceFunction &eta);

  void localAssembleEtaToBd_velocity(local_edge_assembling &lea);
  void localAssembleEtaToBd_pressure(local_edge_assembling &lea);
  
  
  /**
   * add etaToBd*eta to the right hand side 'IrhsNSE'
   * This imposes boundary data "eta" on the interface during the iteration
   */
  void addEta(const InterfaceFunction& eta, double factor=1.0);
  
  
  /** 
   * local assembling routines for all integrals over the interface contributing 
   * to the matrix. Depending on the type of problem different integrals are
   * assembled.
   */
  void localAssembleMat_velocity_velocity(local_edge_assembling &lea);
  void localAssembleMat_velocity_pressure(local_edge_assembling &lea);
  void localAssembleMat_pressure_pressure(local_edge_assembling &lea);


  /**
   * prepare everything needed to call 'map_solution_to_interface' after solving
   * 
   * This function calls 'find_global_DOF_interface', 
   * 'create_structure_of_map_solution_to_interface', 
   * 'assemble_Dirichlet_map_to_interface' and 
   * 'assemble_Neumann_map_to_interface'.
   * (the last two only if necessary)
   */
  void prepare_map_solution_to_interface(const InterfaceFunction &eta);
  
  void create_structure_of_map_solution_to_interface(unsigned int);
  void assemble_Dirichlet_map_to_interface(local_matrices & m, bool D_RR,
                                           double a = 1.0);
  void assemble_Dirichlet_map_to_interface_D_RR(local_matrices & m, double a);
  void assemble_Neumann_map_to_interface(local_matrices &m, double a = 1.0);
  public:
  
  // constructor
  DarcyProblem(const TDomain& domain, Example_StokesDarcy2D* ex,
               InterfaceCondition t, std::vector<TInnerInterfaceJoint *>& in,
               bool rt, int reference_id);

  // destructor
  ~DarcyProblem();
  
  
  // assemble linear terms
  void assemble()
  {
    if(mixed)
      darcy_mixed->assemble();
    else
      darcy_primal->assemble();
  }
  
  /**
   * copy the right hand side which has been assembled in the Darcy subdomain
   * to the member variable 'IrhsNSE'. This means saving all contributions to 
   * the right hand side not originating from the interface conditions. This 
   * has to be called before the interface iteration. 
   */
  void copy_rhs();
  
  /**
   * Assemble all integrals over the interface. These only cover the integrals 
   * which contribute to the matrix. Depending on the 'typeOfProblem' 
   * different integrals are assembled. As an option one can specify the 
   * BlockMatrix2D into which the integrals should be assembled. "m==NULL" 
   * means m=&matrixDARCY (default)
   */
  void mat_DARCY_iIntegrals(BlockMatrix * m = NULL);
  
  /** find all Darcy dofs which correspond to interface dofs */
  void find_global_DOF_interface(const InterfaceFunction &eta);
  
  void prepare_input_output_maps(const InterfaceFunction &eta)
  { 
    prepare_map_solution_to_interface(eta);
    create_etaToBd(eta); // do this after 'prepare_map_solution_to_interface'
    /*etaToBd->PrintFull("etaToBd", 2);
    map_sol2eta->PrintFull("map_sol2eta", 2);
    if(map_sol2eta_rhs != NULL)
      map_sol2eta_rhs->print("map_sol2eta_rhs");*/
  }
  
  /** Description of the function solve:
   * Assemble interface integrals and solve for the Darcy part with given 
   * boundary function eta
   */
  void solve(const InterfaceFunction &eta);

  /** Description of the function map_solution_to_interface:
   * After computing a solution one can map it to an interface function eta. 
   * Depending on what type of boundary conditions on the interface are 
   * prescribed, this map means different things. E.g. if this is a Neumann 
   * problem on the interface, this map returns Dirichlet data (and vice versa).
   * 
   * If C is the mapping of the solution to the interface variable eta, then 
   * this function computes eta += a * C(solution)
   * 
   * If 'old' is set to true, the old solution is used instead of the current 
   */
  void map_solution_to_interface(InterfaceFunction &eta, double a = 1.0,
                                 bool old = false);
  
  /** find periodic boundaries dofs. 
   * This fills the map<int,int> 'periodic_dofs' such that a call to 
   * 'getPeriodicDOF(int)' now makes sense
   */
  void findPeriodicDOFs();
  /** enable periodic boundaries in some matrix which has as many rows as the
   * Darcy matrix. This is used for 
   * - D, the Darcy matrix itself (default, i.e., with no arguments)
   * - C, coupling matrix used for the direct solution
   * - E, representing (eta_p,psi), which is added to the rhs during iteration
   * 
   * The method makePeriodicBoundary() (without arguments must be called first)
   * If the second argument is true then there will be ones on the diagonal and 
   * -1 in the off-diagonal entry coupling with this row. 
   */
  void makePeriodicBoundary(std::shared_ptr<TMatrix> mat = nullptr,
                            bool darcyMat = false);
  
  /** check if the given dof is a periodic dof. If no, -1 is returned. If yes, 
   * the dof to which this dof is coupled (periodicDOF) is returned. The row 
   * 'dof' should then be deleted. If mat==getMat().squareBlock(0), the row
   * 'dof' should be replaced by a 1 on the diagonal and -1 on 'periodicDOF'.
   */ 
  int getPeriodicDOF(int dof) const
  {
    std::map<int,int>::const_iterator it = periodic_dofs.find(dof);
    if (it == periodic_dofs.end()) return -1;
    else return it->second;
  }
  
  void checkPeriodicDOFs();
  
  // getters and setters
  std::vector<TInnerInterfaceJoint*>& getInterface() const
  { return interface; }
  
  InterfaceCondition getTypeOf_bci() const
  { return typeOf_bci; }
  
  void setTypeOf_bci(InterfaceCondition type)
  { typeOf_bci = type; }
  
  // true if mixed formulation is used
  bool getRT() const 
  { return mixed; }
  
  const TFEFunction2D& getP() const
  { return ((mixed) ? (darcy_mixed->get_pressure()) : 
                      (darcy_primal->get_function())); }
  
  TFEFunction2D* getPOld() const
  { return p_Darcy_old; }
  
  const TFEFunction2D& getU() const
  { if(mixed) return darcy_mixed->get_velocity(); 
    else ErrThrow("there is no velocity FEFunction"); throw; }
  
  TFEFunction2D* getUOld() const
    { return u_Darcy_old; }
  
  double* getSolOld() const
  { return solution_old->get_entries(); }
  
  BlockVector& get_solution()
  { return (mixed) ? darcy_mixed->get_solution() 
                   : darcy_primal->get_solution(); }
  
  const BlockVector& get_solution() const
  { return (mixed) ? darcy_mixed->get_solution() 
                   : darcy_primal->get_solution(); }
  
  const double* getSol() const
  { return this->get_solution().get_entries(); }

  std::shared_ptr<TMatrix> getCouplingMatrix() const {return stokesToDarcy; }
  void setCouplingMatrix(std::shared_ptr<TMatrix> m) { stokesToDarcy = m; }
  
  std::shared_ptr<TMatrix> get_map_sol2eta() const { return map_sol2eta; }
  std::shared_ptr<TMatrix> get_etaToBd() const { return etaToBd; }
  
  BlockMatrix& getMat() 
  { if(mixed) return darcy_mixed->get_matrix(); 
    else return darcy_primal->get_matrix(); }
  
  std::shared_ptr<TMatrix> getComposedMatForBigSystem()
  {
    if(!mat_aux)
      mat_aux = this->getMat().get_combined_matrix();
    return mat_aux; 
  }
  
  BlockVector & get_rhs()
  {return (mixed) ? darcy_mixed->get_rhs()
                   : darcy_primal->get_rhs(); }
  
  const BlockVector & get_rhs() const 
  {return (mixed) ? darcy_mixed->get_rhs()
                   : darcy_primal->get_rhs(); }
  
  const double* getRhs() const
  { return this->get_rhs().get_entries(); }
  
  double* getRhs()
  { return this->get_rhs().get_entries(); }
  
  void setDirectSol(double *s)
  {
    *solution_direct = s;
  }
  
  BlockVector * get_direct_sol() const
  { return this->solution_direct; }
  
  TFEFunction2D* getPDirect() const
  { return p_Darcy_direct; }
  
  TFEFunction2D* getUDirect() const
  { return u_Darcy_direct; }
  
  const Example2D & get_example() const
  {
    return (mixed) ?  (const Example2D&) darcy_mixed->get_example() 
                   :  (const Example2D&) darcy_primal->get_example();
  }
  
  DoubleFunct2D *const * get_exact() const 
  {
    const std::vector<DoubleFunct2D*> & e = this->get_example().get_exact();
    return &e[0];
  }
  
  InterfaceFunction*& get_homogeneous_solution() {return eta_hom;}
  
  TCollection* get_collection() const 
  { return getP().GetFESpace2D()->GetCollection(); }
  
  void print_assemble_input_output()
  {
    OutPut("Darcy " << typeOf_bci << "\t");
    OutPut((assemble_on_input ? "" : "no ") << "assembling on input\t");
    OutPut((assemble_on_return ? "" : "no ") << "assembling on output\n");
  }
};

/** ************************************************************************ */
/** Description of the function ComputeDarcyOutFlowOnInterface
 * One has to make sure that 
 * int_{\partial\Omega_S} u.n = \int_{\Omega_S} div u.
 * I.e. one has to set 
 * - int_Interface u.n       (mixed formuation in Darcy subdomain)
 * - int_Interface K phi_n   (primal formuation in Darcy subdomain)
 * to a fixed number.
 * 
 * We project K phi_n = u.n such that the above divergence constraint is 
 * fulfilled. This is only necessary for the Dirichlet-Dirichlet problem and 
 * if all other boundary parts in the Stokes subdomain are of Dirichlet-type.
 * 
 * This function computes 
 * - int_Interface u.n       (mixed formuation in Darcy subdomain)
 * - int_Interface K phi_n   (primal formuation in Darcy subdomain) 
 * for the Darcy-velocity. With this one can project  
 */
//double DarcyOutFlowOnInterface(DarcyProblem &d);

#endif //DARCYPROBLEM_H
