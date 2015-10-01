#ifndef STOKESPROBLEM_H
#define STOKESPROBLEM_H

class StokesProblem;

#include <NSE2D.h>
#include <auxiliaryFunctions.h>
#include <InterfaceFunction.h>
#include <map>

/** Description of the class StokesProblem
 * More or less everything related to the Stokes subproblem is done within
 * an instance of this class. Mainly that is
 * - Assembling of interface terms into matrix and right hand side
 * - solve the Stokes subproblem
 */
class StokesProblem : public NSE2D
{
  private:
  // member variables:
  
  // right hand side which stores the Data of the problem, it is assembled only
  // once. Then this DrhsNSE is added to the right hand side for the solver in
  // every iteration step
  BlockVector *DrhsNSE;
  
  // matrix which maps a darcy solution into the Stokes space, e.g. (phi,v.n)
  // it is assembled only if a direct solution is computed
  std::shared_ptr<TMatrix> darcyToStokes; 
  
  // solutions (previous iterate)
  BlockVector* solution_old;
  TFEVectFunct2D *u_NSE_old;
  TFEFunction2D* p_NSE_old;
  
  // solutions (from big coupled system)
  BlockVector* solution_direct;
  TFEVectFunct2D *u_NSE_direct;
  TFEFunction2D* p_NSE_direct;
  
  // The type of boundary condition on the interface (bci) determines which 
  // integrals over the interface should be assembled (normal direction).
  InterfaceCondition typeOf_bci;
  
  // vector containing all interface edges
  std::vector<TInnerInterfaceJoint *>& interface;
  // describing which interface conditions are used for the tangential component
  //  0 : u.t + alpha t.T.n = 0    (Beavers-Joseph-Saffman)
  //  1 : u.t = 0                  (no slip)
  int icond;
  
  // for C-RR the iteration matrix is different from the matrix used in the 
  // direct method. 'mat_aux' stores the one from the direct method
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
  
  // for the riverbed example with "nonhomogeneous periodic boundary conditions" 
  // this map contains all pairs of velocity dofs which have to be identified.
  std::map<int,int> periodic_dofs;
  
  // The operator represented by this class is affine linear. In order to use 
  // the iterative solvers, we need the linear part of this operator. To 
  // achieve this we save the homogeneous (nonzero) solution. So that the 
  // linear part of this operator can be computed by substracting this 
  // homogeneous solution.
  InterfaceFunction* eta_hom;
  
  // some degrees of freedom (dofs) on the interface directly correspond to some 
  // velocity dofs in this Stokes subdomain (that means restricted to the 
  // interface they have the same values). This variable stores this information 
  // and is written in 'find_global_DOF_interface(eta)'.
  // global_DOF_interface.size() is the overall number of Stokes velocity dofs. 
  // For one such dof 'i' the size of the vector global_DOF_interface[i] is the  
  // number of interface dofs it is eqal to. This is usually one. However if the 
  // interface function eta is discontinuous (as is needed for D-RR), then it is
  // the number of cells this dof is in. In that case the i-th dof is equal to
  // the 'global_DOF_interface[i][0].second'-th dof on the interface restricted 
  // to the cell 'global_DOF_interface[i][0].first'. 
  std::vector<std::vector<std::pair<TBaseCell *, int> > > global_DOF_interface;
  // for all dofs on the interface we need a normal vector. It is unclear how to
  // define this for a dof on a vertex on the interface, if the neighboring 
  // normals are different.
  std::vector< std::vector<n_t_vector> > normals;
  
  
  // When data from one subdomain is passed to the other, usually some integrals
  // have to be evaluated. If Neumann (or Robin) data is returned these are 
  // integrals over (parts of) the domain, if Dirichlet data is returned this is
  // an integral over the interface. While in the first case this has to be done
  // in the originating subdomain, the second case can be assembled in either
  // the originating or the receiving subdomain. 
  // If we solve the coupled problem or the fixed point formulation, we do all 
  // the necessary assemblings in the originating subdomain, such that the 
  // receiving subdomain only needs to add some numbers to the rhs-vector. 
  // However in the Stecklov-Poincare formulation things get trickier. Here it 
  // depends on which of the two formulations is used (depends on 
  // StoDa_StokesFirst). The following to booleans simply save this information, 
  // so we don't have to check the database so often.  
  bool assemble_on_input;
  bool assemble_on_return;
  
  
  // member methods:
  
  /**
   * create 'etaToBd'. Here the structure of 'etaToBd' is created and  
   * the function 'Assemble_etaToBd' is called.
   */
  void create_etaToBd(const InterfaceFunction &eta);
  void Assemble_etaToBd(const InterfaceFunction &eta);
  
  /** 
   * local assembling routines for all integrals over the interface contributing 
   * to 'EtaToBd'. Depending on the type of problem different integrals are
   * assembled.
   */
  void localAssembleEtaToBd_velocity(local_edge_assembling &lea);
  void localAssembleEtaToBd_pressure(local_edge_assembling &lea);
  
  /**
   * add etaToBd*eta to the right hand side 'DrhsNSE'
   * This imposes boundary data "eta" on the interface during the iteration
   */
  void addEta(const InterfaceFunction& eta, double factor = 1.0);
  
  
  /* local assembling routines for all integrals over the interface contributing 
   * to the matrix. Depending on the type of problem different integrals are
   * assembled.
   */
  void localAssembleMat_velocity_velocity(local_edge_assembling &lea);
  void localAssembleMat_velocity_pressure(local_edge_assembling &lea);
  void localAssembleMat_pressure_pressure(local_edge_assembling &lea);
  
  /* prepare everything needed to call 'map_solution_to_interface' after solving
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

  /* constructor */
  StokesProblem(const TDomain & domain, Example_NSE2D* ex,
                InterfaceCondition t, std::vector<TInnerInterfaceJoint *>& in,
                int c);
  
  /* destructor */
  ~StokesProblem();
  
  /**
   * Assemble all integrals over the interface. These only cover the integrals 
   * which contribute to the matrix. Depending on the 'typeOfProblem' 
   * different integrals are assembled. As an option one can specify the 
   * BlockMatrix2D into which the integrals should be assembled. "m==NULL" 
   * means m=&matrixNSE (default)
   */
  void mat_NSE_iIntegrals(BlockMatrixNSE2D* m = nullptr);
  
  /** find all Stokes dofs which correspond to interface dofs */
  void find_global_DOF_interface(const InterfaceFunction &eta);
  
  void prepare_input_output_maps(const InterfaceFunction &eta)
  {
    prepare_map_solution_to_interface(eta);
    create_etaToBd(eta);
    /*etaToBd->PrintFull("etaToBd", 2);
    map_sol2eta->PrintFull("map_sol2eta", 2);
    if(map_sol2eta_rhs != NULL)
     map_sol2eta_rhs->print("map_sol2eta_rhs");*/
    /*if(typeOf_bci == Neumann)
    {
      WriteMat((char*)"etaToBd_N.txt", etaToBd);
      WriteMat((char*)"map_sol2eta_N.txt", map_sol2eta);
      if(map_sol2eta_rhs != NULL)
        WriteVec((char*)"map_sol2eta_rhs_N.txt", map_sol2eta_rhs->full(),
                 map_sol2eta_rhs->length());
      WriteMat((char*)"Mat_N.txt", Stokes_2D::getMatrix()->composeBlockMatrix());
    }
    else if(typeOf_bci == Dirichlet)
    {
      WriteMat((char*)"etaToBd_D.txt", etaToBd);
      WriteMat((char*)"map_sol2eta_D.txt", map_sol2eta);
      if(map_sol2eta_rhs != NULL)
        WriteVec((char*)"map_sol2eta_rhs_D.txt", map_sol2eta_rhs->full(),
                 map_sol2eta_rhs->length());
      WriteMat((char*)"Mat_D.txt", Stokes_2D::getMatrix()->composeBlockMatrix());
    }*/
  }
  
  /**
   * copy the right hand side which has been assembled in the Stokes subdomain
   * to the member variable 'DrhsNSE'. This means saving all contributions to 
   * the right hand side not originating from the interface conditions. This 
   * has to be called before the interface iteration. 
   */
  void copy_rhs();
  
  /** 
   * Add boundary data function eta to right hand side and then solve
   */
  void solve(const InterfaceFunction& eta);
  
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
   * Stokes matrix. This is used for 
   * - S, the Stokes matrix itself  (default, i.e., no arguments)
   * - C, coupling matrix used for the direct solution
   * - E, representing (eta_f,u.n), which is added to the rhs during iteration
   * 
   * The method makePeriodicBoundary() (without arguments must be called first)
   * If the second argument is true then there will be ones on the diagonal and 
   * -1 in the off-diagonal entry coupling with this row. If additionally the 
   * third argument is true, then there are zeros put where otherwise 1 and -1
   * are put (as described in the previous sentence). This enlarges the 
   * structure of the matrix 'mat', it is then possible to add two such 
   * matrices. This is currently implemented only due to periodic boundary
   * values (riverbed example).  
   */
  void makePeriodicBoundary(std::shared_ptr<TMatrix> mat = nullptr,
                            bool stokesMat = false, bool p = false);
  
  void checkPeriodicDOFs();
  
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
  
  /* reorder rows, so that for each interface dof, one row corresponds to the
   * normal component and one corresponds to the tangetial component. The 
   * second version of this function is made for the coupling matrix used in 
   * the big system.
   * Otherwise interfaces not aligned with an axis won't work. To do this we
   * will call TMatrix::changeRows() on each matrix with velocity test space
   * 
   * The optional argument 'only_normal' can be set to only alter the row in
   * normal direction leaving the tangential direction unchanged (as is needed 
   * for etaToBd) 
   */
  void reorder_rows(BlockMatrixNSE2D* m = nullptr);
  void reorder_rows(std::shared_ptr<TMatrix> m, bool only_normal = false);
  
  
  /* getters and setters */ 
  
  // return the interface condition (u.t=0  or BJS)
  int getIcond() const 
  { return icond; }
  
  // get normal vector of a specific dof. If there are two you can set i=1
  // at a node which is located at a kink return a vector n which is neither one
  // of the two possible normals, n1, n2, but a (uniqe) linear combination such 
  // that n.n1 = 1 and n.n2 = 1. Note that this 'generalized normal' n does not
  // have length 1.
  n_t_vector& normal(unsigned int dof, unsigned int i = 10);
   
  
  // return true if this dof is located at a kink, otherwise false 
  bool kink(unsigned int dof) const
  { return normals[dof].size() >= 2; }
  
  std::vector<TInnerInterfaceJoint*>& getInterface() const
  { return interface; }
  
  InterfaceCondition getTypeOf_bci() const
  { return typeOf_bci; }
  
  void setTypeOf_bci(InterfaceCondition type)
  { typeOf_bci = type; }
  
  TFEFunction2D* getPOld() const
  { return p_NSE_old; }
  
  TFEVectFunct2D* getUOld() const 
  { return u_NSE_old; }
  
  double * getSolOld() const { return solution_old->get_entries(); }
  
  std::shared_ptr<TMatrix> getComposedMatForBigSystem()
  {
    if(!mat_aux)
      mat_aux = this->get_matrix().get_combined_matrix();
    return mat_aux;
  }
  void setComposedMatForBigSystem(std::shared_ptr<TMatrix> m) { mat_aux = m; }
  
  std::shared_ptr<TMatrix> getCouplingMatrix() const {return darcyToStokes; }
  void setCouplingMatrix(std::shared_ptr<TMatrix> m) { darcyToStokes = m; }
  
  std::shared_ptr<TMatrix> get_map_sol2eta() const { return map_sol2eta; }
  std::shared_ptr<TMatrix> get_etaToBd() const { return etaToBd; }
  
  void setDirectSol(double *s)
  {
    *solution_direct = s;
  }
  
  BlockVector * get_direct_sol() const
  { return this->solution_direct; }

  TFEVectFunct2D* getUDirect() const
  { return u_NSE_direct; }
  
  TFEFunction2D* getPDirect() const
  { return p_NSE_direct; }
  
  InterfaceFunction*& get_homogeneous_solution() {return eta_hom;}
  
  TCollection* get_collection() const
  { return this->get_velocity_space().GetCollection(); }
  
  void print_assemble_input_output() const
  {
    OutPut("Stokes " << typeOf_bci << "\t");
    OutPut((assemble_on_input ? "" : "no ") << "assembling on input\t");
    OutPut((assemble_on_return ? "" : "no ") << "assembling on output\n");
  }
  void print_normals();
};

#endif  // STOKESPROBLEM_H
