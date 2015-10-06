#ifndef INTERFACEFUNCTION_H
#define INTERFACEFUNCTION_H
class InterfaceFunction;

#include <auxiliaryFunctions.h>
#include <StokesProblem.h>
#include <DarcyProblem.h>
#include <StokesDarcy2D.h>

/** ************************************************************************* */
/** Description of the class InterfacFunction
 * An object of this class is a representation of a function on the interface.
 * This is more or less only implemented in the special way which is needed 
 * here. We can construct such an interface function using a vector of 
 * interface edges (vector<TInnerInterfaceJoint *>). The function then can 
 * be evaluated at any point (on the interface). Other finite element 
 * functions on either one subdomain (TFEFunction2D*) can be interpolated 
 * into an interface function. 
 */
class InterfaceFunction : public BlockVector
{
 private:
  std::vector<TInnerInterfaceJoint *> interface;
  /* code for the used space on the interface
   * 2 : continuous and piecewise quadratic
   * 1 : continuous and piecewise linear   (not implemented, maybe never will)
   * 0 : piecewise constant                (not implemented, maybe never will)
   * -1: piecewise linear                  (not implemented, maybe never will)
   * -2: piecewise quadratic */
  int spaceType;
  /* array determining the degrees of freedom for each edge, i.e. the positions
   * in the array 'blocks[0]' (see BlockVector.h) which corresponds to this 
   * edge. It has length 3*interface.size() for (continuous) piecewise 
   * quadratics. I.e. DOF[3*i+j] tells you the position of the j-th (0<=j<=2) 
   * degree of freedom on the i-th (0<=i<=interface.size()) edge within the 
   * array 'blocks[0]'. */
  std::vector<unsigned int> DOF;
  /* store information on the neighboring structure. For edge i (interface[i]) 
   * its neighbors have indices adjacency[2*i] and adjacency[2*i+1]. We call 
   * the neighbors first and second neighbor. A value of i for the i-th edge 
   * means there is no neighbor. The first neightbor (At adjacency[2*i]) is at 
   * (x,y), where (x,y) are from interface[i]->GetParams(x,y,vecx,vecy); Then 
   * the second neighbor (at adjacency[2*i+1]) is at the other edge point 
   * (x+vecx,y+vecy). */
  std::vector<unsigned int> adjacency;
  /* store normal and tangential vector for each dof.  
   */
  std::vector<n_t_vector> normal_tangential;
  /* depending on which space this functions is supposed to belong to, this
   * vector contains the dofs which are set to 0. These are dofs at the 
   * (Dirichlet) boundary. This is needed to restrict this interface function
   * to the correct interface space (depending on StokesFirst) in case of the
   * Stecklov-Poincare interface equation
   */
  std::vector<unsigned int> *dofs_to_be_restricted;
  
  /* create the arrays 'DOF', 'adjacency', 'normal_tangential' and return 
   * the length of this vector
   */
  unsigned int make_adjacency_and_DOF();
 public:
  /* constructor */
  InterfaceFunction(const std::vector<TInnerInterfaceJoint *>& in, int s);
  
  // empty constructor, used for iterative template solvers
  // The InterfaceFunction is not useful after this constructor alone!
  InterfaceFunction();
  
  // copy constructor, just copy all entries
  InterfaceFunction(const InterfaceFunction& eta);
  
  // constructor, given a certain length (needed for gmres)
  InterfaceFunction(const unsigned int a);
  
  /* destructor */
  ~InterfaceFunction();
  
  /* get the function value at a specified point (x,y), terminates the program
   * if the point is not on the interface. */
  double GetValue(const double x, const double y) const;
  /* get the function value at a specified point (x,y) which is on the i-th 
   * interface edge, throws 'false' if point is not on this edge */
  double GetValueLocal(const int i, const double x, const double y) const;
  /* interpolate a function on the interface and add it to this interface
   * function. The function lives on either one of the subdomains. The 
   * parameter 'comp' determines if the function or certain derivatives 
   * are taken into account. The following values are poosible:
   * 0: value of function 'f' at interface   (default)
   * 1: normal component of gradient of function 'f' at interface
   * 2: tangential component of gradient of function 'f' at interface
   * 3: first component of gradient of function 'f' at interface
   * 4: second component of gradient of function 'f' at interface
   *-1: for Raviart-Thomas elements, normal component of f 
   * 
   * The additional factor is optional.
   * */
  void add(TFEFunction2D *f, int comp=0, double factor=1.0);
  /* interpolate a vector function on the interface and add it to this 
   * interface function. The function lives on either one of the subdomains. 
   * Which component of the vector is used, can be controlled via the 
   * (optional) parameter 'comp'. The following values are possible:
   * 1: normal component        (default)
   * 2: tangential component
   * 3: first component
   * 4: second component
   * 5: n.DD(v).n
   * 
   * The additional factor is optional.
   */
  void add(TFEVectFunct2D* v, int comp=1, double factor=1.0);
  /* add other interface function to this one (multiplied by a factor).  
   * It is checked if 'this' and eta belong to the same space.
   * The additional factor is optional.
   * 
   * "this += factor*eta"
   */
  void add(const InterfaceFunction * const eta, double factor=1.0);
  /* return the normal and tangential vector for a given dof */
  n_t_vector get_normal(unsigned int dof) const
  { return normal_tangential[dof]; }
  
  /* set this function to zero on the boundary of the interface */
  void zeroAtBoundary();
  // restrict correctly according to StokesFirst, see the explanation of the 
  // other restrict method (a template version) for what 'invert' does
  void restrict(const StokesDarcy2D& sd, bool invert = false);
  // restrict to the subspace which is a trace space of the 'subproblem' (either
  // (StokesProblem or DarcyProblem). This replaces the dof values with zeros,
  // who touch a Dirichlet (outer) boundary in the respective subdomain.  
  // If the boolean 'invert' is true, then all values are set to zero, except
  // those dofs touching a Dirichlet boundary.
  // this function uses the vector 'dofs_to_be_restricted' which is created in
  // case it doesn't exist yet. Once it exist, you can use restrict() without 
  // parameters as well
  template < class subproblem >
  void restrict(const subproblem&, bool invert = false);
  void restrict(bool invert = false);
  
  // compute the integral of this InterfaceFunction
  double integral() const;
  
  // Add an appropriate function on the interface to this InterfaceFunction such
  // that the resulting function has integral 'a'. This makes sure that the
  // subproblems are solvable, i.e. that the divergence constraint is not 
  // violated.
  void set_integral(double a);
  
  /* update eta_p
   * Depending on the updating strategy (e.g. C-RR, D-RR, Neumann-Neumann) the 
   * update of the interface function is done. Then the interface function is
   * ready to be passed over to the Darcy part as boundary data.
   */
  void update(StokesProblem &s, DarcyProblem &d,InterfaceFunction *eta_f=NULL);
  /* update  eta_f
   * Depending on the updating strategy (e.g. C-RR, D-RR, Neumann-Neumann) the 
   * update of the interface function is done. Then the interface function is
   * ready to be passed over to the Stokes part as boundary data.
   */
  void update(DarcyProblem &d, StokesProblem &s,InterfaceFunction *eta_p=NULL);
  
  /* Print out all arrays, for testing */
  void PrintInfo(std::string name = "") const;
  /* Print out all x,y,f(x,y) triples */
  void PrintVals(std::string name = "") const;
  /*
   * This is not implemented yet
   * Write a gnuplot file "a.gp". Type "gnuplot a.gp" so see the output. This 
   * interface function is plotted against the parameter t in [0,1] which 
   * parametrises the interface.   
   */
  void PrintGnuplotFile(char* a) const;
  
  int getSpaceType() const { return spaceType; }
  
  bool isBoundaryEdge(unsigned int i) const
  { return (adjacency[2*i]==i || adjacency[2*i+1]==i); }
  
  /* get the DOF of the i-th local dof on 'edge' */
  int getDOF(const int edge, const int i) const { return DOF.at(3*edge + i); }
  
  // operator= more or less only calls the BlockVector-version of operator=
  InterfaceFunction& operator=(const InterfaceFunction& r)
  {
    this->BlockVector::operator=(r);
    if(spaceType == -4711) 
    {
      // empty constructor has been called on this object --> copy all variables
      // or constructor InterfaceFunction(const double a)
      spaceType = r.getSpaceType();
      interface = r.interface;
      DOF = r.DOF;
      adjacency = r.adjacency;
      if(r.dofs_to_be_restricted != NULL)
      {
        dofs_to_be_restricted = new std::vector<unsigned int>();
        *dofs_to_be_restricted = *(r.dofs_to_be_restricted);
      }
      else
        dofs_to_be_restricted = NULL;
    }
    return *this;
  }
  InterfaceFunction& operator=(const double& a)
  {
    this->BlockVector::operator=(a);
    return *this;
  }
};

#endif // INTERFACEFUNCTION_H
