#ifndef AUXILIARYFUNCTIONS_H
#define AUXILIARYFUNCTIONS_H

#include <cmath>

#include <Domain.h>
#include <InnerInterfaceJoint.h>
#include <BlockVector.h>


#include <solver.h>
#include <preconditioner.h>

typedef preconditioner <BlockMatrix, BlockVector> block_prec;
typedef solver <BlockMatrix, BlockVector, block_prec> block_solver;

// OutPut to OutFile and console, depending on verbosity level
#define Out(x,v) { if(TDatabase::ParamDB->SC_VERBOSE > v) OutPut(x);}


struct n_t_vector
{
  double nx;
  double ny;
  double tx;
  double ty;
  void set(double NX, double NY, double TX, double TY)
  {
    nx = NX; 
    ny = NY; 
    tx = TX; 
    ty = TY;
  };
  void set(double NX, double NY)
  {
    nx = NX; 
    ny = NY; 
    tx = -NY; 
    ty = NX;
  };
  n_t_vector(){ nx = 0; ny = 0; tx = 0; ty = 0;}
  n_t_vector(double NX, double NY){ nx = NX; ny = NY; tx = -NY; ty = NX;}
  ~n_t_vector(){ /*OutPut("destructor of n_t_vector\n");*/ }
  friend  bool operator== (n_t_vector &lhs,n_t_vector &rhs)
  {
    double diff = std::abs(lhs.nx - rhs.nx) + std::abs(lhs.ny - rhs.ny)
                  + std::abs(lhs.tx - rhs.tx) + std::abs(lhs.ty - rhs.ty);
    //bool equal = (lhs.nx == rhs.nx && lhs.ny == rhs.ny && lhs.tx == rhs.tx
    //    && lhs.ty == rhs.ty); // this produces false positives
    return diff < 1e-14;  
  }
  friend  bool operator!= (n_t_vector &lhs,n_t_vector &rhs)
  { return !(lhs==rhs); }
  friend std::ostream& operator<<(std::ostream& out, const n_t_vector n);
  double norm() const { return sqrt(nx*nx + ny*ny); }
};

struct test_ansatz_function
{
  double u, ux, uy; // ansatz function and first two derivatives
  double v, vx, vy; // test function and first two derivatives
  // ansatz and test function for vector valued basis functions, derivatives
  // are then not needed (yet)
  double u2,v2; 
  test_ansatz_function(double U, double UX, double UY, double V, double VX,
                       double VY) :
      u(U), ux(UX), uy(UY), v(V), vx(VX), vy(VY), u2(0.0), v2(0.0)
  {
  };
  test_ansatz_function(){u = ux = uy = v = vx = vy = u2 = v2 = 0.0;}
  void setTest(double V, double VX, double VY, double V2 = 1e10)
  {
    v = V;
    vx = VX;
    vy = VY;
    v2 = V2;
  }
  void setAnsatz(double U, double UX, double UY, double U2 = 1e10)
  {
    u = U;
    ux = UX;
    uy = UY;
    u2 = U2;
  }
};

struct local_matrices
{
  /** @brief a vector of matrices
   * 
   * the number of matrices (size of this vector) is usually 9. Each entry is 
   * a map from row-indices to a map. The second map maps column indices to 
   * entries (doubles). 
   */
  std::vector< std::map<int, std::map<int,double> > >m;
  local_matrices() : m(std::vector< std::map<int, std::map<int,double> > >(9))
  {
    //m = std::vector< std::map<int, std::map<int,double> > >(9);
  }
  void add(std::map<int, std::map<int,double> > a){ m.push_back(a); }
  void add(std::map<int, std::map<int,double> > a11, 
           std::map<int, std::map<int,double> > a12,
           std::map<int, std::map<int,double> > a21, 
           std::map<int, std::map<int,double> > a22,
           std::map<int, std::map<int,double> > c,
           std::map<int, std::map<int,double> > b1,
           std::map<int, std::map<int,double> > b2,
           std::map<int, std::map<int,double> > b1t, 
           std::map<int, std::map<int,double> > b2t)
  {
    m.push_back(a11);
    m.push_back(a12);
    m.push_back(a21);
    m.push_back(a22);
    m.push_back(c);
    m.push_back(b1);
    m.push_back(b2);
    m.push_back(b1t);
    m.push_back(b2t);
  }
  
  void info(size_t verbose) const;
};

/** @brief store data to do the local assembling 
 *  
 * This is done for one edge, one quadrature point, one ansatz dof, one test 
 * dof. 
 */
struct local_edge_assembling
{
  n_t_vector nt;
  test_ansatz_function at;
  local_matrices m;
  double hE; // length of edge
  double qw; // quadrature weight
  int ansatzDOF;
  int testDOF;
};


enum ProblemType{Neumann_Neumann, Robin_Robin, Robin_RobinJS,
                 Dirichlet_DirichletSTAB,
                 Dirichlet_Dirichlet,
                 Neumann_Dirichlet, Dirichlet_Neumann};

// boundary condition on interface (similar to BoundCond in "Constants.h"
enum InterfaceCondition{Neumann, Dirichlet, Robin, weakRobin, DirichletSTAB};

// printing the enumerations with the correct name rather than a number
std::ostream& operator<<(std::ostream& out, const ProblemType value);
std::ostream& operator<<(std::ostream& out, const InterfaceCondition value);


/**
 * The intention of this function is to ensure that the current set of 
 * parameters in the database is valid, i.e. that the program will run without
 * errors after this function has terminated.
 * 
 * If parameters in the database are not valid, the program will exit or correct
 * the value, giving a warning. 
 * 
 * A second intention is to write out a little text describing the used 
 * parameters. E.g it tells which methods were chosen and such things.
 */
void check_all_parameters();



/** Description of the function GetInnerInterfaceJoints:
  Given a domain, this function counts the number of all joints which are 
  marked as 'InnerInterfaceJoint'. Such a joint is characterized by different
  reference IDs of the two adjacent cells.
  
  This function should be called on the coarsest grid first because it is 
  rather slow. Then it has to be called as many times as there are refinements
  to get a valid vector called "interface" of interface joints.
    
  If this is the first time this function is called a new vector "interface"
  will be created. 
  If the domain has been refined once, there already exists such a vector, but 
  it is outdated. This function updates this vector using the children of each 
  interface joint. This doubles the length of the vector, because every joint
  has two children.
*/
void GetInnerInterfaceJoints(
    std::vector< const TInnerInterfaceJoint* >& interface, const TDomain& 
Domain);

/** compute |x-y|/|x| in the discrete l^2 norm */
double relDiff(int length, const double* x, const double* y);
/** compute |x-y| in the discrete l^2 norm */
double relabsDiff(int length, const double * x, const double * y);

void WriteMat(char* filename, TMatrix* m);
void WriteVec(char* filename, double *v, int n);


/** ************************************************************************* */

/*
 * get the normal vector to this interface edge "iEdge" pointing outward of 
 * this cell. If you want the normal pointing into this cell, set "outerNormal"
 * to false. With this function you can get the tangent vector as 
 * (tx,ty)=(-ny,nx)
 */
void getNormal(TBaseCell *cell, const TInnerInterfaceJoint * iEdge, 
               double &nx, double &ny, bool outerNormal=true);
/* for a given vertex, find all edges in a collection which share this vertex. 
 * The vector 'edges' must be of length 0. */
void get_edges_of_vertex(const TCollection * const coll,
                         const TVertex * const v, std::vector<TJoint *>& edges);
/* for a given edge, find first (i=0) or second (i=1) vertex. */ 
TVertex* get_vertex_of_edge(const TJoint *const edge, int i);




// y+=m*x
// x has length m->GetN_Columns(), y has length m->GetN_Rows()
void AddMatrixVectorProduct(TMatrix* m,double *x, double *y);
// copy the i-th row of 'm' to the map 'target_index'-th row of 'rows', only 
// non-zero entries are copied. if 'target_index==-1', 'target_index=i is used
void copy_row(int i, const TMatrix* m,
              std::map<int, std::map<int,double>> &rows, double a,
              int target_index = -1);

#endif // AUXILIARYFUNCTIONS_H
