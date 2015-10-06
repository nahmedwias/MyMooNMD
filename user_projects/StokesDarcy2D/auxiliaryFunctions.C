#include <algorithm>
#include <stdlib.h>
#include <auxiliaryFunctions.h>
#include <LinAlg.h>
#include <InterfaceFunction.h>

using namespace std;

/** ************************************************************************ */
void local_matrices::info(size_t verbose) const
{
  // the member 'm' has either 1 or nine matrices 
  OutPut("local_matrices object with " << this->m.size() << " matrices\n");
  for(auto it = this->m.begin(), end = this->m.end();
      it != end && it->size() != 0; ++it)
  {
    OutPut("local matrix " << std::distance(this->m.begin(), it) << endl);
    for(auto i = it->begin(), i_end = it->end(); i != i_end; ++i)
    {
      for(auto j = i->second.begin(), j_end = i->second.end(); j != j_end; ++j)
      {
        OutPut(i->first << "  " << j->first << "\t" << j->second << endl);
      }
    }
  }
}

/** ************************************************************************ */
std::ostream& operator<<(std::ostream& out, const ProblemType value)
{
  const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
  switch(value)
  {
    PROCESS_VAL(Neumann_Neumann);     
    PROCESS_VAL(Robin_Robin);     
    PROCESS_VAL(Robin_RobinJS);
    PROCESS_VAL(Dirichlet_DirichletSTAB);
    PROCESS_VAL(Dirichlet_Dirichlet);
    PROCESS_VAL(Neumann_Dirichlet);
    PROCESS_VAL(Dirichlet_Neumann);
    default: s = "unkown problem type"; break;
  }
#undef PROCESS_VAL
  return out << s;
}

/** ************************************************************************ */
std::ostream& operator<<(std::ostream& out, const InterfaceCondition value)
{
  const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
  switch(value)
  {
    PROCESS_VAL(Neumann);     
    PROCESS_VAL(Dirichlet);     
    PROCESS_VAL(Robin);
    PROCESS_VAL(weakRobin);
    PROCESS_VAL(DirichletSTAB);
    default: s = "unkown interface condition"; break;
  }
#undef PROCESS_VAL
  return out << s;
}

/** ************************************************************************ */
std::ostream& operator<<(std::ostream& out, const n_t_vector n)
{
  //out << "n = (" << n.nx << "," << n.ny << "),  t = (" << n.tx << "," << n.ty 
  //    << ")";
  //return out;
  return out << "n = (" << n.nx << "," << n.ny << ")";
}

/** ************************************************************************ */
void check_all_parameters()
{
  const int pt = TDatabase::ParamDB->StoDa_problemType;
  const int updating_strategy = TDatabase::ParamDB->StoDa_updatingStrategy;
  const int icond = TDatabase::ParamDB->StoDa_interfaceType;
  const int solution_strategy = TDatabase::ParamDB->StoDa_solutionStrategy;
  
  const int stokes_first = TDatabase::ParamDB->StoDa_StokesFirst;
  if(stokes_first != 0 && stokes_first != 1)
    ErrThrow("set 'StoDa_StokesFirst' to either 0 or 1");
  
  if(icond != 0 && icond != 1)
    ErrThrow("Set 'StoDa_interfaceType' to either 0 (Beavers-Joseph-Saffman) "
             + "or 1 (zero tangential velocity) to control the tangential "
             + "condition in the (Navier-) Stokes subdomain");
  
  if(TDatabase::ParamDB->StoDa_weakGamma == 0.0)
    ErrThrow("set 'StoDa_weakGamma' to some positive number for weak "
           + "(Dirichlet) interface conditions (tangential and normal "
           + "direction). Set it to a negative number to enforce the normal "
           + "condition exactly. For Dirichlet problems on the interface in " 
           + "that case also the tangential component will be set exactly if "
           + "'StoDa_interfaceType' is set to 1 (u.t = 0).");
  
  OutPut("\n\n\nsolving a Stokes--Darcy problem in 2D\n");
  OutPut("The coupled problem is a ");
  switch(pt)
  {
    case 0:
      OutPut("Neumann-Neumann-");
      TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
      break;
    case 1:
      OutPut( "Robin-Robin-");
      TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
      break;
    case 2:
      OutPut("weak Robin-Robin-method\n");
      ErrMsg("This is not yet correctly implemented!");
      exit(0);
      break;
    case 3:
      OutPut("stabilized Dirichlet-Dirichlet-method\n");
      ErrMsg("This is not yet correctly implemented!");
      exit(0);
      break;
    case 4:
      OutPut("Dirichlet-Dirichlet-method\n");
      break;
    default:
      ErrMsg("so far only Neumann-Neumann(0) and Robin-Robin(1) is supported");
      exit(0);
      break;
  }
  OutPut("coupled problem with\n");
  OutPut(
      (icond == 0 ? "Beavers-Joseph-Saffmann condition" : "zero tangential component"));
  OutPut(" on the interface\n");
  
  switch(solution_strategy)
  {
    case 0:
      OutPut("No iteration, the coupled system is solved directly\n")
      break;
    case 1: // solve coupled system iteratively using Jacobi or Gauss-Seidel
    case -1:
      OutPut("The coupled system will be solved iteratively using a\n");
      OutPut((TDatabase::ParamDB->StoDa_algorithm==1?"Gauss--Seidel":"Jacobi"));
      OutPut(" type method");
      if(TDatabase::ParamDB->StoDa_algorithm == 1)
      {
        OutPut(" solving ");
        OutPut((stokes_first ? "Stokes" : "Darcy"));
        OutPut(" first\n");
      }
      else
        OutPut(endl);
      if(updating_strategy == 1)
      {
        if(solution_strategy == 1 && pt == 1)
        {
          ErrMsg("trying to solve a regular Robin-Robin problem.\nThere is no " 
                 << "(big) solution of the coupled system avaiable.\nPlease " 
                 << "choose StoDa_solutionStrategy: -1.");
          exit(0);
        }
      }
      if(updating_strategy == 3)
      {
        if(pt != 0) // Neumann-Neumann
        {
          ErrMsg("For C-RR choose StoDa_problemType 0 (Neumann-Neumann)");
          exit(0);
        }
        else 
          OutPut("The C-RR method is employed\n");
      }
      else if(updating_strategy == 4)
      {
        // D-RR
        if(pt != 1) // Robin-Robin
        {
          ErrMsg("For D-RR choose StoDa_problemType 1 (Robin-Robin)");
          exit(0);
        }
        else
          OutPut("The D-RR method is employed\n");
      }
      else if(updating_strategy != 1 && updating_strategy != 2)
      {
        ErrMsg("unknown 'StoDa_updatingStrategy'");
        exit(0);
      }
      break;
    case 2: // solve fixed point formulation
    case -2:
      OutPut("The fixed point equation will be solved iteratively");
      
      if(updating_strategy == 1)
      {
        OutPut(endl);
        if(solution_strategy == 2 && pt == 1)
        {
          ErrMsg("trying to solve a regular Robin-Robin problem.\nThere is no "
                 << "(big) solution of the coupled system avaiable.\nPlease " 
                 << "choose StoDa_solutionStrategy: -2.");
          exit(0);
        }
      }
      else if(updating_strategy == 2)
      {
        OutPut(" using a Newton iteration\n");
        ErrMsg("Newton for the fixed point equation not yet implemented");
        exit(0);
      }
      else if(updating_strategy == 4)
      {
        OutPut(endl);
      }
      else
      {
        OutPut(endl);
        ErrMsg("unknown 'StoDa_updatingStrategy'");
        exit(0);
      }
      break;
    case 3: // solve Stecklov-Poincare formulation
    case -3:
      OutPut("The Stecklov-Poincare equation will be solved using a\n");
      switch(updating_strategy)
      {
        case 1:
          OutPut("preconditioned Richardson iteration\n");
          break;
        case 2:
          OutPut("preconditioned conjugate gradient (cg) method\n");
          break;
        case 3:
          OutPut("preconditioned conjugate gradient square (cgs) method\n"); 
          break;
        case 4:
          OutPut("preconditioned generalized minimal residual (gmres) ");
          OutPut("method\n");
          break;
        case 5:
          OutPut("preconditioned bi-conjugate gradient (BiCg) method\n");
          ErrMsg("Newton for the Stecklov-Poincare equation not implemented");
          exit(0);
          break;
        case 6: // Newton
          OutPut("Newton iteration\n");
          ErrMsg("Newton for the Stecklov-Poincare equation not implemented");
          exit(0);
          break;
        default:
          ErrMsg("unknown 'StoDa_updatingStrategy'");
          exit(0);
          break;
      }
      if(pt != 0)
      {
        ErrMsg("only Neumann-Neumann problems are supported so far, set\n" <<
               "'StoDa_problemType' to 0.");
        exit(0);
      }
      if(TDatabase::ParamDB->StoDa_weakGamma > 0.0)
      {
        ErrMsg("WARNING: weak gamma positive trying to solve a " << 
               "Stecklov-Poincare equation. You probably want to set it " <<
               "negative.");
      }
      if(TDatabase::ParamDB->StoDa_theta_f != 1.0 && stokes_first)
      {
        ErrMsg("WARNING: damping activated for Stecklov-Poincare. Is this what"
               << " you want?");
      }
      if(TDatabase::ParamDB->StoDa_theta_p != 1.0 && !stokes_first)
      {
        ErrMsg("WARNING: damping activated for Stecklov-Poincare. Is this what"
               << " you want?");
      }
      break;
    default:
      ErrMsg("unsupported 'StoDa_solutionStrategy'");
      exit(0);
      break;
  }
  if(solution_strategy < 0)
    OutPut("The coupled system is not solved directly\n");
  
  OutPut("\n\n");
}

/** ************************************************************************ */
void GetInnerInterfaceJoints(vector<TInnerInterfaceJoint *>& interface,
                             const TDomain & Domain)
{
  if(interface.empty())
  {
    // no such vector has been initialzed so far --> make new vector.
    
    // find interface joints (cheap on coarsest mesh)
    // first loop: find number of interface joints
    // second loop: fill vector of interface joints
    
    // using the entire domain will lead to double appearences of every 
    // interface edge in the list of interface joints. We only use one part of 
    // the domain. We could as well use the other subdomain
    TCollection *coll_1 = Domain.GetCollection(It_Finest, 0, 1);
    const int N_Cells = coll_1->GetN_Cells();
    // loop over all cells of this subdomain
    int N_iJoints = 0;
    for(int i = 0; i < N_Cells; i++)
    {
      TBaseCell *cell = coll_1->GetCell(i);
      const int N_Joints = cell->GetN_Edges();
      // loop over all edges of this cell
      for(int j = 0; j < N_Joints; j++)
      {
        if(cell->GetJoint(j)->GetType() == InnerInterfaceJoint)
        { // this edge is an interface joint, this was set in ReadGeo.C
          N_iJoints++;
        }
      }
    }
    // possible improvment: Since we know how often the interface will be 
    // refined, we could reserve enough memory for all interface edges on the
    // finest grid here. Then no reallocation/copy is necessary. However this
    // vector is rather short and therefore reallocation/ copy is cheap.
    interface.reserve(N_iJoints);
    // loop over all cells of this subdomain, fill vector "interface"
    for(int i = 0; i < N_Cells; i++)
    {
      TBaseCell *cell = coll_1->GetCell(i);
      const int N_Joints = cell->GetN_Edges();
      // loop over all edges of this cell
      for(int j = 0; j < N_Joints; j++)
      {
        TJoint *joint = cell->GetJoint(j);
        if(joint->GetType() == InnerInterfaceJoint)
        {
          interface.push_back((TInnerInterfaceJoint*) joint);
        }
      }
    }
    delete coll_1;
  }
  else
  {
    // update vector of interface joints after refinement.
    // every edge has 2 children
    const int N_iJoints = interface.size();
    interface.reserve(2 * N_iJoints);
    // loop over all (old) interfaces, fill vector 'interface'
    for(int i = 0; i < N_iJoints; i++)
    {
      TInnerInterfaceJoint *iijoint = interface[i];
      interface[i] = iijoint->GetChild(0);
      interface.push_back(iijoint->GetChild(1));
    }
  }
  // maybe sort the interface joints
}

/** ************************************************************************ */
void getNormal(TBaseCell * cell, TInnerInterfaceJoint * edge, double &nx,
               double &ny, bool outerNormal)
{
  edge->GetNormal(nx, ny);
  // now check if this normal points in the right direction:
  const int N_Joints = cell->GetN_Edges(); // number of edges in this cell
  // corrdinates of edge points and one more of this cell
  double x0, x1, x2, y0, y1, y2;
  const int iEdge = edge->GetIndexInNeighbor(cell);
  cell->GetVertex(iEdge)->GetCoords(x0, y0);
  cell->GetVertex((iEdge + 1) % N_Joints)->GetCoords(x1, y1);
  // vertex opposite to this edge (works for quadrilaterals,too) 
  cell->GetVertex((iEdge + 2) % N_Joints)->GetCoords(x2, y2);
  // make sure that normal points outwards
  const bool pointingOutward = (nx * (x2 - x0) + ny * (y2 - y0) < 0);
  if((!pointingOutward && outerNormal) || (pointingOutward && !outerNormal))
  {
    nx *= -1;
    ny *= -1;
  }
}

/** ************************************************************************ */
void get_edges_of_vertex(const TCollection * const coll,
                         const TVertex * const v, std::vector<TJoint *>& edges)
{
  if(edges.size() != 0)
  {
    ErrMsg("please provide an empty vector 'edges'");
    exit(1);
  }
  // loop over all cells
  int n_cells = coll->GetN_Cells();
  for(int i = 0; i < n_cells; i++)
  {
    TBaseCell * cell = coll->GetCell(i);
    int n_edges = cell->GetN_Joints(); // == number of vertices in 2D
    // loop over all vertices/edges
    for(int j = 0; j < n_edges; j++)
    {
      if(cell->GetVertex(j) == v)
      {
        // vertex found on edge j, it belongs to two edges
        edges.push_back(cell->GetJoint(j));
        edges.push_back(cell->GetJoint((j-1+n_edges)%n_edges));
      }
    }
  }
  // remove dublicates:
  // first sort, then remove consecutive dublicates
  sort(edges.begin(), edges.end());
  const vector<TJoint *>::iterator it = unique(edges.begin(), edges.end());
  // resize to new size
  edges.resize(distance(edges.begin(), it));
}

/** ************************************************************************ */
TVertex* get_vertex_of_edge(const TJoint *const edge, int i)
{
  if(!edge)
  {
    ErrMsg("No valid TJoint");
    exit(1);
  }
  if(i != 0 && i != 1)
  {
    ErrMsg("the edge index must be either 0 or 1. You provided " << i);
    exit(1);
  }
  // find the first neighbor of this Joint 
  TBaseCell * neighbor = edge->GetNeighb(0);
  if(!neighbor)
  {
    neighbor = edge->GetNeighb(1);
    if(!neighbor)
    {
      ErrMsg("This joint has no neighbor");
      exit(1);
    }
  }
  
  // find local index of this 'edge' in 'neighbor'
  int loc_index = -1;
  int n_edges = neighbor->GetN_Joints();
  for(int j = 0; j < n_edges; j++)
  {
    if(neighbor->GetJoint(j) == edge)
    {
      loc_index = j;
      break;
    }
  }
  if(loc_index == -1)
  {
    ErrMsg("this edge could not be found in this cell");
    exit(1);
  }
  
  // in 2D the vertex and the edge indices are the same. That means the i-th 
  // local edge has local vertices i and i+1 (modulo number of edges).
  TVertex * v = neighbor->GetVertex((loc_index+i)%n_edges);
  return v;
}

/** ************************************************************************ */
void copy_row(int i, const TMatrix* m,
              std::map<int, std::map<int, double> >& rows,
              double a, int target_index)
{
  if(target_index == -1)
    target_index = i;
  int * rowPtr = m->GetRowPtr();
  int * colPtr = m->GetKCol();
  double * entries = m->GetEntries();
  if(i < 0 || i >= m->GetN_Rows())
  {
    ErrMsg("the specified row does not exist");
    exit(1);
  }
  if(a == 0.0)
    return; // nothing needs to be done
  
  // loop over all entries in the i-th row
  for(int r = rowPtr[i]; r < rowPtr[i+1]; r++)
  {
    if(entries[r] != 0.0)
      (rows[target_index])[colPtr[r]] += a * entries[r];
  }
}

/** ************************************************************************ */
void WriteMat(char* filename, TMatrix* m)
{
  std::ofstream mat;
  mat.open(filename);
  int n_Eqn = m->GetN_Rows();
  int* row = m->GetRowPtr();
  int* col = m->GetKCol();
  double * vals = m->GetEntries();
  for(int i = 0; i < n_Eqn; i++)
  {
    int begin = row[i];
    int end = row[i + 1];
    for(int j = begin; j < end; j++)
      mat << setw(6) << i + 1 << "\t" << col[j] + 1 << "\t" << setw(30)
          << vals[j] << endl;
  }
  mat.close();
}

/** ************************************************************************ */
void WriteVec(char* filename, double *v, int n)
{
  std::ofstream vec;
  vec.open(filename);
  for(int i = 0; i < n; i++)
    vec << setw(6) << i + 1 << "\t" << setw(30) << v[i] << endl;
  vec.close();
}

/** ************************************************************************ */
double relDiff(int length, const double* x, const double *y)
{
  double ret = 0;
  double norm_x = Ddot(length,x,x);
  for(int i=0;i<length;i++)
  {
    ret += POW(x[i]-y[i],2);
  }
  return sqrt(ret/norm_x);
}

/** ************************************************************************ */
double relabsDiff( int length, const double* x, const double* y )
{
  double ret = 0;
  double norm_x = MAX(Ddot(length,x,x),1);
  for(int i=0;i<length;i++)
    ret += POW(x[i]-y[i],2);
  return sqrt(ret/norm_x);
}
