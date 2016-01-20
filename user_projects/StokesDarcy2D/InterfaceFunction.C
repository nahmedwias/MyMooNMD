#include <algorithm>
#include <InterfaceFunction.h>
#include <BoundEdge.h>

/** ************************************************************************ */
InterfaceFunction::InterfaceFunction(
    const std::vector<const TInnerInterfaceJoint*>& in, int s) :
    BlockVector(), interface(in), spaceType(s)
{
  if(spaceType != 2 && spaceType != -2)
    ErrThrow("Constructor for InterfaceFunction: only piecewise quadratic ",
             "functions (continuous or discontinuous) are suported so far!");
  
  DOF.resize((abs(spaceType) + 1) * interface.size(), 0); // initialize to zero
  adjacency.resize(2 * interface.size(), 0); // initialize to zero
  normal_tangential.resize(DOF.size());
  unsigned int length = make_adjacency_and_DOF();
  
  // length of the array values is 2*interface.size()+1 for P2 and
  // 3* interface.size() for P2-discontinuous
  //int length = 2*interface.size()+1;
  //if(spaceType==-2) length = 3*interface.size();
  
  // initialize BlockVector part of this class
  this->BlockVector::entries.resize(length, 0.);
  this->BlockVector::lengths.resize(1, length);
  this->BlockVector::actives.resize(1, length);
  dofs_to_be_restricted = NULL;
}

/** ************************************************************************ */
InterfaceFunction::InterfaceFunction() :
    BlockVector()
{
  interface.resize(0);
  DOF.resize(0);
  adjacency.resize(0);
  normal_tangential.resize(0);
  spaceType = -4711; // some unrealistic value
  dofs_to_be_restricted = NULL;
}

/** ************************************************************************ */
InterfaceFunction::InterfaceFunction(const InterfaceFunction& eta) :
    BlockVector((const BlockVector&)eta)
{
  spaceType = eta.spaceType;
  interface = eta.interface;
  DOF = eta.DOF;
  normal_tangential = eta.normal_tangential;
  adjacency = eta.adjacency;
  if(eta.dofs_to_be_restricted)
  {
    dofs_to_be_restricted = new std::vector<unsigned int>();
    *dofs_to_be_restricted = *(eta.dofs_to_be_restricted); // copy
  }
  else
    dofs_to_be_restricted = NULL;
}

/** ************************************************************************ */
InterfaceFunction::InterfaceFunction( const unsigned int a ) :
    BlockVector(a)
{
  spaceType = -4711;
  interface.resize(0);
  DOF.resize(0);
  normal_tangential.resize(0);
  adjacency.resize(0);
  dofs_to_be_restricted = NULL;
}

/** ************************************************************************ */
InterfaceFunction::~InterfaceFunction()
{
  delete dofs_to_be_restricted;
}

/** ************************************************************************ */
double InterfaceFunction::GetValueLocal(const int i, const double x,
                                        const double y) const
{
  TInnerInterfaceJoint const * const IJoint = interface[i];
  double t = -1;
  // check if the Point (x,y) is on this edge
  // if so, set t such that (x,y) = (x0,y0) + t*(vecx,vecy) 
  {
    double nx, ny;
    IJoint->GetNormal(nx, ny);
    double x0, y0, vecx, vecy;
    IJoint->GetParams(x0, y0, vecx, vecy);
    // length of this edge (squared)
    const double length = vecx * vecx + vecy * vecy;
    // the vector (x,y)-(x0,y0) should be parallel to the edge 
    // so the inner product with the normal should be zero
    if(abs((x - x0) * nx + (y - y0) * ny) > 1.0e-14)
      throw false;
    // the Point (x,y) is on the line defined by this edge. Now check if it  
    // lies between the two edge points (x0,y0) and (x0+vecx,y0+vecy)
    t = (x - x0) * vecx + (y - y0) * vecy;
    if(t < 0.0 || length < t)
      throw false;
    // transform t to the interval [-1,1]
    t *= 2 / length;
    t -= 1.0;
  }
  // Point (x,y) is on this edge!
  
  // degrees of freedom associated to this edge are the three entries 
  // DOF[3*i], DOF[3*i+1], and DOF[3*i+2]
  double value = 0;
  // for piecewise quadratic functions in 1D there are three basis functinos
  // on [-1,1], the basis functions are t*(t-1)/2, (1-t*t), and t*(t+1)/2
  value += this->BlockVector::entries[DOF[3 * i]] * t * (t - 1) / 2;
  value += this->BlockVector::entries[DOF[3 * i + 1]] * (1 - t * t);
  value += this->BlockVector::entries[DOF[3 * i + 2]] * t * (t + 1) / 2;
  
  return value;
}

/** ************************************************************************ */
double InterfaceFunction::GetValue(const double x, const double y) const
{
  // this function is untested!!!!
  // loop over all interface edges, if (x,y) belongs to one then return
  for(unsigned int i = 0; i < interface.size(); i++)
  {
    double a;
    try
    {
      a = this->GetValueLocal(i, x, y);
    }
    catch(bool r)
    {
      continue;
    }
    return a;
  }
  ErrThrow("Value not on the interface!!!!");
}

/** ************************************************************************ */
unsigned int InterfaceFunction::make_adjacency_and_DOF()
{
  if(spaceType != 2 && spaceType != -2)
  {
    ErrThrow("Constructor for InterfaceFunction: only piecewise quadratic ",
             "functions (continuous or discontinuous) are suported so far!");
  }
  const int n_dof_per_edge = abs(spaceType) + 1; // number of dofs per edge
  
  // initilize all values in 'adjacency' such that every edge has no neighbors.
  // then change all entries according to existing neighbors.
  for(unsigned int i = 0; i < interface.size(); i++)
  {
    adjacency[2 * i] = i;
    adjacency[2 * i + 1] = i;
  }
  // finding neighbors
  for(unsigned int i = 0; i < interface.size(); i++)
  {
    // check if both edge points already have a neighbor assigned
    if(adjacency[2 * i] != i && adjacency[2 * i + 1] != i)
      continue;
    // one edge point and the vector pointing to the other edge point
    double x, y, vecx, vecy;
    interface[i]->GetParams(x, y, vecx, vecy);
    for(unsigned int j = 0; j < interface.size(); j++)
    {
      if(i == j)
        continue;
      // one edge point and the vector pointing to the other edge point
      double jx, jy, jvecx, jvecy;
      interface[j]->GetParams(jx, jy, jvecx, jvecy);
      // check if one edge point of the i-th edge is an edge point of the j-th
      if(x == jx && y == jy)
      {
        adjacency[2 * i] = j;
        adjacency[2 * j] = i;
      }
      else if(x == jx + jvecx && y == jy + jvecy)
      {
        adjacency[2 * i] = j;
        adjacency[2 * j + 1] = i;
      }
      else if(x + vecx == jx && y + vecy == jy)
      {
        adjacency[2 * i + 1] = j;
        adjacency[2 * j] = i;
      }
      else if(x + vecx == jx + jvecx && y + vecy == jy + jvecy)
      {
        adjacency[2 * i + 1] = j;
        adjacency[2 * j + 1] = i;
      }
    }
  }
  
  // write normal and tangential vector
  // for a continuous interface function the interface has to be straight 
  for(unsigned int i = 0; i < interface.size(); i++)
  {
    double nx, ny;
    TBaseCell *s_cell = interface[i]->GetNeighbour(0);
    if(s_cell->GetReference_ID() != 2) // assuming Stokes-subdomain has id 2
      s_cell = interface[i]->GetNeighbour(1);
    getNormal(s_cell, interface[i], nx, ny);
    for(int j = 0; j < n_dof_per_edge; j++)
      normal_tangential[n_dof_per_edge * i + j].set(nx, ny);
  }
  
  // count the number of degrees of freedom (DOF)
  unsigned int countDOF = 0;
  // index of a neighboring edge
  if(spaceType == 2) // continuous space
  {
    for(unsigned int i = 0; i < interface.size(); i++)
    {
      // first DOF:
      // check if this DOF has been assigned from the neighboring edge already
      // which can only happen if the normals on both edges coincide
      unsigned int j = adjacency[2 * i];
      if(j < i)
         //&& normal_tangential[n_dof_per_edge * i]
         //   == normal_tangential[n_dof_per_edge * j])
      {
        // this dof has been assigned from the neighboring edge already
        // now it is either the first (offset=0) or the third (offset=2) dof in  
        // this neighboring edge
        unsigned int offset = 0;
        if(adjacency[2 * j + 1] == i)
          offset = n_dof_per_edge - 1; // == abs(spaceType)
        DOF[n_dof_per_edge * i] = DOF[n_dof_per_edge * j + offset];
      }
      else
        DOF[n_dof_per_edge * i] = countDOF++;
      // second DOF:
      DOF[n_dof_per_edge * i + 1] = countDOF++;
      // third DOF:
      // check if this DOF has been assigned from the neighboring edge already
      j = adjacency[2 * i + 1];
      if(j < i)
         //&& normal_tangential[n_dof_per_edge * i + 2]
         //   == normal_tangential[n_dof_per_edge * j])
      {
        // this dof has been assigned from the neighboring edge already
        // now it is either the first (offset=0) or the third (offset=2) dof in  
        // this neighboring edge 
        unsigned int offset = 0;
        if(adjacency[2 * j + 1] == i)
          offset = n_dof_per_edge - 1; // == abs(spaceType)
        DOF[n_dof_per_edge * i + 2] = DOF[n_dof_per_edge * j + offset];
      }
      else
        DOF[n_dof_per_edge * i + 2] = countDOF++;
    }
  }
  else if(spaceType == -2) // discontinuous space (three dofs on each edge)
  {
    for(unsigned int i = 0; i < interface.size(); i++)
    {
      for(int j = 0; j < n_dof_per_edge; j++)
        DOF[n_dof_per_edge * i + j] = countDOF++;
    }
  }
  else
  {
    ErrThrow("InterfaceFunction::make_adjacency_and_DOF(): unsopported space ",
             "type!");
  }
  return countDOF;
}

/** ************************************************************************ */
void InterfaceFunction::PrintInfo(std::string name) const
{
  // number of degrees of freedom on each edge:
  unsigned int n = abs(spaceType) + 1; //3;
                   
  Output::print<1>("\n");
  if(name != "")
    Output::print<1>("Information on InterfaceFunction ", name);
  Output::print<1>("Number of edges ", interface.size(), ", number of DOFs ",
                   this->length());
  Output::print<1>("\nAdjacency:");
  for(unsigned int i = 0; i < interface.size(); i++)
  {
    Output::print<1>(i, ":\t", adjacency[2*i], "\t", adjacency[2*i+1]);
  }
  Output::print<1>("\nDOFs");
  std::ostringstream os;
  for(unsigned int i = 0; i < interface.size(); i++)
  {
    os << i << ":\t";
    for(unsigned int j = 0; j < n; j++)
      os << DOF[n*i+j] << "\t";
    os << "\n";
  }
  Output::print<1>(os.str());
  Output::print<1>("\nnormals");
  for(unsigned int i = 0; i < interface.size(); i++)
  {
    for(unsigned int j = 0; j < n; j++)
    {
      Output::print<1>("edge ", i, " loc_dof ", j, " :\t",
                       normal_tangential[n*i+j]);
    }
  }
  Output::print<1>("\nValues");
  //for(unsigned int i = 0; i < values.size(); i++)
  for(unsigned int i = 0; i < this->length(); i++)
    Output::print<1>(i, "\t", this->BlockVector::entries[i]);
}

/** ************************************************************************ */
void InterfaceFunction::PrintVals(std::string name) const
{
  // number of degrees of freedom on each edge:
  unsigned int n = abs(spaceType) + 1; //3
  if(name == "")
  {
    Output::print<1>("\nValues");
  }
  else
  {
    Output::print<1>("\nValues of ", name);
  }
  for(unsigned int i = 0; i < interface.size(); i++)
  {
    double x, y, vecx, vecy;
    interface[i]->GetParams(x, y, vecx, vecy);
    // do not print values at the two edge point twice for continuous 
    // (spaceType>0) interface functions 
    if(adjacency[2 * i] >= i || spaceType < 0)
      Output::print<1>(x, "\t", y, "\t", this->BlockVector::entries[DOF[n*i]]);
    Output::print<1>(x + 0.5*vecx, "\t", y + 0.5*vecy, "\t",
                     this->BlockVector::entries[DOF[n*i+1]]);
    if(adjacency[2 * i + 1] >= i || spaceType < 0)
      Output::print<1>(x + vecx, "\t", y + vecy, "\t",
                       this->BlockVector::entries[DOF[n*i+2]]);
  }
  Output::print<1>("ValuesEnd");
}

/** ************************************************************************ */
void InterfaceFunction::PrintGnuplotFile(char* a) const
{
  std::ofstream FILE(a);
  
  // how many intervals are exported for each interface edge 
  const unsigned int n_sub_intervals = 10; // choose an even number
  FILE << "# use this file with gnuplot via: gnuplot <file_name>";
  FILE << "# number of subintervals = " << n_sub_intervals << "\n";
  FILE << "# number of interface edges = " << interface.size() << "\n";
  
  //FILE << "set xrange [0:1]\n";
  FILE << "set title \""<< a << "\"\n";
  FILE << "unset key\n";
  FILE << "plot \"-\" using 1:2 w lines\n" << setprecision(12);
  
  // compute length of interface:
  double length = 0.0;
  for(unsigned int i = 0; i < interface.size(); i++)
    length += interface[i]->GetLength();
  // find begining of interface (lower left)
  int startingEdge = -1;
  // starting left or right (first or third dof)
  bool start_left = true;
  // coordinates of vertex which is furthest to the left of all vertices on 
  // the interface. If there are more than one, than the vertex with the 
  // smallest y-value is chosen. 
  double x_min = 1e10, y_min = 1e10;
  for(unsigned int i = 0; i < interface.size(); i++)
  {
    double x, y, vecx, vecy;
    interface[i]->GetParams(x, y, vecx, vecy);
    if(x < x_min)
    {
      x_min = x;
      y_min = y;
      startingEdge = i;
      start_left = true;
    }
    if(x + vecx < x_min)
    {
      x_min = x + vecx;
      y_min = y + vecy;
      startingEdge = i;
      start_left = false;
    }
    if(x == x_min && y < y_min)
    {
      y_min = y;
      startingEdge = i;
      start_left = true;
    }
    if(x + vecx == x_min && y + vecy < y_min)
    {
      y_min = y + vecy;
      startingEdge = i;
      start_left = false;
    }
  }
  
  // loop over all interface edges. find parameter t. get value there, print
  unsigned int currentEdge = startingEdge;
  double t = 0.0; // starting parameter t in [0,1]
  
  FILE << "\t" << t << "\t"
       << this->BlockVector::entries[DOF[3*currentEdge + (start_left ? 0 : 2)]] 
       << endl;
  while(t < 1.0)
  {
    // increment of t for this edge
    double t_incr = interface[currentEdge]->GetLength() / length;
    // loop over all points on this edge which are printed to the gnuplot file
    // n_sub_intervals+1 points on each edge, n_sub_intervals pieces
    for(unsigned int i = 1; i <= n_sub_intervals; i++) 
    {
      double s = -1.0 + i*2.0/n_sub_intervals; // local parameter on reference edge [-1,1]
      // degrees of freedom associated to this edge are the three entries 
      // DOF[(abs(spaceType)+1)*currentEdge  ], 
      // DOF[(abs(spaceType)+1)*currentEdge+1], and 
      // DOF[(abs(spaceType)+1)*currentEdge+2]
      double value = 0;
      const unsigned int offset = start_left ? 0 : 2;
      // for piecewise quadratic functions in 1D there are three basis functinos
      // the basis functions are t*(t-1)/2, (1-t*t), and t*(t+1)/2
      value += this->entries[DOF[3 * currentEdge + offset]] * s * (s - 1) / 2;
      value += this->entries[DOF[3 * currentEdge + 1]] * (1 - s * s);
      value += this->entries[DOF[3 * currentEdge + 2 - offset]] * s * (s + 1)/2;
      FILE << "\t" << t + i*t_incr/n_sub_intervals << "\t" << value << endl;
    }
    t += t_incr;
    // find next edge
    unsigned int next_edge = adjacency[2 * currentEdge + (start_left ? 1 : 0)];
    if(currentEdge == next_edge)
    {
      if(t < 1.0 - t_incr/(2*n_sub_intervals))
      {
        ErrThrow("no neighbor here ", +t);
      }
      else
        break;
    }
    if(adjacency[2*next_edge + 1] == currentEdge)
      start_left = false;
    currentEdge = next_edge;
    //FILE << "\n"; // make one line for each interval (small skips in between)
  }
  FILE << "\te\n";
  FILE << "pause -1 \n";//\"Hit any key to continue\"\n";
  FILE.close();
}

/** ************************************************************************ */
void InterfaceFunction::add(TFEFunction2D *f, int comp, double factor)
{
  // number of degrees of freedom on each edge:
  const unsigned int n = 3; //abs(spaceType)+1; 
  // id of the space on which 'f' is defined. All cells have same id
  const int ID =
      f->GetFESpace2D()->GetCollection()->GetCell(0)->GetReference_ID();
  // coordinates of first edge point and vector pointing to the second
  double x, y, vecx, vecy;
  // normal and tangential component of the edge (computed only if needed)
  double nx, ny, tx, ty;
  for(unsigned int i = 0; i < interface.size(); i++)
  {
    // get first ege point and vector pointing to the second
    interface[i]->GetParams(x, y, vecx, vecy);
    // get cell of this edge in the feSpace of 'f'
    TBaseCell *cell = interface[i]->GetNeighbour(0);
    if(cell->GetReference_ID() != ID)
      cell = interface[i]->GetNeighbour(1);
    // index of this cell in feSpace of 'f'
    const int cellInd = cell->GetCellIndex();
    // get normal (only if it is needed)
    if(comp == 1 || comp == 2 || comp == -1)
      getNormal(cell, interface[i], nx, ny);
    if(comp == 2)
    { // we take the normal and rotate it counterclockwise
      tx = -ny;
      ty = nx;
    }
    // evaluate 'f' and if necessary its first derivatives at both edge points 
    // and midpoint (needed for quadratic functions on the interface)
    // loop over these three points (degrees of freedom) on this edge
    for(unsigned int j = 0; j < 3; j++)
    {
      // local factor which is used to projecto onto a continuous space
      // this results in a averaging of a (possibly discontinuous) function f 
      double locfactor = factor;
      // the current point has coordinates (loc_x , loc_y)
      const double loc_x = x + j * 0.5 * vecx;
      const double loc_y = y + j * 0.5 * vecy;
      // this degree of freedom is also a degree of freedom in a neighbor
      if(j == 0 && adjacency[2 * i] != i && spaceType > 0)
        locfactor *= 0.5;
      if(j == 2 && adjacency[2 * i + 1] != i && spaceType > 0)
        locfactor *= 0.5;
      // value of 'f' at required points
      // actually we need an array, but only the first entry is accessed in 
      // f->FindValueLocal(...), so the address of this double value works, too
      double fval = 0.0;
      if(comp == 0) // evaluate f at the current point and store value in "fval"
        f->FindValueLocal(cell, cellInd, loc_x, loc_y, &fval);
      else if(comp == -1)
      { // normal component for Raviart-Thomas elements
        double vals[2];
        f->FindValueLocal(cell, cellInd, loc_x, loc_y, vals);
        fval = vals[0] * nx + vals[1] * ny;
      }
      else
      {
        // array to store function value and both first derivatives 
        double vals[3];
        // evaluate f and its derivatives and store them in the array "vals"
        f->FindGradientLocal(cell, cellInd, loc_x, loc_y, vals);
        // calulate the corresponding component according to "comp"
        switch(comp)
        {
          case 1:
            fval = vals[1] * nx + vals[2] * ny;
            break;
          case 2:
            fval = vals[1] * tx + vals[2] * ty;
            break;
          case 3:
            fval = vals[1];
            break;
          case 4:
            fval = vals[2];
            break;
          default:
            ErrThrow("InterfaceFunction::add: unsupported component. Exiting");
            break;
        }
      }
      this->BlockVector::entries[DOF[n * i + j]] += locfactor * fval;
    }
  }
}

/** ************************************************************************ */
void InterfaceFunction::add(TFEVectFunct2D* v, int comp, double factor)
{
  // number of degrees of freedom on each edge:
  const unsigned int n = 3; //abs(spaceType)+1; 
  // the two components of this vector valued function
  TFEFunction2D* v1 = v->GetComponent(0);
  TFEFunction2D* v2 = v->GetComponent(1);
  // id of the space on which 'v' is defined. All cells have same id
  const int ID = v1->GetFESpace2D()->GetCollection()->GetCell(0)
      ->GetReference_ID();
  // coordinates of first edge point and vector pointing to the second
  double x, y, vecx, vecy;
  // normal and tangential component of the edge
  double nx, ny, tx, ty;
  // value of the two components of 'v' at required points
  // actually we need an array, but only the first entry is accessed in 
  // v1->FindValueLocal(...), so that the address of this double value works,
  // too
  double val1, val2;
  // temporary value which will be added into the values-array of this 
  // ineterface function
  double temp = 0;
  // derivatives in x and y direction for the two components of v
  double der1[3], der2[3];
  for(unsigned int i = 0; i < interface.size(); i++)
  {
    // get cell of this edge in the feSpace of 'v'
    TBaseCell *cell = interface[i]->GetNeighbour(0);
    if(cell->GetReference_ID() != ID)
      cell = interface[i]->GetNeighbour(1);
    // index of this cell in feSpace of 'v'
    const int cellInd = cell->GetCellIndex();
    // get first ege point and vector pointing to the second
    interface[i]->GetParams(x, y, vecx, vecy);
    // get normal (only if it is needed)
    if(comp == 1 || comp == 2)
      getNormal(cell, interface[i], nx, ny);
    if(comp == 2)
    { // we take the normal and rotate it counterclockwise
      tx = -ny;
      ty = nx;
    }
    
    // evaluate 'v' and if necessary its first derivatives at both edge points 
    // and midpoint (needed for quadratic functions on the interface)
    // loop over these three points (degrees of freedom) on this edge
    for(unsigned int j = 0; j < 3; j++)
    {
      // local factor which is used to projecto onto a continuous space
      // this results in a averaging of a (possibly discontinuous) function v 
      double locfactor = factor;
      // the current point has coordinates (loc_x , loc_y)
      const double loc_x = x + j * 0.5 * vecx;
      const double loc_y = y + j * 0.5 * vecy;
      // this degree of freedom is also a degree of freedom in a neighbor
      if(j == 0 && adjacency[2 * i] != i && spaceType > 0)
        locfactor *= 0.5;
      if(j == 2 && adjacency[2 * i + 1] != i && spaceType > 0)
        locfactor *= 0.5;
      
      if(comp != 4) // no need to evaluate 1st component if only 2nd is required
        v1->FindValueLocal(cell, cellInd, loc_x, loc_y, &val1);
      if(comp != 3) // no need to evaluate 2nd component if only 1st is required
        v2->FindValueLocal(cell, cellInd, loc_x, loc_y, &val2);
      if(comp == 5)
      {
        v1->FindGradientLocal(cell, cellInd, loc_x, loc_y, der1);
        v2->FindGradientLocal(cell, cellInd, loc_x, loc_y, der2);
      }
      // calulate the corresponding component according to "comp"
      switch(comp)
      {
        case 1: // normal component
          temp = val1 * nx + val2 * ny;
          break;
        case 2: // tangential component
          temp = val1 * tx + val2 * ty;
          break;
        case 3: // first component
          temp = val1;
          break;
        case 4: // second component
          temp = val2;
          break;
        case 5: // n.DD(v).n
          temp = nx * nx * der1[1] + nx * ny * (der1[2] + der2[1])
                 + ny * ny * der2[2];
          break;
        default:
          ErrThrow("InterfaceFunction::add(...): The parameter comp ",
                   "should be either 1,2,3,4, or 5. You have chosen ", comp,
                   ". Exiting");
      }
      this->BlockVector::entries[DOF[n * i + j]] += locfactor * temp;
    }
  }
  delete v1;
  delete v2;
}

/** ************************************************************************ */
void InterfaceFunction::add(const InterfaceFunction * const eta, double factor)
{
  if(spaceType != eta->getSpaceType())
  {
    ErrThrow("adding interface functions of different space types is not yet ",
             "supported");
  }
  // add factor*eta to this interface function
  // '-1' to take full vector
  this->BlockVector::add(eta->get_entries(), -1, factor); 
}

/** ************************************************************************ */
void InterfaceFunction::zeroAtBoundary()
{
  for(unsigned int i = 0; i < interface.size(); i++)
  {
    if(adjacency[2 * i] == i) // no first neighbor (there is a boundary here)
      this->BlockVector::entries[DOF[3 * i]] = 0;
    if(adjacency[2 * i + 1] == i) //no second neighbor (there is a boundary here)
      this->BlockVector::entries[DOF[3 * i + 2]] = 0;
  }
}

/** ************************************************************************ */
void InterfaceFunction::restrict(const StokesDarcy2D& sd, bool invert)
{
  if(sd.Stokes_first())
    restrict(*(sd.darcy()), invert);
  else
    restrict(*(sd.stokes()), invert);
}

/** ************************************************************************ */
template < class subproblem >
void InterfaceFunction::restrict(const subproblem& sp, bool invert)
{
  if(dofs_to_be_restricted == NULL)
  {
    // create vector 'dofs_to_be_restricted'
    // vector of size zero
    dofs_to_be_restricted = new std::vector<unsigned int>();
    
    // loop over all interface edges, find edges which are touching the boundary
    for(unsigned int iedge = 0; iedge < interface.size(); iedge++)
    {
      if(!isBoundaryEdge(iedge))
        continue; // not an edge touching the boundary
      // this edges touches the boundary, now find interface dof on boundary
      
      // loop over all (2) vertices of this edge. This is only necessary if 
      // this edge touches the boundary with both its vertices (which is rare)
      for(int ivert = 0; ivert < 2; ivert++)
      {
        if(adjacency[2 * iedge + ivert] != iedge)
          continue; // this vertex is not at the boundary, it must be the other
        unsigned int iDOF = DOF[(abs(spaceType) + 1) * iedge
                                + ivert * abs(spaceType)];
        
        // not pretty, but we need this
        double nx, ny, tx, ty;
        getNormal(interface[iedge]->GetNeighb(0), interface[iedge], nx, ny);
        interface[iedge]->GetTangent(tx, ty);
        // this edge is (or is not) directed clockwise in the first neighbor
        int invDir = 0; // invert direction not necessary
        if(tx==ny && ty==-nx)
          invDir = 1; // invert direction necessary
        //else if(tx == -ny && ty == nx)
        
        // function telling us which boundary conditions are prescribed
        BoundCondFunct2D* bc_func = sp.get_example().boundary_conditions[0];
        TBoundEdge* bound_edge = NULL;
        TCollection* coll = sp.get_collection();
        
        // find vertex to which lies on interface and boundary
        TVertex * vert = get_vertex_of_edge(interface[iedge], (ivert+invDir)%2);
        
        // find the boundary edge among all edges sharing this vertex
        std::vector<TJoint *> edges_sharing_vert;
        get_edges_of_vertex(coll, vert, edges_sharing_vert);
        
        // find the edge which is a boundary edge
        for(unsigned int j = 0; j < edges_sharing_vert.size(); j++)
        {
          if(edges_sharing_vert[j]->GetType() == BoundaryEdge)
          {
            bound_edge = (TBoundEdge *) edges_sharing_vert[j];
            break;
          }
        }
        if(bound_edge == NULL)
          // this should not happen
          ErrThrow("could not find boundary edge touching the interface edge, ",
                   "which does touch the boundary!");
        
        double t0, t1;
        bound_edge->GetParameters(t0, t1);
        // boundary component on which this dof is
        int comp = bound_edge->GetBoundComp()->GetID();
        
        // boundary condition at midpoint of this boundary edge
        BoundCond cond;
        bc_func(comp, (t0 + t1) / 2.0, cond);
        //Output::print<1>("boundary edge ", comp, " ", cond, " ", vert);
        if(cond == DIRICHLET)
        {
          dofs_to_be_restricted->push_back(iDOF);
        }
      }
      // sorting is not strictly necessary, but will speed up the restriction
      // if invert==true
      std::sort(dofs_to_be_restricted->begin(), dofs_to_be_restricted->end());
    }
    Output::print<1>("number of dofs on the interface which are restricted: ",
                     dofs_to_be_restricted->size());
  }
  // restrict this interface function to the correct space
  this->restrict(invert);
}

/** ************************************************************************ */
void InterfaceFunction::restrict(bool invert)
{
  if(dofs_to_be_restricted == NULL)
    ErrThrow("you have to call restrict with an object of type StokesProblem ",
             "DarcyProblem or StokesDarcy2D");
  
  if(invert == false)
  {
    for(unsigned int i = 0; i < dofs_to_be_restricted->size(); i++)
    {
      this->BlockVector::entries[(*dofs_to_be_restricted)[i]] = 0;
    }
  }
  else
  {
    // remove all entries but the ones specified by 'dofs_to_be_restricted'
    unsigned int idof = 0;
    for(unsigned int i = 0; i < this->length(); i++)
    {
      if(idof < dofs_to_be_restricted->size()
         && dofs_to_be_restricted->at(idof) == i)
      {
        idof++;
        continue;
      }
      else
        this->BlockVector::entries[i] = 0;
    }
  }
}

/** ************************************************************************ */
double InterfaceFunction::integral() const 
{
  double integral = 0.0;
  if(spaceType != 2 && spaceType != -2)
    ErrThrow("only piecwise quadratic interface functions are supported");
  
  for(unsigned int i = 0; i < interface.size(); i++)
  {
    const double length = interface[i]->GetLength();
    // the three interface functions are, t in [-1,1]   (integral)
    // t*(t-1)/2       ( 1/3 )
    // (1-t*t)         ( 4/3 )
    // t*(t+1)/2       ( 1/3 )
    integral += this->BlockVector::entries[DOF[3 * i]] * length / 6.0;
    integral += this->BlockVector::entries[DOF[3 * i + 1]] * length * 4.0 / 6.0;
    integral += this->BlockVector::entries[DOF[3 * i + 2]] * length / 6.0;
  }
  return integral;
}

/** ************************************************************************ */
void InterfaceFunction::set_integral(double a)
{
  // this function only changes values not corresponding to Dirichlet dofs!
  // we need to know the integral value of this, but neglecting the Dirichlet
  // dofs, so we make a copy and restrict that copy
  InterfaceFunction eta(*this);
  eta = *this; // copy 'this' to eta
  eta.restrict();
  double integral = eta.integral(); // integral of 'this' without Dirichlet dofs
  
  eta = 1; // constant function
  eta.restrict(); // constant function except for Dirichlet dofs
  if(fabs(integral) > 1e-14)
  {
    eta *= integral / eta.integral();
    // now eta has the same integral as 'this' and is constant (except for 
    // Dirichlet dofs)
    *this -= eta;
    // now 'this' has zero integral
  }
  // 'this' has zero integral (neglecting Dirichlet dofs)
  if(a == 0.0)
    return; // only continue if a nonzero integral is desired
  // eta has nonzero integral and is constant (except for Dirichlet dofs)
  eta *= a / eta.integral();
  // eta is now constant with integral 'a'
  *this += eta;
  // now 'this' has integral 'a' (except for Dirichlet dofs)
}

