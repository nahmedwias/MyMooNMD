/**
 * @brief A test program to test some H(div) finite elements
 *
 * If the nodal functionals are \f$N_i\f$ and the basis functions \f$b_i\f$ 
 * then, this function checks that
 * \f[ N_i(b_j) = \delta_{ij} \f]
 * with the Kronecker \f$\delta\f$.
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#ifdef __3D__
#include <FEDatabase3D.h>
#endif // 3D
#include <list>

// =======================================================================
// main program
// =======================================================================

bool CheckFEOnCell(const TFE2D &hdiv_fe, TBaseCell & cell, 
                   TCollection & coll)
{
  // get the points on the reference element, at which we have to evaluate the 
  // basis functions in order to evaluate the nodal functional
  int N_Points;
  const double *xi, *eta;
  hdiv_fe.GetNodalFunctional2D()->GetPointsForAll(N_Points, xi, eta);
  
  // the basis functions object
  TBaseFunct2D* bf = hdiv_fe.GetBaseFunct2D();
  
  // evaluate the basis functions at these points (xi[i],eta[i])
  // dimension of the basis function (usually 1, for H(div) elements it is 2)
  int baseVectDim = bf->GetBaseVectDim();
  // number of basis functions
  int nDof = hdiv_fe.GetN_DOF();
  // number of basis functions, this is the length of the array needed to 
  // evaluate the basis functions (therefore the factor baseVectDim)
  int nBaseFunc = nDof * baseVectDim;
  
  // the id of the reference transformation
  RefTrans2D refTransID = hdiv_fe.GetRefTransID();
  TFEDatabase2D::SetCellForRefTrans(&cell, refTransID);
  
  double uorig[N_Points][2*nDof];
  double AllPointValues[N_Points][nBaseFunc];
  for(int k = 0; k < N_Points; k++)
  {
    bf->GetDerivatives(D00, xi[k], eta[k], AllPointValues[k]);
    // Piola transform
    TFEDatabase2D::GetOrigValues(refTransID, xi[k], eta[k], bf, &coll, 
                                 (TGridCell *) &cell, AllPointValues[k], 
                                 nullptr, nullptr, uorig[k], nullptr, nullptr);
  }
  
  double PointValues[N_Points * baseVectDim];
  double FunctionalValues[nDof];
  bool ret = true;
  for(int k = 0; k < nDof; k++)
  {
    // if this dof is on an edge it must be multiplied with -1 if this cell 
    // index is larger than that of the neighbor at this edge
    int factor = 1;
    // the joint this dof is on (-1 for inner dofs)
    int joint = hdiv_fe.GetFEDesc2D()->GetJointOfThisDOF(k);
    if(joint != -1) // if this dof is not an inner dof
    {
      factor = cell.GetNormalOrientation(joint);
    }
    
    for(int l = 0; l < N_Points; l++)
    {
      for(int i = 0; i < baseVectDim; ++i)
      {
        PointValues[l + i * N_Points] = factor * uorig[l][k + i*nDof];
      }
    }
    
    hdiv_fe.GetNodalFunctional2D()->GetAllFunctionals(&coll, &cell, 
                                                      PointValues, 
                                                      FunctionalValues);
    
    for(int i = 0; i < nDof; i++)
    {
      if( fabs(FunctionalValues[i]) < 1e-10 )
      {
        FunctionalValues[i] = 0;
      }
      //Output::print(k, " ", i, " ", FunctionalValues[i]);
      if( i == k && fabs(FunctionalValues[i]-1) > 1e-8 )
      {
        Output::print("basis function: ", k, " nodal functional: ", i, " ", 
                      FunctionalValues[i]);
        ret = false;
      }
      if( i != k && fabs(FunctionalValues[i]-0) > 1e-8 )
      {
        Output::print("basis function: ", k, " nodal functional: ", i, " ",
                      FunctionalValues[i]);
        ret = false;
      }
    }
    if(ret == false)
      break;
  }
  
  return ret;
}




// check the properties on a given mesh and given elements
bool check(TDomain & domain, const std::list<FE2D>& elements)
{
  unsigned int nRefinements = 3;
  // refine grid up to the coarsest level
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  
  // get the collection of cells on the finest mesh
  TCollection * coll = domain.GetCollection(It_Finest, 0);
  
  // call TFE2D::CheckNFandBF on all these finite elements
  for(auto e : elements)
  {
    Output::print("starting with ", e, " on the reference cell");
    TFE2D * hdiv_fe = TFEDatabase2D::GetFE2D(e);
    // this checks only on the reference cell and does not return true or false
    hdiv_fe->CheckNFandBF();
  }
  
  unsigned int nCells = coll->GetN_Cells();
  // compute global normals on all cells
  for(unsigned int c = 0; c < nCells; ++c)
  {
    coll->GetCell(c)->SetNormalOrientation();
  }
  
  for(auto e : elements)
  {
    Output::print("starting with ", e, " on a mesh");
    TFE2D * hdiv_fe = TFEDatabase2D::GetFE2D(e);
    for(unsigned int c = 0; c < nCells; ++c)
    {
      bool succesful = CheckFEOnCell(*hdiv_fe, *coll->GetCell(c), *coll);
      if(!succesful) // this test failed
      {
        Output::print("test failed with element ", e, " on cell ", c);
        return false;
      }
    }
  }
  delete coll;
  
  return true;
}


#ifdef __3D__
// check the properties on a given mesh and given elements
bool check(TDomain &, const std::list<FE3D>& elements)
{
  // call TFE3D::CheckNFandBF on all these finite elements
  for(auto e : elements)
  {
    Output::print("starting with ", e, " on the reference cell");
    auto * hdiv_fe = TFEDatabase3D::GetFE3D(e);
    // this checks only on the reference cell and does not return true or false
    hdiv_fe->CheckNFandBF();
  }
  
  /// @todo implement the EvalAll routines for the 3D Hdiv elements
  /*
  unsigned int nRefinements = 2;
  // refine grid up to the coarsest level
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  
  // get the collection of cells on the finest mesh
  TCollection * coll = domain.GetCollection(It_Finest, 0);
  
  unsigned int nCells = coll->GetN_Cells();
  // compute global normals on all cells
  for(unsigned int c = 0; c < nCells; ++c)
  {
    coll->GetCell(c)->SetNormalOrientation();
  }
  
  for(auto e : elements)
  {
    Output::print("starting with ", e, " on a mesh");
    auto * hdiv_fe = TFEDatabase3D::GetFE3D(e);
    for(unsigned int c = 0; c < nCells; ++c)
    {
      bool succesful = CheckFEOnCell(*hdiv_fe, *coll->GetCell(c), *coll);
      if(!succesful) // this test failed
      {
        Output::print("test failed with element ", e, " on cell ", c);
        return false;
      }
    }
  }
  delete coll;
  */
  
  return true;
}

#endif // 3D



int main(int, char**)
{
#ifdef __2D__
  //  declaration of databases
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.add("refinement_n_initial_steps", (size_t) 1,"");

  Output::print("************************************************************");
  Output::print("\ntest with quads");
  {
    // the domain is initialised with default description and default
    // initial mesh
    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles"});
    TDomain domain(db);
    
    std::list<FE2D> elements = { N_RT0_2D_Q_M, N_RT1_2D_Q_M, N_RT2_2D_Q_M, 
                                 N_RT3_2D_Q_M, N_BDM1_2D_Q_M, N_BDM2_2D_Q_M,
                                 N_BDM3_2D_Q_M };
    
    if(!check(domain, elements))
      return 1;
  }
  
  Output::print("************************************************************");
  Output::print("test with triangles");
  {
    // the domain is initialised with default description and default
    // initial mesh
    db["geo_file"]= "TwoTriangles";
    TDomain domain(db);
    
    std::list<FE2D> elements = { N_RT0_2D_T_A, N_RT1_2D_T_A, N_RT2_2D_T_A,
                                 N_RT3_2D_T_A, N_BDM1_2D_T_A, N_BDM2_2D_T_A,
                                 N_BDM3_2D_T_A };
    
    if(!check(domain, elements))
      return 1;
  }
  
#else // 2D -> 3D
  //  declaration of databases
  TDatabase Database;
  TFEDatabase3D FEDatabase3d;
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  
  Output::print("************************************************************");
  Output::print("\ntest with hexahedra");
  {
    // the domain is initialised with default description and default
    // initial mesh
    db.add("boundary_file", "Default_UnitCube", "");
    db.add("geo_file", "Default_UnitCube_Hexa", ""); // Default_UnitCube_Tetra
    TDomain domain(db);
    
    std::list<FE3D> elements = { N_RT0_3D_H_A, N_RT1_3D_H_A };
    
    if(!check(domain, elements))
      return 1;
  }
#endif // 3D
  
  return 0;
}
