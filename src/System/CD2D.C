#include <CD2D.h>
#include <LocalAssembling2D.h>
#include <Domain.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <Solver.h>
#include <MultiGrid2D.h>
#include <MainUtilities.h> // L2H1Errors

CD2D::CD2D(TDomain *domain, const Example_CD2D* e)
    : matrix(1, NULL), rhs(1, NULL), function(1, NULL), example(
        e != NULL ? e : new Example_CD2D()), multigrid(NULL)
{
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain->GetCollection(It_Finest, 0);
  
  // create finite elememt space, access through 'function'
  int ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  TFESpace2D* space = new TFESpace2D(coll, (char*) "scalar_space",
                                     (char*) "description",
                                     this->example->get_bc(0), ORDER, NULL);
  int n_dof = space->GetN_DegreesOfFreedom();
  
  // create right hand side and solution
  this->rhs[0] = new double[n_dof];
  double* sol = new double[n_dof]; // access to solution through 'function'
  // set solution and right hand side vectors to zero
  memset(sol, 0, n_dof * SizeOfDouble);
  memset(rhs[0], 0, n_dof * SizeOfDouble);
  
  this->function[0] = new TFEFunction2D(space, (char*) "solution",
                                        (char*) "solution", sol, n_dof);
  
  this->matrix[0] = new TSystemMatScalar2D(space);
  this->matrix[0]->Init(this->example->get_bc(0), this->example->get_bd(0));
  
  // print out some information
  double h_min, h_max;
  coll->GetHminHmax(&h_min, &h_max);
  OutPut("N_Cells    : " << setw(12) << coll->GetN_Cells() << endl);
  OutPut("h (min,max): " << setw(12) << h_min << " " << setw(12) <<h_max<<endl);
  OutPut("dof all    : " << setw(12) << n_dof << endl);
  OutPut("dof active : " << setw(12) << space->GetN_ActiveDegrees() << endl);
  
  // done with the conrtuctor in case we're not using multigrid
  if(TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR != 5 
    || TDatabase::ParamDB->SOLVER_TYPE != 1)
    return;
  // else multigrid
  
  // create spaces, functions, matrices on coarser levels
  double *param = new double[2]; // memory leak
  param[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
  param[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
  this->multigrid = new TMultiGrid2D(1, 2, param);
  // number of refinement levels for the multigrid
  int LEVELS = TDatabase::ParamDB->LEVELS;
  if(LEVELS > domain->get_ref_level() + 1)
    LEVELS = domain->get_ref_level() + 1;
  
  this->function.resize(LEVELS, nullptr);
  this->matrix.resize(LEVELS, nullptr);
  this->rhs.resize(LEVELS, nullptr);
  
  // the matrix and rhs side on the finest grid are already constructed 
  // now construct all matrices, rhs, and solutions on coarser grids
  for(int i = 0; i < LEVELS - 1; i++)
  {
    unsigned int grid = i + domain->get_ref_level() + 1 - LEVELS;
    TCollection *coll = domain->GetCollection(It_EQ, grid);
    // index of the corresponding matrix, rhs, and solution in their respective
    // vectors
    unsigned int index = LEVELS - 1 - i;
    
    space = new TFESpace2D(coll, (char*) "p", (char*) "p", example->get_bc(0),
                           ORDER, NULL);
    n_dof = space->GetN_DegreesOfFreedom();
    this->matrix.at(index) = new TSystemMatScalar2D(space);
    this->matrix.at(index)->Init(this->example->get_bc(0),
                                 this->example->get_bd(0));
    
    this->rhs[index] = new double[n_dof];
    sol = new double[n_dof]; // access to solution through 'function'
    // set solution and right hand side vectors to zero
    memset(sol, 0, n_dof * SizeOfDouble);
    memset(this->rhs[index], 0, n_dof * SizeOfDouble);
    
    this->function[index] = new TFEFunction2D(space, (char*) "solution",
                                              (char*) "solution", sol, n_dof);
    
    TMGLevel2D *multigrid_level = new TMGLevel2D(
        i, this->matrix[index]->get_square_matrix(), rhs[index], sol, 2, NULL);
    this->multigrid->AddLevel(multigrid_level);
  }
  // add last multigrid level (on finest mesh)
  TMGLevel2D *multigrid_level = new TMGLevel2D(
      LEVELS - 1, this->matrix[0]->get_square_matrix(), rhs[0], 
      this->function[0]->GetValues(), 2, NULL);
  this->multigrid->AddLevel(multigrid_level);
}

CD2D::~CD2D()
{
  // delete matrix
  for(auto mat : this->matrix)
    delete mat;
  for(auto r : this->rhs)
    delete[] r;
  for(auto f : this->function)
  {
    delete f->GetFESpace2D();
    delete f;
  }
  delete multigrid;
}

void CD2D::assemble()
{
  LocalAssembling2D_type t = CD2D_Galerkin;
  switch(TDatabase::ParamDB->DISCTYPE)
  {
    case 1:
      t = CD2D_Galerkin;
      break;
    case 2:
      t = CD2D_SUPG;
      break;
    case 6:
      t = CD2D_GLS;
      break;
    default:
      ErrMsg(
          "currently DISCTYPE " << TDatabase::ParamDB->DISCTYPE << " is not " << "supported by the class CD2D")
      ;
      throw("unsupported DISCTYPE");
  }
  
  // this loop has more than one iteration only in case of multigrid
  for(unsigned int grid = 0, n_grids = matrix.size(); grid < n_grids; ++grid)
  {
    // create a local assembling object which is needed to assemble the matrix
    LocalAssembling2D la(t, &(this->function[grid]),
                         this->example->get_coeffs());
    // assemble the system matrix with given local assembling, solution and rhs 
    this->matrix[grid]->Assemble(la, this->function[grid]->GetValues(),
                                 this->rhs[grid]);
  }
}

void CD2D::solve()
{
  double t = GetTime();
  
  TSquareMatrix2D *SqMat[1] = { this->getMatrix()->get_square_matrix() };
  Solver((TSquareMatrix **)SqMat, NULL, this->rhs[0], 
         this->function[0]->GetValues(), MatVect_Scalar, Defect_Scalar, 
         this->multigrid, this->getSize(), 0);
  
  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    t = GetTime() - t;
    OutPut(" solving of a CD2D problem done in " << t << " seconds\n");
  }
}

void CD2D::output(int i)
{
  if(!TDatabase::ParamDB->WRITE_VTK && !TDatabase::ParamDB->MEASURE_ERRORS)
    return;
  
  // print the value of the largest and smallest entry in the finite element 
  // vector
  this->function[0]->PrintMinMax();
  
  // write solution to a vtk file
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    // last argument in the following is domain, but is never used in this class
    TOutput2D Output(1, 1, 0, 0, NULL);
    Output.AddFEFunction(this->function[0]);
    std::string filename(TDatabase::ParamDB->OUTPUTDIR);
    filename += "/" + std::string(TDatabase::ParamDB->BASENAME);
    if(i >= 0)
      filename += "_" + std::to_string(i);
    filename += ".vtk";
    Output.WriteVtk(filename.c_str());
  }
  
  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(TDatabase::ParamDB->MEASURE_ERRORS)
  {
    double errors[4];
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
    TFESpace2D* space = this->function[0]->GetFESpace2D();
    
    this->function[0]->GetErrors(this->example->get_exact(0), 3, AllDerivatives,
                                 4, SDFEMErrors, this->example->get_coeffs(),
                                 &aux, 1, &space, errors);
    
    OutPut("L2     : " << errors[0] << endl);
    OutPut("H1-semi: " << errors[1] << endl);
    OutPut("SD     : " << errors[2] << endl);
    OutPut("L_inf  : " << errors[3] << endl);
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)
}
