#include <CD3D.h>
//#include <LocalAssembling3D.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <Output3D.h>
#include <LinAlg.h>

#include <Solver.h>
#include <DirectSolver.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <JacobiIte.h>
#include <MultiGridScaIte.h>

#include <MultiGrid3D.h>
#include <MainUtilities.h> // L2H1Errors

CD3D::CD3D(TDomain *domain, const Example_CD3D* e)
    : matrix(1, NULL), rhs(1, NULL), function(1, NULL), example(
        e != NULL ? e : new Example_CD3D()), multigrid(NULL)
{
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain->GetCollection(It_Finest, 0);
  
  // create finite elememt space, access through 'function'
  int ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  TFESpace3D* space = new TFESpace3D(coll, (char*) "scalar_space",
                                     (char*) "description",
                                     this->example->get_bc(0), ORDER);
  int n_dof = space->GetN_DegreesOfFreedom();
  
  // create right hand side and solution
  this->rhs[0] = new double[n_dof];
  double* sol = new double[n_dof]; // access to solution through 'function'
  // set solution and right hand side vectors to zero
  memset(sol, 0, n_dof * SizeOfDouble);
  memset(rhs[0], 0, n_dof * SizeOfDouble);
  
  this->function[0] = new TFEFunction3D(space, (char*) "solution",
                                        (char*) "solution", sol, n_dof);
  
  this->matrix[0] = new BlockMatrixCD3D(space);
  this->matrix[0]->Init(this->example->get_coeffs(), this->example->get_bc(0),
                        this->example->get_bd(0));
  
  // print out some information
  double h_min, h_max;
  coll->GetHminHmax(&h_min, &h_max);
  OutPut("N_Cells    : " << setw(12) << coll->GetN_Cells() << endl);
  OutPut("h (min,max): " << setw(12) << h_min << " " << setw(12) <<h_max<<endl);
  OutPut("dof all    : " << setw(12) << n_dof << endl);
  OutPut("dof active : " << setw(12) << space->GetN_ActiveDegrees() << endl);
  
  // done with the conrtuctor in case we're not using multigrid
  if(TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5
    && TDatabase::ParamDB->SOLVER_TYPE == 1)
  {
    // multigrid
    
    // create spaces, functions, matrices on coarser levels
    double *param = new double[2]; // memory leak
    param[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
    param[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
    this->multigrid = new TMultiGrid3D(1, 2, param);
    // number of refinement levels for the multigrid
    int LEVELS = TDatabase::ParamDB->LEVELS;
    if(LEVELS > domain->get_ref_level() + 1)
      LEVELS = domain->get_ref_level() + 1;
    
    this->function.resize(LEVELS, nullptr);
    this->matrix.resize(LEVELS, nullptr);
    this->rhs.resize(LEVELS, nullptr);
    
    // some parameter used to construct the multigrid object
    int n_aux = 2;
    if(TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR
       || TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR)
      n_aux = 4;
    
    // the matrix and rhs side on the finest grid are already constructed 
    // now construct all matrices, rhs, and solutions on coarser grids
    for(int i = 0; i < LEVELS - 1; i++)
    {
      unsigned int grid = i + domain->get_ref_level() + 1 - LEVELS;
      TCollection *coll = domain->GetCollection(It_EQ, grid);
      // index of the corresponding matrix, rhs, and solution in their
      // respective vectors
      unsigned int index = LEVELS - 1 - i;
      
      space = new TFESpace3D(coll, (char*) "p", (char*) "p", example->get_bc(0),
                             ORDER);
      n_dof = space->GetN_DegreesOfFreedom();
      this->matrix.at(index) = new BlockMatrixCD3D(space);
      this->matrix.at(index)->Init(this->example->get_coeffs(),
                                   this->example->get_bc(0),
                                   this->example->get_bd(0));
      
      this->rhs[index] = new double[n_dof];
      sol = new double[n_dof]; // access to solution through 'function'
      // set solution and right hand side vectors to zero
      memset(sol, 0, n_dof * SizeOfDouble);
      memset(this->rhs[index], 0, n_dof * SizeOfDouble);
      
      this->function[index] = new TFEFunction3D(space, (char*) "solution",
                                                (char*) "solution", sol, n_dof);
      
      TMGLevel3D *multigrid_level = new TMGLevel3D(
          i, this->matrix[index]->get_square_matrix(), rhs[index], sol, n_aux,
          NULL);
      this->multigrid->AddLevel(multigrid_level);
    }
    // add last multigrid level (on finest mesh)
    TMGLevel3D *multigrid_level = new TMGLevel3D(
        LEVELS - 1, this->matrix[0]->get_square_matrix(), rhs[0], 
        this->function[0]->GetValues(), n_aux, NULL);
    this->multigrid->AddLevel(multigrid_level);
  }
  
  #ifdef _MPI
  double t1,t2,tdiff;
  Comm = MPI_COMM_WORLD;
  unsigned int n_levels = this->matrix.size();
  
  ParMapper = new TParFEMapper3D*[n_levels]; 
  ParComm   = new TParFECommunicator3D*[n_levels];
  
  
  int MaxSubDomainPerDof = 4; //! todo this is known in main only
  
  int out_rank = TDatabase::ParamDB->Par_P0;
  int rank;
  MPI_Comm_rank(Comm, &rank);
  
  if(TDatabase::ParamDB->timeprofiling)  t1 = MPI_Wtime();
  for(unsigned int i = 0; i < n_levels; i++)
  {
    this->function[i]->GetFESpace3D()->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
    
    ParMapper[i] = new TParFEMapper3D(
        1, this->function[i]->GetFESpace3D(),
        this->matrix[i]->get_square_matrix()->GetRowPtr(),
        this->matrix[i]->get_square_matrix()->GetKCol());
    ParComm[i]   = new TParFECommunicator3D(ParMapper[i]);
  }
  
  if(TDatabase::ParamDB->SOLVER_TYPE == DIRECT)
  {
    TSquareMatrix3D* a = this->matrix[0]->get_square_matrix();
    DS = new TParDirectSolver(ParComm[0], NULL,
                              &a, NULL);
  }
  
  if(TDatabase::ParamDB->timeprofiling)
  {
    t2 = MPI_Wtime();
    tdiff = t2-t1;
    int out_rank = TDatabase::ParamDB->Par_P0;
    int rank;
    MPI_Comm_rank(Comm, &rank);
    MPI_Reduce(&tdiff, &t1, 1, MPI_DOUBLE, MPI_MIN, out_rank, Comm);
    if(rank == out_rank)
    {
      printf("Time taken for FeSpace SubDomain dof mapping for all levels is %e\n", t1);
    }
  }
  #endif

  #ifdef _OMPONLY
  if(TDatabase::ParamDB->SOLVER_TYPE == DIRECT 
      && TDatabase::ParamDB->DSType == 1)
  {
    DS = new TParDirectSolver(this->matrix[0]->get_square_matrix());
  }
  #endif 
}

CD3D::~CD3D()
{
  // delete matrix
  for(auto mat : this->matrix)
    delete mat;
  for(auto r : this->rhs)
    delete[] r;
  for(auto f : this->function)
    delete f;
  delete multigrid;
}

void CD3D::assemble()
{
  //LocalAssembling3D_type t = CD3D;
  
  // this loop has more than one iteration only in case of multigrid
  for(unsigned int grid = 0, n_grids = matrix.size(); grid < n_grids; ++grid)
  {
    // create a local assembling object which is needed to assemble the matrix
    //LocalAssembling3D la(t, &(this->function[grid]),
    //                     this->example->get_coeffs());
    // assemble the system matrix with given local assembling, solution and rhs 
    //this->matrix[grid]->Assemble(la, this->function[grid]->GetValues(),
    //                             this->rhs[grid]);
    this->matrix[grid]->Assemble(this->example->get_coeffs(),
                                 this->function[grid]->GetValues(), rhs[grid]);
  }
  

  //have to shift this in pardirectsolver    
  #ifdef _OMPONLY     
  if(TDatabase::ParamDB->SOLVER_TYPE == DIRECT && TDatabase::ParamDB->DSType == 1)
    DS->AssembleMatrix(this->matrix->get_square_matrix());
  #endif
}

void CD3D::solve()
{
  double t = GetTime();
  
  TSquareMatrix3D *SqMat[1] = { this->getMatrix()->get_square_matrix() };
  
  int N_Levels = TDatabase::ParamDB->LEVELS;
  if(TDatabase::ParamDB->SC_MG_TYPE_SCALAR)
    ++N_Levels;
  
  switch(TDatabase::ParamDB->SOLVER_TYPE)
  {
    case AMG_SOLVE:
      Solver((TSquareMatrix *)this->getMatrix()->get_square_matrix(), rhs[0],
             this->function[0]->GetValues());
      break;
    case GMG:
    {
      int n_dof = this->function[0]->GetFESpace3D()->GetN_DegreesOfFreedom();
      double *Itmethod_sol;
      double *Itmethod_rhs;
      
      TItMethod *iterative_method;
      TItMethod *preconditioner;
      
      // build preconditioner
      switch(TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR)
      {
        case 1:
          preconditioner = new TJacobiIte(MatVect_Scalar, Defect_Scalar, NULL,
                                          0, n_dof, 1
                                          #ifdef _MPI
                                          ,ParComm[0]
                                          #endif
                                          );
          
          Itmethod_sol = this->function[0]->GetValues();
          Itmethod_rhs = rhs[0];
          break;
        case 5:
          preconditioner = new TMultiGridScaIte(MatVect_Scalar, Defect_Scalar,
                                                NULL, 0, n_dof, this->multigrid,
                                                1);
          Itmethod_sol = new double[n_dof];
          Itmethod_rhs = new double[n_dof];
          memcpy(Itmethod_sol, this->function[0]->GetValues(),
                 n_dof * SizeOfDouble);
          memcpy(Itmethod_rhs, rhs[0], n_dof * SizeOfDouble);
          break;
        default:
          ErrMsg("Unknown preconditioner !!!");
          throw("Unknown preconditioner !!!");
      }
      switch(TDatabase::ParamDB->SC_SOLVER_SCALAR)
      {
        case 11:
          iterative_method = new TFixedPointIte(MatVect_Scalar, Defect_Scalar,
                                                preconditioner, 0, n_dof, 1
                                                #ifdef _MPI
                                                , ParComm[0]
                                                #endif
                                                );
          break;
        case 16:
          iterative_method = new TFgmresIte(MatVect_Scalar, Defect_Scalar,
                                            preconditioner, 0, n_dof, 1
                                            #ifdef _MPI
                                            , ParComm[0]
                                            #endif
                                            );
          break;
        default:
          ErrMsg("Unknown solver !!!");
          throw("unknown solver");
      }
      
      // solve linear system
      iterative_method->Iterate((TSquareMatrix **)SqMat, NULL, Itmethod_sol,
                                Itmethod_rhs);
      
      if(TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
      {
        memcpy(this->function[0]->GetValues(), Itmethod_sol,
               n_dof * SizeOfDouble);
        memcpy(rhs[0], Itmethod_rhs, n_dof * SizeOfDouble);
        delete [] Itmethod_sol;
        delete [] Itmethod_rhs;
      }
      delete iterative_method;
      delete preconditioner;
      break;
    }
    case DIRECT:
      #ifdef _MPI
      DS->Solve(this->function[0]->GetValues(), rhs[0], true);
      #endif
      
      #ifdef _OMPONLY
      if(TDatabase::ParamDB->DSType == 1)
        DS->Solve(this->function[0]->GetValues(), rhs[0], true);
      else
      { 
        ErrMsg("Select Proper Solver");
        exit(0);
      }
      #endif
      
      #ifdef _SEQ
      DirectSolver((TSquareMatrix*)this->getMatrix()->get_square_matrix(),
                   rhs[0], this->function[0]->GetValues());
      #endif
      break;
    default:
      OutPut("Unknown Solver" << endl);
      exit(4711);
  }
  
  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    t = GetTime() - t;
    OutPut(" solving of a CD3D problem done in " << t << " seconds\n");
  }
}

void CD3D::output(int i)
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
    TOutput3D Output(1, 1, 0, 0, NULL);
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
    TAuxParam3D aux(1, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
    MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };
    TFESpace3D* space = this->function[0]->GetFESpace3D();
    
    this->function[0]->GetErrors(this->example->get_exact(0), 4, AllDerivatives,
                                 2, L2H1Errors, this->example->get_coeffs(),
                                 &aux, 1, &space, errors);
    
    OutPut("L2     : " << errors[0] << endl);
    OutPut("H1-semi: " << errors[1] << endl);
    //OutPut("SD     : " << errors[2] << endl);
    //OutPut("L_inf  : " << errors[3] << endl);
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)
}

#ifdef _MPI
unsigned int CD3D::get_total_dof() const
{
  int n_dof = this->getSize();
  int n_total_dof = 0;
  int out_rank = TDatabase::ParamDB->Par_P0;
  MPI_Reduce(&n_dof, &n_total_dof, 1, MPI_INT, MPI_SUM, out_rank, Comm);
  return (unsigned int) n_total_dof;
}

unsigned int CD3D::get_n_total_cells() const
{
  int n_cells = this->getSpace()->GetN_Cells();
  int n_total_cells = 0;
  int out_rank = TDatabase::ParamDB->Par_P0;
  MPI_Reduce(&n_cells, &n_total_cells, 1, MPI_INT, MPI_SUM, out_rank, Comm);
  return (unsigned int) n_total_cells;
}
#endif
