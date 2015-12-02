#include <CD2D.h>
#include <Database.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <Solver.h>
#include <MultiGrid2D.h>
#include <MainUtilities.h> // L2H1Errors
#include <AlgebraicFluxCorrection.h>


/** ************************************************************************ */
CD2D::System_per_grid::System_per_grid(const Example_CD2D& example,
                                       TCollection& coll)
 : fe_space(&coll, (char*)"space", (char*)"cd2d fe_space", example.get_bc(0),
            TDatabase::ParamDB->ANSATZ_ORDER, nullptr),
   matrix(this->fe_space, example.get_bd(0)),
   rhs(this->matrix, true),
   solution(this->matrix, false),
   fe_function(&this->fe_space, (char*)"c", (char*)"c",
               this->solution.get_entries(), this->solution.length())
{
  
}

/** ************************************************************************ */
CD2D::CD2D(const TDomain& domain, int reference_id)
 : CD2D(domain, Example_CD2D(), reference_id)
{
}

/** ************************************************************************ */
CD2D::CD2D(const TDomain& domain, Example_CD2D example, int reference_id)
 : systems(), example(example), multigrid(nullptr)
{
  this->set_parameters();
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
  
  // when using afc, create system matrices as if all dofs were active
  if(TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION > 0)
  {
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1;
  }

  // create finite element space and function, a matrix, rhs, and solution
  this->systems.emplace_back(this->example, *coll);
  
  
  // print out some information
  TFESpace2D & space = this->systems.front().fe_space;
  double h_min, h_max;
  coll->GetHminHmax(&h_min, &h_max);
  Output::print<1>("N_Cells    : ", setw(12), coll->GetN_Cells());
  Output::print<2>("h (min,max): ", setw(12), h_min, " ", setw(12), h_max);
  Output::print<1>("dof all    : ", setw(12), space.GetN_DegreesOfFreedom());
  Output::print<2>("dof active : ", setw(12), space.GetN_ActiveDegrees());
  
  
  // done with the conrtuctor in case we're not using multigrid
  if(TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR != 5 
    || TDatabase::ParamDB->SOLVER_TYPE != 1)
    return;
  // else multigrid
  
  // create spaces, functions, matrices on coarser levels
  double *param = new double[2]; // memory leak
  param[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
  param[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
  this->multigrid.reset(new TMultiGrid2D(1, 2, param));
  // number of refinement levels for the multigrid
  int LEVELS = TDatabase::ParamDB->LEVELS;
  if(LEVELS > domain.get_ref_level() + 1)
    LEVELS = domain.get_ref_level() + 1;
  
  // the matrix and rhs side on the finest grid are already constructed 
  // now construct all matrices, rhs, and solutions on coarser grids
  for(int i = LEVELS - 2; i >= 0; i--)
  {
    unsigned int grid = i + domain.get_ref_level() + 1 - LEVELS;
    TCollection *coll = domain.GetCollection(It_EQ, grid, reference_id);
    this->systems.emplace_back(example, *coll);
  }
  
  // create multigrid-level-objects, must be coarsest first
  unsigned int i = 0;
  for(auto it = this->systems.rbegin(); it != this->systems.rend(); ++it)
  {
    TMGLevel2D *multigrid_level = new TMGLevel2D(
      i, it->matrix.get_matrix(), it->rhs.get_entries(), 
      it->solution.get_entries(), 2, NULL);
    i++;
    this->multigrid->AddLevel(multigrid_level);
  }
}

/** ************************************************************************ */
CD2D::~CD2D()
{
  // delete the collections created during the contructor
  for(auto & s : this->systems)
    delete s.fe_space.GetCollection();
}

/** ************************************************************************ */
void CD2D::set_parameters()
{
  //////////////// Algebraic flux correction ////////////
  if(TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION != 0)
  {//some kind of afc enabled
    //make sure that galerkin discretization is used
    if (TDatabase::ParamDB->DISCTYPE !=	1)
    {//some other disctype than galerkin
      TDatabase::ParamDB->DISCTYPE = 1;
      Output::print("DISCTYPE changed to 1 (GALERKIN) because Algebraic Flux ",
                    "Correction is enabled.");
    }
  }
}

/** ************************************************************************ */
void CD2D::assemble()
{
  LocalAssembling2D_type t = LocalAssembling2D_type::ConvDiff;

  // this loop has more than one iteration only in case of multigrid
  for(auto & s : this->systems)
  {
    TFEFunction2D * pointer_to_function = &s.fe_function;
    // create a local assembling object which is needed to assemble the matrix
    LocalAssembling2D la(t, &pointer_to_function, this->example.get_coeffs());
    // assemble the system matrix with given local assembling, solution and rhs 
    s.matrix.Assemble(la, s.solution, s.rhs);
  }

  // when using afc, do it now
  if(TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION > 0)
  {
    performAlgebraicFluxCorrection();
  }

}

/** ************************************************************************ */
void CD2D::solve()
{
  double t = GetTime();
  System_per_grid& s = this->systems.front();
  TSquareMatrix2D *SqMat[1] = { s.matrix.get_matrix() };
  Solver((TSquareMatrix **)SqMat, NULL, s.rhs.get_entries(), 
         s.solution.get_entries(), MatVect_Scalar, Defect_Scalar, 
         this->multigrid.get(), this->get_size(), 0);
  
  t = GetTime() - t;
  Output::print<2>(" solving of a CD2D problem done in ", t, " seconds");
}

/** ************************************************************************ */
void CD2D::output(int i)
{
  if(!TDatabase::ParamDB->WRITE_VTK && !TDatabase::ParamDB->MEASURE_ERRORS)
    return;
  
  // print the value of the largest and smallest entry in the finite element 
  // vector
  TFEFunction2D & fe_function = this->systems.front().fe_function;
  fe_function.PrintMinMax();
  
  // write solution to a vtk file
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    // last argument in the following is domain, but is never used in this class
    TOutput2D Output(1, 1, 0, 0, NULL);
    Output.AddFEFunction(&fe_function);
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
    double errors[5];
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
    const TFESpace2D* space = fe_function.GetFESpace2D();
    
    fe_function.GetErrors(this->example.get_exact(0), 3, AllDerivatives, 4,
                          SDFEMErrors, this->example.get_coeffs(), &aux, 1, 
                          &space, errors);
    
    Output::print<1>("L2     : ", errors[0]);
    Output::print<1>("H1-semi: ", errors[1]);
    Output::print<1>("SD     : ", errors[2]);
    Output::print<1>("L_inf  : ", errors[3]);
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)
}

/** ************************************************************************ */
void CD2D::performAlgebraicFluxCorrection()
{
  for(auto & s : this->systems) // do it on all levels
  {
    //determine which kind of afc to use
    switch (TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION)
    {
      case 1: //FEM-TVD
      {
        //get pointers to the relevant objects
        TFESpace2D& feSpace = s.fe_space;
        TSquareMatrix2D * matrix = s.matrix.get_matrix();
        double* solEntries = s.solution.get_entries();
        double* rhsEntries = s.rhs.get_entries();

        // fill an array "neumannToDirichlet" with those rows, that got internally treated as
        // Neumann although they are Dirichlet
        size_t nNeumannToDirichlet = feSpace.GetN_Dirichlet();
        int* neumannToDirichlet = new int[nNeumannToDirichlet];
        int dirichletDofStartIndex = feSpace.GetDirichletBound();
        int* dirichletDofStartPtr = &feSpace.GetGlobalNumbers()[dirichletDofStartIndex];
        for (size_t i = 0; i < nNeumannToDirichlet ;++i){
          neumannToDirichlet[i]= dirichletDofStartPtr[i];
        }

        // Number of dofs.
        int nDofs = feSpace.GetN_DegreesOfFreedom();

        // array of entries for matrix D
        double* entriesMatrixD = new double[matrix->GetN_Entries()]();

        // apply FEM-TVD
        AlgebraicFluxCorrection::FEM_TVD_ForConvDiff(
            matrix, nDofs, nDofs,
            entriesMatrixD,
            solEntries,rhsEntries,
            nNeumannToDirichlet, neumannToDirichlet, 1);

        //...and finally correct the entries in the Dirchlet rows
        AlgebraicFluxCorrection::correctDirichletRows(*matrix);

        //clean up
        delete[] entriesMatrixD;
        delete[] neumannToDirichlet;
        break;
      }
      default:
      {
        ErrThrow("The chosen ALGEBRAIC_FLUX_CORRECTION scheme is unknown to class CD2D.");
      }
    }
  }
}
