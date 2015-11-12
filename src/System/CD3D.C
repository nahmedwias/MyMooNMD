#include <CD3D.h>
//#include <LocalAssembling3D.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <Output3D.h>
#include <LinAlg.h>

//#include <Solver.h>
#include <DirectSolver.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <JacobiIte.h>
#include <MultiGridScaIte.h>

#include <MultiGrid3D.h>
#include <MainUtilities.h> // L2H1Errors


#ifdef _MPI
	CD3D::SystemPerGrid::SystemPerGrid(const Example_CD3D& example,
	                                       TCollection& coll, int maxSubDomainPerDof)
	 : feSpace_(&coll, (char*)"space", (char*)"cd3d fe_space", example.get_bc(0),
	            TDatabase::ParamDB->ANSATZ_ORDER),
	   matrix_(feSpace_, example.get_bd(0)), //system block matrix
	   rhs_(matrix_, true), // suitable right hand side vector filled with zeroes
	   solution_(matrix_, false), // suitable solution vector filled with zeroes
	   feFunction_(&feSpace_, (char*)"c", (char*)"c", solution_.get_entries(), solution_.length()),
	   parMapper_(), // will be reset shortly
	   parComm_() // will be reset shortly

	{
	  //inform the fe space about the maximum number of subdomains per dof
		feSpace_.SetMaxSubDomainPerDof(maxSubDomainPerDof);

		//Must be reset here, because feSpace needs special treatment
		//This includes copy assignment - all because there is no good
		// way to communicate Maximum number of subdomains per dof to FESpace...
		parMapper_ = TParFEMapper3D(3, &feSpace_, matrix_.get_matrix()->GetRowPtr(), //magic number 3 is space dimension
				   	  matrix_.get_matrix()->GetKCol());
		parComm_ = TParFECommunicator3D(&parMapper_);
	}
#else
	/** ************************************************************************ */
	CD3D::SystemPerGrid::SystemPerGrid(const Example_CD3D& example,
	                                       TCollection& coll)
	 : feSpace_(&coll, (char*)"space", (char*)"cd3d fe_space", example.get_bc(0),
	            TDatabase::ParamDB->ANSATZ_ORDER),
	   matrix_(feSpace_, example.get_bd(0)), //system block matrix
	   rhs_(matrix_, true), // suitable right hand side vector filled with zeroes
	   solution_(matrix_, false), // suitable solution vector filled with zeroes
	   feFunction_(&feSpace_, (char*)"c", (char*)"c", solution_.get_entries(), solution_.length())
	{

	}
#endif

/** ************************************************************************ */
CD3D::CD3D(const TDomain& domain, const Example_CD3D& example
#ifdef _MPI
,int maxSubDomainPerDof
#endif
)
    : systems_(), example_(example), multigrid_(nullptr)
{

  // Create the collection of cells from the finest grid of the domain
#ifdef _MPI
	int mpiRank;
	int mpiSize;
	MPI_Comm& comm = TDatabase::ParamDB->Comm ;
	MPI_Comm_rank(comm, &mpiRank); 	//find out which rank I am
	MPI_Comm_rank(comm, &mpiSize);	//find out how many processes are in the game

	// create collection of mesh cells
	TCollection* cellCollection;
//	if(TDatabase::ParamDB->MapperType == 2)
//	{
//		//CB I am not entirely sure, but I think this is the way to go for old master/slave mapping type
//		// - exclude halo cells from the collection.
//		cellCollection = domain.GetOwnCollection(It_Finest, 0, mpiRank);
//	}
//	else
//	{
		//collection consisting of own and halo cells
		cellCollection = domain.GetCollection(It_Finest, 0);
//	}

	// create finite element space and function, a matrix, rhs, and solution
	systems_.emplace_back(example_, *cellCollection, maxSubDomainPerDof);

#else
	TCollection *cellCollection = domain.GetCollection(It_Finest, 0);
	// create finite element space and function, a matrix, rhs, and solution
	systems_.emplace_back(example_, *cellCollection);

    // print out some information
    TFESpace3D & space = this->systems_.front().feSpace_;
    double hMin, hMax;
    cellCollection->GetHminHmax(&hMin, &hMax);
    OutPut("N_Cells    : " << setw(12) << cellCollection->GetN_Cells() << endl);
    OutPut("h (min,max): " << setw(12) << hMin << " " << setw(12) <<hMax<<endl);
    OutPut("dof all    : " << setw(12) << space.GetN_DegreesOfFreedom() << endl);
    OutPut("dof active : " << setw(12) << space.GetN_ActiveDegrees() << endl);
#endif



  // done with the constructor in case we're not using multigrid
  if(TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR != 5
    || TDatabase::ParamDB->SOLVER_TYPE != 1)
  {
      return;
  }
  else
  {
	  ErrMsg("Multigrid not yet implemented, do not use SC_PRECONDITIONER_SCALAR: 5 or SOLVER_TYPE: 1 " << endl);
	  exit(-1);
	  // else multigrid TODO Enable this class to use iterative solver with mg prec and pure multigrid solve.
  }


}

/** ************************************************************************ */
void CD3D::assemble()
{
  LocalAssembling3D_type laType = LocalAssembling3D_type::CD3D;
  // this loop has more than one iteration only in case of multigrid
  for(auto & s : systems_)
  {
    TFEFunction3D * feFunctionPtr = &s.feFunction_;
    // create a local assembling object which is needed to assemble the matrix
    LocalAssembling3D laObject(laType, &feFunctionPtr, example_.get_coeffs());
    // assemble the system matrix with given local assembling, solution and rhs
    s.matrix_.assemble(laObject, s.solution_, s.rhs_);
  }

}

/** ************************************************************************ */
void CD3D::solve()
{

  //hold a reference to last grid system (should be the finest)
  SystemPerGrid& syst = systems_.front() ;

  //CB Try if we can get a preconditioned iterative method working.

  //Hold some variable which will be needed in building preconditioner and solver.
  int nDof = syst.feSpace_.GetN_DegreesOfFreedom();
  TItMethod* preconditioner;
  TItMethod* iterativeSolver;

  TSquareMatrix3D *sqMat[1] = { syst.matrix_.get_matrix() };
  double* solutionEntries = syst.solution_.get_entries();
  double* rhsEntries =syst.rhs_.get_entries();

  #ifdef _MPI
  TParFECommunicator3D* parComm = &syst.parComm_;
  #endif _MPI

//  if (TDatabase::ParamDB->SOLVER_TYPE == 2)
//  { //Direct solver(s)
//	  TSquareMatrix3D *SqMat = syst.matrix_.get_matrix();
//	  DirectSolver(SqMat, syst.rhs_.get_entries(), syst.solution_.get_entries());
//  }
//  else if (TDatabase::ParamDB->SOLVER_TYPE == 1)
//  {   // iterative solver(s)

	  // Determine and build preconditioner.
	  switch (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR)
	  {
	  case 1: //Jacobi Iteration
	  {
		  #ifdef _MPI
		  preconditioner = new TJacobiIte(MatVect_Scalar, Defect_Scalar, NULL, 0, nDof, 1, parComm);
		  #else
		  preconditioner = new TJacobiIte(MatVect_Scalar, Defect_Scalar, NULL, 0, nDof, 1);
		  #endif
		  break;
	  }
	  case 5: //TMultiGridScaIte
		  ErrMsg("SC_PRECONDITIONER_SCALAR: " << TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR << " has to be implemented." << endl);
		  break;
	  default:
		  ErrMsg("Unknown SC_PRECONDITIONER_SCALAR: " << TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR);
	  } //end building preconditioner

	  switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
	  {
	  case 11: //fixed point iteration TODO What is the actual name of the method?
	  {
		  #ifdef _MPI
	    iterativeSolver = new TFixedPointIte(MatVect_Scalar, Defect_Scalar, preconditioner, 0, nDof, 1, parComm);
		  #else
	    iterativeSolver = new TFixedPointIte(MatVect_Scalar, Defect_Scalar, preconditioner, 0, nDof, 1);
		  #endif
		  break;
	  }
	  case 16:
		  ErrMsg("SC_SOLVER_SCALAR: " << TDatabase::ParamDB->SC_SOLVER_SCALAR << " (FGMRES) has to be implemented." << endl);
		  break;
	  default:
		  ErrMsg("Unknown solver !!!" << endl);
	  }

	  //call the solver
	  iterativeSolver->Iterate((TSquareMatrix**) sqMat, nullptr, solutionEntries, rhsEntries);

	  delete iterativeSolver;
//  }
//  else
//  {
//	  ErrThrow("Unknown SOLVER_TYPE. Choose either 1 (iterative) or 2 (direct).");
//  }
}

void CD3D::output(int i)
{
	if(!TDatabase::ParamDB->WRITE_VTK && !TDatabase::ParamDB->MEASURE_ERRORS)
	{
		return;
	}

	//hold a reference to the finest grid system
    SystemPerGrid& syst = systems_.back() ;

	// print the value of the largest and smallest entry in the FE vector
    syst.feFunction_.PrintMinMax();

	// write solution to a vtk file
	if(TDatabase::ParamDB->WRITE_VTK)
	{
		// last argument in the following is domain, but is never used in this class
		TOutput3D Output(1, 1, 0, 0, NULL);
		Output.AddFEFunction(&syst.feFunction_);
#ifdef _MPI
		char SubID[] = "";
		Output.Write_ParVTK(MPI_COMM_WORLD, 0, SubID);
#else
		std::string filename(TDatabase::ParamDB->OUTPUTDIR);
		filename += "/" + std::string(TDatabase::ParamDB->BASENAME);
		if(i >= 0)
			filename += "_" + std::to_string(i);
		filename += ".vtk";
		Output.WriteVtk(filename.c_str());
#endif
	}

	// measure errors to known solution
	// If an exact solution is not known, it is usually set to be zero, so that
	// in such a case here only integrals of the solution are computed.
	if(TDatabase::ParamDB->MEASURE_ERRORS)
	{
		double errors[4];
		TAuxParam3D aux(1, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
		MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };
		TFESpace3D* space = syst.feFunction_.GetFESpace3D();

		syst.feFunction_.GetErrors(example_.get_exact(0), 4, AllDerivatives,
								   2, L2H1Errors, example_.get_coeffs(),
								   &aux, 1, &space, errors);
		#ifdef _MPI
		// usual code block to gather information about this
		// processes role in the mpi communicator
		MPI_Comm globalComm = TDatabase::ParamDB->Comm;
		int mpiRank, mpiSize;
		MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
		MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
		bool iAmOutRank= (mpiRank == TDatabase::ParamDB->Par_P0);

		double errorsReduced[4]; //memory for global (across all processes) error

		MPI_Allreduce(errors, errorsReduced, 2, MPI_DOUBLE, MPI_SUM, globalComm);
		for(i=0;i<2;i++)
			errors[i] = sqrt(errorsReduced[i]);
		if(iAmOutRank){
			OutPut(endl);
			OutPut( "L2: " << sqrt(errorsReduced[0]) << endl);
			OutPut( "H1-semi: " << sqrt(errorsReduced[1]) << endl);
		}
		#else
		OutPut(endl);
		OutPut( "L2: " << errors[0] << endl);
		OutPut( "H1-semi: " << errors[1] << endl);
		#endif
	} // if(TDatabase::ParamDB->MEASURE_ERRORS)
}
