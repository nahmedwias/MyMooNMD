#include <NSE3D.h>
#include <Database.h>



NSE3D::SystemPerGrid::SystemPerGrid(const Example_NSE3D& example,
                                    TCollection& coll
#ifdef _MPI
                                    , int maxSubDomainPerDof
#endif
) :  velocitySpace_(&coll, (char*)"u", (char*)"nse3d velocity", example.get_bc(0), //bd cond at 0 is x velo bc
                    TDatabase::ParamDB->VELOCITY_SPACE),
     pressureSpace_(&coll, (char*)"p", (char*)"nse3d pressure", example.get_bc(3), //bd condition at 3 is pressure bc
                    TDatabase::ParamDB->PRESSURE_SPACE),
     //matrix_(velocitySpace_, pressureSpace_, example.get_bd()), //TODO: use new BlockFEMatrix!
     rhs_(matrix_, true),
     solution_(matrix_, false),
     u_(&velocitySpace_, (char*)"u", (char*)"u", solution_.block(0),
        solution_.length(0), 3),
     p_(&pressureSpace_, (char*)"p", (char*)"p", solution_.block(3),
        solution_.length(3))

#ifdef _MPI
     //default construct parallel infrastructrue, will be reset in the body
     , parMapperVelocity_(),
     parMapperPressure_(),
     parCommVelocity_(),
     parCommPressure_()
#endif

{
#ifdef _MPI

  if(TDatabase::ParamDB->NSTYPE == 1 || TDatabase::ParamDB->NSTYPE == 3)
  {
    ErrThrow("Parallel solve needs correctly stored BT-blocks, so NSTYPE 2 or 4.");
  }

  velocitySpace_.SetMaxSubDomainPerDof(maxSubDomainPerDof);
  pressureSpace_.SetMaxSubDomainPerDof(maxSubDomainPerDof);

  //Must be reset here, because feSpace needs special treatment
  // This includes copy assignment - all because there is no good
  // way to communicate Maximum number of subdomains per dof to FESpace...
  parMapperVelocity_ = TParFEMapper3D(3, &velocitySpace_, //magic number
                                      matrix_.block(0,0)->GetRowPtr(), // The parMapper for velocity expects row and column
                                      matrix_.block(0,0)->GetKCol());  // pointer from the A matrices
  parMapperVelocity_ = TParFEMapper3D(1, &pressureSpace_,
                                      matrix_.block(0,3)->GetRowPtr(), // The parMapper for velocity expects row and column
                                      matrix_.block(0,3)->GetKCol());  // pointer from the BTransposed matrices - I don't know why exactly!
  // TODO The upper thing will lead to extreme trouble if there are no real BT blocks,
  // i.e. NSTYPE 1 or 3. Catch that case in an input-check method AND in the constructor.

  parCommVelocity_ = TParFECommunicator3D(&parMapperVelocity_);
  parCommPressure_ = TParFECommunicator3D(&parMapperPressure_);
#endif
}


NSE3D::NSE3D(std::list<TCollection* > collections, const Example_NSE3D& example
#ifdef _MPI
             , int maxSubDomainPerDof
#endif
) : systems_(), example_(example), multigrid_(nullptr)
{
  // start with only non-multigrid case

  // Check at least if the collections list contains exactly one Collection.
  if(collections.size() != 1 )
  {
    ErrThrow("Non-multigrid: Expected exactly one collection!");
  }

  // Get the one given collection.
  TCollection& cellCollection = *collections.front();

#ifdef _MPI
      // create finite element space and function, a matrix, rhs, and solution
      systems_.emplace_back(example_, cellCollection, maxSubDomainPerDof);
#else
      // create finite element space and function, a matrix, rhs, and solution
      systems_.emplace_back(example_, cellCollection);
#endif

}


