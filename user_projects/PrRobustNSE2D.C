#include <PrRobustNSE2D.h>
#include <Assemble2D.h>
#include <NSE2D.h>
#include <Domain.h>
#include <Database.h>
/** ************************************************************************ */
PrRobustNSE2D::SystemPerGrid::SystemPerGrid(const Example_NSE2D& example, 
TCollection& coll, unsigned int order, const TFESpace2D& space, 
                   const TFESpace2D& presSpace)
: projection_space(space.GetCollection(), (char*)"u", 
                   (char*)"projection velocity", example.get_bc(0),
                   TDatabase::ParamDB->PROJECTION_SPACE, nullptr),
 projMat(&space, &projection_space, &presSpace),// V_h, X_h
 rhsXh(projMat, false)
{
  // nothing to do here
}
/** ************************************************************************ */
PrRobustNSE2D::PrRobustNSE2D(const TDomain& domain, 
       const Example_NSE2D& _example, unsigned int reference_id)
 : NSE2D(domain, _example, reference_id), Systems()
{
  // create the collection of cells from the domain 
  TCollection *coll = domain.GetCollection(It_Finest,0,reference_id);

  int proj_order = TDatabase::ParamDB->PROJECTION_METHOD;
  
  this->Systems.emplace_back(example, *coll, proj_order, 
                             this->NSE2D::get_velocity_space(), 
                             this->NSE2D::get_pressure_space());
  unsigned int nDofProj = 
     this->get_projection_space().GetN_DegreesOfFreedom();
  Output::print<1>("dof        : ",  setw(10), nDofProj);
  Output::print<1>("active dof : ",  setw(10), this->get_projection_space().GetActiveBound());
}
/** ************************************************************************ */
void PrRobustNSE2D::assembleMatrixRhs()
{
  // assembling of the matrices A and B's only
  TFEFunction2D *fe_functions[3] = 
     { this->NSE2D::systems.front().u.GetComponent(0), 
       this->NSE2D::systems.front().u.GetComponent(1), 
       &this->NSE2D::systems.front().p };
  
  LocalAssembling2D la(NSE2D_Galerkin,fe_functions, 
                       this->NSE2D::get_example().get_coeffs());
  // assemble the system matrix
  this->NSE2D::systems.front().matrix.Assemble(la, 
                                      this->NSE2D::systems.front().rhs);
  if(TDatabase::ParamDB->DISCTYPE != RECONSTRUCTION)
    return;
  // assemble the right hand side 
  this->assemblerhs();

  this->Systems.front().projMat.apply(
                this->Systems.front().rhsXh,this->NSE2D::systems.front().rhs);  
}
/** ************************************************************************ */
void PrRobustNSE2D::assemblerhs()
{
  TFEFunction2D *fe_functions[2] = 
     { this->NSE2D::systems.front().u.GetComponent(0), 
       this->NSE2D::systems.front().u.GetComponent(1) };
  // local assemble      
  LocalAssembling2D larhs(RECONSTR_GALERKIN_Rhs, fe_functions,
                          this->NSE2D::get_example().get_coeffs());
  
  BoundCondFunct2D *bc[2] = { this->get_projection_space().GetBoundCondition(),
                              this->get_projection_space().GetBoundCondition() };
  
  BoundValueFunct2D * const * const BoundValue 
       = this->NSE2D::get_example().get_bd();
  BoundValueFunct2D *bv[2] = {BoundValue[2], BoundValue[2] };
  
  const TFESpace2D *v_space = &this->get_projection_space();
  const TFESpace2D * pointer_to_space[1] = { v_space };
  
  // assemble the right hand side separately  
  this->Systems.front().rhsXh.reset();
  
  // setting the space for sign in GetSignOfThisDOF();
  unsigned int vel_space = TDatabase::ParamDB->VELOCITY_SPACE;
  TDatabase::ParamDB->VELOCITY_SPACE = TDatabase::ParamDB->PROJECTION_SPACE;
  
  double *rhs_blocks[1] = { this->Systems.front().rhsXh.get_entries() };
  
  // assemble right hand side only
  Assemble2D_VectFE(1, pointer_to_space, 0, nullptr, 0, 
                    nullptr, 1, rhs_blocks, pointer_to_space, larhs, 
                    bc, bv);
  TDatabase::ParamDB->VELOCITY_SPACE = vel_space;
}
/** ************************************************************************ */
void PrRobustNSE2D::solve()
{
  this->NSE2D::solve();
  this->NSE2D::output(0);  
}

/** ************************************************************************ */