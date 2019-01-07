#include "TimeConvectionDiffusionPOD.h"
#include "Database.h"

#ifdef __2D__
#include "Assemble2D.h"
#include "SquareMatrix2D.h"
#else
#include "Assemble3D.h"
#include "SquareMatrix3D.h"
#endif

template <int d>
ParameterDatabase TimeConvectionDiffusionPOD<d>::default_tcd_pod_database()
{
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TimeConvectionDiffusionPOD parameter database");
  
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(ParameterDatabase::default_output_database());
  db.merge(AlgebraicFluxCorrection::default_afc_database());
  db.merge(LocalAssembling<d>::default_local_assembling_database(), true);
  db.merge(TimeDiscretization::default_TimeDiscretization_database(), true);
  return db;
}
/* ************************************************************************* */
template <int d>
TimeConvectionDiffusionPOD<d>::TimeConvectionDiffusionPOD(
    TCollection& coll, const ParameterDatabase& param_db)
 : TimeConvectionDiffusionPOD<d>(coll,
                                 param_db, Example_TimeCD(param_db))
{
}
/* ************************************************************************* */
template <int d>
TimeConvectionDiffusionPOD<d>::TimeConvectionDiffusionPOD(
  TCollection& coll, const ParameterDatabase& param_db, 
  const Example_TimeCD& ex)
 : POD(param_db), 
 fe_space(&coll, "space", "time_cd_pod space", ex.get_bc(0),
          TDatabase::ParamDB->ANSATZ_ORDER),
 db(param_db), time_stepping_scheme(param_db), example(ex),
 outputWriter(param_db)
{
#ifdef __3D__
  this->gramian_matrix = BlockFEMatrix::CD3D(fe_space);
#else
  this->gramian_matrix = BlockFEMatrix::CD2D(fe_space);
#endif
  this->pod_mode = BlockVector(this->gramian_matrix, false);
  this->fe_function = FEFunction(&this->fe_space, (char*)"c", (char*)"c",
                     this->pod_mode.get_entries(), this->pod_mode.length());
  
  this->set_parameters();
  
  this->output_problem_size_info();
  
  this->outputWriter.add_fe_function(&this->fe_function);
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusionPOD<d>::set_parameters()
{
  //set problem_type to Time_CD if not yet set
  if(!rom_db["problem_type"].is(2))
  {
    if (rom_db["problem_type"].is(0))
    {
      rom_db["problem_type"] = 2;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to Time_CD."
          "It is now reset to the correct value for Time_CD (=2).");
      rom_db["problem_type"] = 2;
    }
  }

  // an error when using ansatz order 0
  if(TDatabase::ParamDB->ANSATZ_ORDER == 0)
  {
    throw std::runtime_error("Ansatz order 0 is no use in convection diffusion "
        "reaction problems! (Vanishing convection and diffusion term).");
  }

  // For assembling of the Gramian matrix for the computation of
  // the POD basis only the "Galerkin" part of the matrices is needed.
  // TDatabase::ParamDB->DISCTYPE = GALERKIN;
  if(!rom_db["space_discretization_type"].is("galerkin"))
    rom_db["space_discretization_type"] = "galerkin";
}
/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusionPOD<d>::assemble_gramian()
{
  Output::print<2>("Assembling the gramian matrix for POD computation...");
  LocalAssembling_type type;
  if (rom_db["pod_inner_product"].get<std::string>() == "l2")
  {
    type = LocalAssembling_type::TCDMassOnly;
  }
  else if (rom_db["pod_inner_product"].get<std::string>() == "eucl")
  {
    Output::print<1>("For given parameter 'pod_inner_product' no assembling needed.");
    return;
  }
  FEFunction *feFunction = &this->fe_function;
  LocalAssembling<d> la(this->db, type, &feFunction, example.get_coeffs());
  
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  int nFESpaces = 1;
  const FESpace *fe_space = &this->fe_space;
  auto block = this->gramian_matrix.get_blocks_uniquely()[0].get();
  
  int nsqMat = 1;
  SquareMatrixD *sqMat[nsqMat];
  sqMat[0] = reinterpret_cast<SquareMatrixD*>(block);
  sqMat[0]->reset();
  int nRhs = 0;
  
  auto * bound_cond = this->fe_space.get_boundary_condition();
  auto * bound_val = this->example.get_bd(0);
  
#ifdef __2D__
  Assemble2D(
#else
  Assemble3D(
#endif
    nFESpaces, &fe_space, nsqMat, sqMat, 0, nullptr, nRhs, nullptr, 
    &fe_space, &bound_cond, &bound_val, la);
}
/* ************************************************************************* */
template <int d> 
void TimeConvectionDiffusionPOD<d>::output_problem_size_info() const
{
  const FESpace* space = &this->fe_space;
  TCollection *coll = space->GetCollection();
#ifndef _MPI
  double hMin, hMax;
  coll->GetHminHmax(&hMin, &hMax);
  int n_cells = coll->GetN_Cells();
  int n_dof = space->GetN_DegreesOfFreedom();
#else // _MPI
  
#endif
  {
    Output::stat("TimeConvectionDiffusion",
                 "Mesh data and problem size on finest grid");
    Output::dash("n cells     : ", setw(13), n_cells);
    Output::dash("h(min, max) : ", setw(13), hMin, " ", setw(13), hMax);
    Output::dash("n dofs      : ", setw(13), n_dof);
#ifndef _MPI
    Output::dash("n dof active: ", setw(13), space->GetActiveBound());
#endif
  }
}
/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusionPOD<d>::compute_pod_basis()
{
  if(this->rom_db["pod_inner_product"] != "eucl")
  {
    this->assemble_gramian();
    //cout<<"norm: "<< gramian_matrix.get_blocks().at(0)->GetNorm()<<endl;exit(0);
    POD::set_gramian(this->gramian_matrix.get_combined_matrix());
  }
  POD::compute_basis();

  if(rom_db["pod_inner_product"].get<std::string>() == "eucl")
  {
    Output::print<1>("POD w.r.t euclidean inner product.");
    POD::write_pod(rom_db["pod_basename"].get<std::string>()+"eucl.");
  }
  else if(rom_db["pod_inner_product"].get<std::string>() == "l2")
  {
    Output::print<1>("POD w.r.t L2 inner product.");
    POD::write_pod(rom_db["pod_basename"].get<std::string>()+"l2.");
  }
  else
  {
    ErrThrow("Unknown value of parameter 'pod_inner_product'. Choose eucl or l2.");
  }
}
/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusionPOD<d>::output()
{
  if (this->fe_space.GetN_DegreesOfFreedom() != POD::get_basis().size1())
  {
    ErrThrow("Current FE space dimension does not coincide with the number "
             "of dof of POD basis.\n"
             "Dimension of FE space : ", this->fe_space.GetN_DegreesOfFreedom(),
             "\nDOF of POD basis      : ", POD::get_basis().size1());
  }

  /* write averages of snapshots into a vtk-file */
  if(this->rom_db["pod_fluct"])
  {
    for(int j=0; j< this->fe_space.GetN_DegreesOfFreedom(); ++j)
    {
      this->pod_mode.get_entries()[j] = POD::get_snaps_avr()(j);
    }
    Output::print<1>("Writing vtk for snapshots average into ",
                      this->rom_db["output_basename"], 0 ,".vtk");
    this->outputWriter.write(0);
  }

  for(int i=0; i < std::min(POD::get_rank(),10); ++i)
  {
    for(int j=0; j<this->fe_space.GetN_DegreesOfFreedom(); ++j)
    {
      this->pod_mode.get_entries()[j] = POD::get_basis()(j,i);
    }
    Output::print<1>("Writing vtk for POD basis into ",
                     this->rom_db["output_basename"], i+1 ,".vtk");
    this->outputWriter.write(i+1);
  }
}

#ifdef __3D__
template class TimeConvectionDiffusionPOD<3>;
#else
template class TimeConvectionDiffusionPOD<2>;
#endif