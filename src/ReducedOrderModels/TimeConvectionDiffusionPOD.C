#include "TimeConvectionDiffusionPOD.h"
#include "Database.h"
#include <sys/stat.h>


#ifdef __2D__
#include "Assemble2D.h"
#include "SquareMatrix2D.h"
#else
#include "Assemble3D.h"
#include "SquareMatrix3D.h"
#endif

#ifdef _MPI
#include "ParFECommunicator3D.h"
#endif

/* ************************************************************************** */
template <int d>
ParameterDatabase TimeConvectionDiffusionPOD<d>::set_pod_basis_database(
                                              const ParameterDatabase& param_db)
{
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("Parameter database for POD basis computation");

  db.merge(LocalAssembling<d>::default_local_assembling_database(), true);

  db.merge(ParameterDatabase::default_output_database(), true);
  db.merge(POD::default_pod_database(), true);
  db.merge(param_db, false);

  // change name to avoid overwritting standard FEM output
  // TODO: unification between naming system (outfile, output_basename,
  //       snaps_basename, pod_basename)
  if(param_db["write_snaps"])
  {
    std::string output_pod_basename = db["output_basename"].get<std::string>()
                                      + "_pod";
    db["output_basename"].set(output_pod_basename, false);
  }
  return db;
}

/* ************************************************************************** */
template <int d>
TimeConvectionDiffusionPOD<d>::TimeConvectionDiffusionPOD(TCollection& coll,
                                              const ParameterDatabase& param_db)
  : TimeConvectionDiffusionPOD<d>(coll, param_db, Example_TimeCD(param_db))
{
}

/* ************************************************************************** */
template <int d>
TimeConvectionDiffusionPOD<d>::TimeConvectionDiffusionPOD(TCollection& coll,
                                              const ParameterDatabase& param_db, 
                                              const Example_TimeCD&    ex)
  : POD(param_db), 
    fe_space(&coll, "space", "time_cd_pod space", ex.get_bc(0),
             TDatabase::ParamDB->ANSATZ_ORDER),
    db(TimeConvectionDiffusionPOD<d>::set_pod_basis_database(param_db)),
    time_stepping_scheme(param_db),
    example(ex),
    outputWriter(db)
{
  // create directory db["pod_directory"]
  std::string directory_name = this->db["pod_directory"].get<std::string>();
  mkdir(directory_name.c_str(), 0777);

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

/* ************************************************************************** */
template <int d>
void TimeConvectionDiffusionPOD<d>::set_parameters()
{
  //set problem_type to Time_CD if not yet set
  if(!db["problem_type"].is(2))
  {
    if (db["problem_type"].is(0))
    {
      db["problem_type"] = 2;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to "
                      "Time_CD. It is now reset to the correct value for "
                      "Time_CD (=2).");
      db["problem_type"] = 2;
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
  if(!db["space_discretization_type"].is("galerkin"))
    db["space_discretization_type"] = "galerkin";
}

/* ************************************************************************** */
template <int d>
void TimeConvectionDiffusionPOD<d>::assemble_gramian()
{
  Output::print<1>("[pod_inner_product = L2]: "
                   "Assembling the gramian matrix for POD computation...");
  //LocalAssembling_type type;
  //type = LocalAssembling_type::TCDMassOnly;

  FEFunction *feFunction = &this->fe_function;
  LocalAssembling<d> la(this->db, LocalAssembling_type::TCDMassOnly,
                        &feFunction, example.get_coeffs());

  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  int nFESpaces = 1;
  const FESpace *fe_space = &this->fe_space;
  /**
     @attention
     we get block [1] because the TCD Mass Matrix functions
     assembles Mass as block[1] (block[0] is reserved for stiffness)
  **/
  auto block = this->gramian_matrix.get_blocks_uniquely()[0].get();

  int nsqMat = 1;
  SquareMatrixD *sqMat[nsqMat];
  sqMat[0] = reinterpret_cast<SquareMatrixD*>(block);
  Output::print("HERE");
  sqMat[0]->reset();
  Output::print("HERE");
  int nRhs = 0;

  auto * bound_cond = this->fe_space.get_boundary_condition();
  auto * bound_val = this->example.get_bd(0);
  Output::print("HERE");
#ifdef __2D__
  Assemble2D(
#else
  Assemble3D(
#endif
    nFESpaces, &fe_space, nsqMat, sqMat, 0, nullptr, nRhs, nullptr,
    &fe_space, &bound_cond, &bound_val, la);
  sqMat[0]->write("test.m");

}

/* ************************************************************************** */
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
  int root = 0; // root process number
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int n_local_master_cells = coll->GetN_OwnCells();
  int n_cells;
  MPI_Reduce(&n_local_master_cells, &n_cells, 1, MPI_DOUBLE, MPI_SUM, root,
             MPI_COMM_WORLD);

  double local_hmin, local_hmax;
  coll->GetHminHmax(&local_hmin, &local_hmax);
  double hMin, hMax;
  MPI_Reduce(&local_hmin, &hMin, 1, MPI_DOUBLE, MPI_MIN, root, MPI_COMM_WORLD);
  MPI_Reduce(&local_hmax, &hMax, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);

  auto par_comm = space->get_communicator();
  int n_dof  = par_comm.get_n_global_dof();
  if(my_rank == root)
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

/* ************************************************************************** */
template <int d>
void TimeConvectionDiffusionPOD<d>::compute_pod_basis()
{
  if(this->db["pod_inner_product"].is("L2"))
  {
    this->assemble_gramian();
    //cout<<"norm: "<< gramian_matrix.get_blocks().at(0)->GetNorm()<<endl;exit(0);
    POD::set_gramian(this->gramian_matrix.get_combined_matrix());
  }
  POD::compute_basis();

  std::string basename = this->db["pod_basename"];
  if(this->db["pod_inner_product"].is("euclidean"))
  {
    Output::print<1>("POD w.r.t euclidean inner product.");
    POD::write_pod(basename + ".euclidean.");
  }
  else if(this->db["pod_inner_product"].is("L2"))
  {
    Output::print<1>("POD w.r.t L2 inner product.");
    POD::write_pod(basename +".L2.");
  }
  else
  {
    ErrThrow("Unknown value of parameter 'pod_inner_product'. "
             "Choose euclidean or L2.");
  }
}

/* ************************************************************************** */
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
  if(this->db["pod_fluctuations_only"])
  {
    for(int j=0; j< this->fe_space.GetN_DegreesOfFreedom(); ++j)
    {
      this->pod_mode.get_entries()[j] = POD::get_snaps_avr()(j);
    }
    Output::print<1>("Writing vtk for snapshots average into ",
                      this->db["output_basename"], 0 ,".vtk");
    this->outputWriter.write(0);
  }

  for(int i=0; i < std::min(POD::get_rank(),10); ++i)
  {
    for(int j=0; j<this->fe_space.GetN_DegreesOfFreedom(); ++j)
    {
      this->pod_mode.get_entries()[j] = POD::get_basis()(j,i);
    }
    Output::print<1>("Writing vtk for POD basis into ",
                     this->db["output_basename"], i+1 ,".vtk");
    this->outputWriter.write(i+1);
  }
}

#ifdef __3D__
template class TimeConvectionDiffusionPOD<3>;
#else
template class TimeConvectionDiffusionPOD<2>;
#endif
