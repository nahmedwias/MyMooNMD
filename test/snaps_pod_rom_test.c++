/** ****************************************************************************
 *
 * @name  snaps_pod_rom_test
 * @brief A test program to test snapshot Collector, POD-basis computation
 *
 * @autor Baptiste Moreau
 * @date  18.03.2019
*******************************************************************************/

#include <Database.h>
#include <Domain.h>
#include <FEDatabase2D.h>
#include <SnapshotsCollector.h>
#include <TimeConvectionDiffusion.h>
#include <TimeConvectionDiffusionPOD.h>
#include <TimeConvectionDiffusionROM.h>
#include <TimeDiscretizations.h>
#include <TimeDiscRout.h>
#include <MainUtilities.h>
#include <ConvDiff.h>
#include <AuxParam2D.h>
#include <list>

/** ***************************************************************************/
void test_mean(ublas::vector<double> mean_snaps)
{
  const int nb_mean = 25;
  double    eps     = 1e-10;

  double mean_exp[nb_mean] = {1.20975262773, -0.0139937900321, 0.00213728930742,
                              0.0153169167998, -1.21218379178, -0.0139937900321,
                              1.20975262773, 0.0153169167998, -1.20603517152,
                              0., 0., 0., 0., 0., 0., -8.1643119943e-17, 0.,
                              -9.99839855107e-33, 1.99967971022e-32,
                              8.1643119943e-17, 8.1643119943e-17,
                              -9.99839855107e-33, 0., 0., -8.1643119943e-17};


  for(int i=0 ; i<nb_mean ; i++)
  {
      if( fabs(mean_snaps[i] - mean_exp[i]) > eps )
      {
        ErrThrow("test snapshot with tcd2d: the average of snapshots is not "
                 "correct ", std::setprecision(12), mean_snaps[i],
                 " != ", std::setprecision(12), mean_exp[i]);
      }
  }
}

/** ***************************************************************************/
void test_eigs(const double* eigs_pod)
{
  const int nb_eigs = 2;
  double    eps     = 1e-10;

  double eigs_exp[nb_eigs] = {0.00764769512177, 8.85214448497e-07};

  for(int i=0 ; i<nb_eigs ; i++)
  {
      if( fabs(eigs_pod[i] - eigs_exp[i]) > eps )
      {
        ErrThrow("test POD with tcd2d: the eigenvalue is not correct ",
                 std::setprecision(12), eigs_pod[i],
                 " != ", std::setprecision(12), eigs_exp[i]);
      }
  }
}

/** ***************************************************************************/
int remove_tmp_test_files(ParameterDatabase db)
{
  int status = 0;
  std::string snaps_dir_name = db["snaps_directory"];
  std::string pod_dir_name   = db["pod_directory"];
  std::string snaps_name    = snaps_dir_name + "/"
                             + db["snaps_basename"].get<std::string>();
  std::string pod_basename  = pod_dir_name + "/"
                            + db["pod_basename"].get<std::string>() + "."
                            + db["pod_inner_product"].get<std::string>() + ".";
  std::string pod_name_pod  = pod_basename + "pod";
  std::string pod_name_mean = pod_basename + "mean";
  std::string pod_name_eigs = pod_basename + "eigs";

  status |= (std::remove(snaps_name.c_str()) == 0) ? 0 : 1;
  status |= (std::remove(snaps_dir_name.c_str()) == 0) ? 0 : 1;

  status |= (std::remove(pod_name_pod.c_str()) == 0) ? 0 : 1;
  status |= (std::remove(pod_name_mean.c_str()) == 0) ? 0 : 1;
  status |= (std::remove(pod_name_eigs.c_str()) == 0) ? 0 : 1;
  status |= (std::remove(pod_dir_name.c_str()) == 0) ? 0 : 1;

  return status;
}

/** ***************************************************************************/
ParameterDatabase set_database()
{
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(Example2D::default_example_database());
    db.merge(LocalAssembling2D::default_local_assembling_database());
    db.merge(TimeDiscretization::default_TimeDiscretization_database());
    db["example"] = 0;
    db["problem_type"] = 2;
    db["reynolds_number"] = 1;

    db["time_discretization"] = "backward_euler";
    db["time_step_length"] = 0.01;
    db["time_end"]=0.1;

    db["write_snaps"]=true;
    db["compute_POD_basis"]=true;
//    db["ROM_method"]=true;

    // declaration of databases
    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles"});
    db.add("refinement_n_initial_steps", (size_t) 2,"");

    // declaration of parameters used for snapshots, POD and ROM
    db.add("snaps_directory", "snapshots_tmp", "");
    db.add("snaps_basename", "snaps_tmp", "");
    db.add("steps_per_snap", (size_t) 5, "");
    db.add("pod_directory", "pod_basis_tmp", "");
    db.add("pod_basename", "pod_tmp", "");
    db.add("pod_rank", (size_t) 0, "");
    db.add("pod_fluctuations_only", true, "");
    db.add("pod_inner_product", "L2", "");

    // test direct solver with Galerkin
    Output::print("\n\nTesting Galerkin\n");
    db.add("solver_type", "direct", "", {"direct", "petsc"});
    db["space_discretization_type"] = "galerkin";

    return db;
}

/** ***************************************************************************/
int main(int, char**)
{
  //============================================================================
  // test with TCD2D
  //============================================================================
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  ParameterDatabase db = set_database();
  TDatabase::ParamDB->ANSATZ_ORDER=1;
  TDomain domain(db);

  // refine grid
  domain.refine_and_get_hierarchy_of_collections(db);

  TimeConvectionDiffusion<2> tcd(domain, db);

  TimeDiscretization& tss = tcd.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.set_time_disc_parameters();

  TDatabase::TimeDB->CURRENTTIME = tss.get_start_time();

  tcd.assemble_initial_time();

  //----------------------------------------------------------------------------
  // solve the PDE using FEM and acquire snapshots
  //----------------------------------------------------------------------------
  if(db["write_snaps"])
  {
    // initialize snapshot writer (only needed if parmoon_db["write_snaps"])
    // the solutions are stored in a special collector which can be later read
    // by the computation of POD the basis
    SnapshotsCollector snaps(db);

    // store initial condition as snapshot
    snaps.write_data(tcd.get_solution());

    while(!tss.reached_final_time_step())
    {
      tss.current_step_++;
      SetTimeDiscParameters(1);

      tss.current_time_ += tss.get_step_length();;
      TDatabase::TimeDB->CURRENTTIME += tss.get_step_length();

      Output::print<1>("\nCURRENT TIME: ", tss.current_time_);
      tcd.assemble();
      tcd.solve();

      // write the snapshots
      if(db["write_snaps"])
      {
        snaps.write_data(tcd.get_solution(), tss.current_step_);
      }
    }
  }

  //----------------------------------------------------------------------------
  // compute the POD basis from snapshots
  //----------------------------------------------------------------------------
  auto collections            = domain.get_grid_collections();
  TCollection& cellCollection = *collections.front();
  TimeConvectionDiffusionPOD<2> tcd_pod(cellCollection, db);

  if(db["compute_POD_basis"])
  {
    tcd_pod.compute_pod_basis();
  }

  //----------------------------------------------------------------------------
  // test snapshots average and eigenvalues from POD
  // delete temporary files used for this test
  //----------------------------------------------------------------------------
  test_mean(tcd_pod.get_snaps_avr());
  test_eigs(tcd_pod.get_pod_eigs());

  if(remove_tmp_test_files(db) != 0)
  {
    ErrThrow("snaps_pod_rom_test: Error deleting temporary files.");
  }

  return 0;
}

