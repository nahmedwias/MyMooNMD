/**
 * @brief A test program for a time dependent CD3D problem
 * algorithm "linear Crank-Nicolson FEM-FCT" as implemented
 * in the namespace AlgebraicFluxCorrection.
 * Uses LeVeque's rotating bodies example.
 *
 * @date 2016/01/13
 * @author Clemens Bartsch, Naveed Ahmed
 */
#include <AlgebraicFluxCorrection.h>
#include <Time_CD3D.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <TimeDiscRout.h>
#include <MainUtilities.h>

#ifdef _MPI
#include <mpi.h>
#include <MeshPartition.h>
double bound = 0;
double timeC = 0;
#endif

// Those are L2, H1 and L^inf the norms of the target solution in space,
// at the given times. The test must reproduce these in order to pass.
std::vector<std::vector<double>> target_norms =
{
  {0.004872777812, 0.115936124, 0.1460876845}, //t=0.1
  {0.01206790999, 0.2407316993, 0.2966152375}, //t=0.2
  {0.02082420123, 0.3809191143, 0.4416315506}, //t=0.3
  {0.03034782753, 0.5272384081, 0.5759509038}, //t=0.4
  {0.04006347334, 0.6711927872, 0.6960808079}, //t=0.5
  {0.04955433812, 0.8067095407, 0.7990662275}, //t=0.6
  {0.05848061731, 0.9289673979, 0.8823724209}, //t=0.7
  {0.06655716561, 1.034141517, 0.9439503241},  //t=0.8
  {0.07354159893, 1.119053591, 0.9822839671},  //t=0.9
  {0.07923718531, 1.181488119, 0.9965685447}   //t=1
};

void check_solution_norms(Time_CD3D &tcd, int m)
{

  MultiIndex3D allDerivatives[4] = { D000, D100, D010, D001 };
  TAuxParam3D aux(1, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
  std::array<double, 5> locError = {};
  const TFESpace3D* space = tcd.get_function().GetFESpace3D();

  tcd.get_function().GetErrors(tcd.get_example().get_exact(0), 4, allDerivatives, 2,
                          L2H1Errors, tcd.get_example().get_coeffs(), &aux, 1, &space,
                          locError.data());

  if(m%10 == 0)
  {
    //Output::print("Norms in step ", m); //This is how the reference values were produced.
    //Output::dash(std::setprecision(10), locError[0], ", ", locError[1], ", ", locError[2]);

    int control_step = m/10 - 1;
    double tol = 1e-9;

    if( fabs(locError[0] - target_norms[control_step][0]) > tol )
      ErrThrow("L2 norm at timestep ", TDatabase::TimeDB->CURRENTTIME,  " is not correct, ",
               locError[0], " != ", target_norms[control_step][0]);

    if( fabs(locError[1]  - target_norms[control_step][1]) > tol )
      ErrThrow("H1 norm at timestep ", TDatabase::TimeDB->CURRENTTIME,  " is not correct, ",
               locError[1], " != ", target_norms[control_step][1]);

    if( fabs(locError[2]  - target_norms[control_step][2]) > tol )
      ErrThrow("L^inf norm at timestep ", TDatabase::TimeDB->CURRENTTIME,  " is not correct, ",
               locError[2], " != ", target_norms[control_step][2]);

    Output::print("Solution norms checked succesfully.");
  }
}


 //test crank-nicolson linear fem-fct scheme
int main(int argc, char* argv[])
{
#ifdef _MPI
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  TDatabase::ParamDB->Comm = comm;
#endif
  {
    TDatabase Database;
    TFEDatabase3D FEDatabase;
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(Example3D::default_example_database());
    db.merge(ParameterDatabase::default_output_database());
    db["problem_type"] = 1;
    db["example"] = 0;
    db["diffusion_coefficient"] = 1e-6;

    db.merge(AlgebraicFluxCorrection::default_afc_database(), true);
    db["algebraic_flux_correction"].set("fem-fct-cn");

    db["space_discretization_type"] = "galerkin";
    TDatabase::ParamDB->ANSATZ_ORDER=1;

    TDatabase::TimeDB->STARTTIME = 0;
    TDatabase::TimeDB->ENDTIME = 1;
    TDatabase::TimeDB->TIMESTEPLENGTH = 0.01;
    SetTimeDiscParameters(0);

    db.add("boundary_file", "Default_UnitCube", "");
    db.add("geo_file", "Default_UnitCube_Hexa", "",
           {"Default_UnitCube_Hexa", "Default_UnitCube_Tetra"});
    TDatabase::ParamDB->DRIFT_Z = 1;
    db.add("refinement_n_initial_steps",(size_t) 4,"",(size_t) 0, (size_t) 5);
    db.add("solver_type", "direct", "", {"direct"});

//    db["output_write_vtk"] = true;
//    db["output_vtk_directory"].set_range(std::set<std::string>({std::string("."),std::string("VTK")}));
//    db["output_vtk_directory"].set("VTK");

    TDomain domain(db);
#ifdef _MPI
  int maxSubDomainPerDof = 0;
#endif
    std::list<TCollection* > gridCollections
      = domain.refine_and_get_hierarchy_of_collections( db
#ifdef _MPI
    , maxSubDomainPerDof
#endif
          );

    // example object
    Example_TimeCD3D example_obj(db);
    // tcd3d system object
#ifdef _MPI
    Time_CD3D tcd(gridCollections, db, example_obj, maxSubDomainPerDof);
#else
    Time_CD3D tcd(gridCollections, db, example_obj);
#endif


    tcd.assemble_initial_time();

    int step = 0;
    //int imag = 0;

    while(TDatabase::TimeDB->CURRENTTIME <
      TDatabase::TimeDB->ENDTIME-1e-10)
    {
      step ++;
      TDatabase::TimeDB->INTERNAL_STARTTIME
         = TDatabase::TimeDB->CURRENTTIME;
      SetTimeDiscParameters(1);

      double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;

      Output::print<1>("\nCURRENT TIME: ",
           TDatabase::TimeDB->CURRENTTIME);
      tcd.assemble();
      tcd.solve();

      //tcd.output(step,imag);

      check_solution_norms(tcd, step);
    }

  }
#ifdef _MPI
  MPI_Finalize();
#endif
}



