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
#include <TimeDiscretizations.h>
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
  {0.005212129535, 0.1176258532, 0.1478578541}, //t=0.1
  {0.01289918241, 0.2492018215, 0.2996827932}, //t=0.2
  {0.02194896677, 0.3929355722, 0.4448496709}, //t=0.3
  {0.0316040805, 0.5398213516, 0.5789690547}, //t=0.4
  {0.0413625835, 0.6829024359, 0.6988235531}, //t=0.5
  {0.05083638957, 0.8166941799, 0.8014768829}, //t=0.6
  {0.05969792488, 0.9366684992, 0.8843946968}, //t=0.7
  {0.06766763184, 1.039072199, 0.9455355656},  //t=0.8
  {0.07451195658, 1.120909102, 0.9833945554},  //t=0.9
  {0.08003741409, 1.180117591, 0.9967426723}   //t=1
};

// These are target norms if NO flux correction scheme (or, indeed, any
// stabilization) is used. They can be used for debugging purposes, i.e.,
// comparison of MPI and SEQ runs.
std::vector<std::vector<double>> target_norms_none =
{
  {0.003943723865, 0.123457981, 0.1406667838},  //t=0.1
  {0.009479834582, 0.2386029396, 0.2887335298}, //t=0.2
  {0.01663589269, 0.3759851755, 0.436005843},   //t=0.3
  {0.02483201597, 0.5325839514, 0.5740274437},  //t=0.4
  {0.03351455826, 0.6991219501, 0.6972416763},  //t=0.5
  {0.04229832163, 0.868334923, 0.8028065696},   //t=0.6
  {0.05090244664, 1.034659303, 0.8887627679},   //t=0.7
  {0.05907653406, 1.193122896, 0.9668985986},   //t=0.8
  {0.06658406325, 1.339123213, 1.036639097},    //t=0.9
  {0.0732177335, 1.468659596, 1.080128}         //t=1
};


void check_solution_norms(Time_CD3D &tcd, int m)
{

  if(m%10 == 0)
  {
    MultiIndex3D allDerivatives[4] = { D000, D100, D010, D001 };
    TAuxParam3D aux(1, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, 0, nullptr);
    std::array<double, 5> locError = {};
    const TFESpace3D* space = tcd.get_function().GetFESpace3D();

    tcd.get_function().GetErrors(tcd.get_example().get_exact(0), 4, allDerivatives, 2,
                            L2H1Errors, tcd.get_example().get_coeffs(), &aux, 1, &space,
                            locError.data());

  #ifdef _MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<double> errorsReduced(4);
    MPI_Reduce(locError.data(), errorsReduced.data(), 2, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&locError[2], &errorsReduced[2], 1, MPI_DOUBLE, MPI_MAX, 0,
               MPI_COMM_WORLD);
    if(rank == 0)
    {//this is only performed on root - just root willl have the correct values
      locError[0] = sqrt(errorsReduced[0]);
      locError[1] = sqrt(errorsReduced[1]);
      locError[2] = errorsReduced[2];
  #endif // _MPI

    Output::print("Norms in step ", m); //This is how the reference values were produced.
    Output::dash(std::setprecision(10), locError[0], ", ", locError[1], ", ", locError[2]);

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

#ifdef _MPI
    }//end if rank == 0
#endif
  }
}


 //test crank-nicolson linear fem-fct scheme
int main(int argc, char* argv[])
{

  TDatabase Database;
  int rank = 0;
#ifdef _MPI
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  TDatabase::ParamDB->Comm = comm;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  {
    TFEDatabase3D FEDatabase;
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(Example3D::default_example_database());
    db.merge(ParameterDatabase::default_output_database());
    db.merge(LocalAssembling3D::default_local_assembling_database());
    db.merge(TimeDiscretization::default_TimeDiscretization_database());
    db["problem_type"] = 1;
    db["example"] = 0;
    db["diffusion_coefficient"] = 1e-6;

    db.merge(AlgebraicFluxCorrection::default_afc_database(), true);
    db["algebraic_flux_correction"].set("fem-fct-cn");

    db["space_discretization_type"] = "galerkin";
    TDatabase::ParamDB->ANSATZ_ORDER=1;

    db["time_end"] = 1;
    TDatabase::TimeDB->TIMESTEPLENGTH = 0.01;
    SetTimeDiscParameters(0);

    db.add("boundary_file", "Default_UnitCube", "");
    db.add("geo_file", "Default_UnitCube_Hexa", "",
           {"Default_UnitCube_Hexa", "Default_UnitCube_Tetra"});
    db.add("refinement_n_initial_steps",(size_t) 4,"",(size_t) 0, (size_t) 5);
    db.merge(Solver<>::default_solver_database(), true);
    db["solver_type"] = "iterative";
    db["iterative_solver_type"] = "fgmres";
    db["preconditioner"] = "jacobi";
    db["residual_tolerance"] = 1.0e-13;
    db["residual_reduction"] =  0.0;

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
    //tcd.output(step);
    double end_time = db["time_end"];
    while(TDatabase::TimeDB->CURRENTTIME < end_time-1e-10)
    {
      step ++;
      TDatabase::TimeDB->INTERNAL_STARTTIME
         = TDatabase::TimeDB->CURRENTTIME;
      SetTimeDiscParameters(1);

      double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;
      if(rank == 0)
        Output::print<1>("CURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
      tcd.assemble();
      tcd.solve();

      //tcd.output(step);
      check_solution_norms(tcd, step);
    }

  }
#ifdef _MPI
  MPI_Finalize();
#endif
}
