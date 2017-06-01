/**
 * A test for the grid transfer tool.
 *
 * The grid transfer tool is a tool that can be use for
 * transferring/interpolating/... FE functions from any grid and space to a
 * single FE space which is chosen at construction time. It is intelligent enough
 * to cache and re-use interpolation information (containment of interpolation
 * points) for all spaces which it was used for once.
 *
 * In this test, the tool is used to transfer functions from different FE spaces on a
 * coarse quad grid to a Q0 space with Dirichlet boundary conditions on a coarser
 * grid. This is the use case which I have in mind.
 * It is checked against a direct interpolation of the analytic functions.
 * NOTE: It shows to give exactly the same results, except for the case where the
 * 'original space' is a fine grid Q0, too, and the exact function is not a
 * polynomial of order leq 1. Why there is a difference in those cases is clear
 * to me - what is puzzling me, is: why is there no difference in the other cases??
 * ANSWER: It is about the cell mid point here. A "correct" Q0 grid transfer operation
 * needs a "correct" value at each cell mid point. This one is lost, when the analytical
 * solution is discretized with Q1. (or something like this...)
 */

#include <Database.h>
#include <Domain.h>
#include <GridTransferTool.h>
#include <MooNMD_Io.h>
#include <FEFunction2D.h>
#include <BlockVector.h>
#include <FEDatabase2D.h>
#include <Database.h>
#include <MainUtilities.h>
#include <PostProcessing2D.h>

void BoundConditionDirichlet(int BdComp, double t, BoundCond &cond)
{
   cond = DIRICHLET;
}

void BoundConditionNeumann(int BdComp, double t, BoundCond &cond)
{
   cond = NEUMANN;
}

void BoundConditionBoth(int BdComp, double t, BoundCond &cond)
{
  if(BdComp == 0 || BdComp == 1)
   cond = DIRICHLET;
  else
   cond = NEUMANN;
}


// some functions dummy given as an analytical representation
void analytic_function_1(double x, double y, double * values)
{
  values[0]=sin(Pi*x);
  values[1]=Pi*cos(Pi*x);
  values[2]=0;
  values[3]=0;
}
void analytic_function_2(double x, double y, double * values)
{
  values[0] = y*y*y - 2*y*x;
  values[1] = 2*y;
  values[2] = 3*y*y - 2*x;
  values[3] = -2;
}

double euclidean_distance(std::vector<double> a, std::vector<double> b)
{
  if(a.size() != b.size())
    ErrThrow("Vector sizes do not match!");

  //yeah I know I should use BLAS instead...
  double d = 0;

  for(int i = 0; i<a.size(); ++i)
    d += (a[i] - b[i])*(a[i] - b[i]);

  return std::sqrt(d);
}

int main(int argc, char **argv)
{
  // declaration of old databases (unfortunately still needed in ParMooN)
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  Output::setVerbosity(2);

  // new database
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(ParameterDatabase::default_output_database(),true);
  db["output_write_vtk"] = true;

  // set a few parameters in the database to construct a unit cube with one
  // cell. then refine uniformly twice
  db.add("boundary_file", "Default_UnitSquare", "");
  db.add("geo_file", "UnitSquare", "",
          {"UnitSquare","TwoTriangles"});
  db.add("refinement_n_initial_steps", (size_t) 4,"", (size_t) 0, (size_t) 5);
  // construct a domain object
  TDomain domain(db);
  //Refinement of grid in 3D for MPI ist different than for the sequential case in 2D of read_write_fe.
  // refine grid
  size_t n_ref = domain.get_n_initial_refinement_steps();
  TCollection* coarser_grid;
  for(unsigned int i = 0; i < n_ref; i++)
  {
    if(i == n_ref - 1)
      coarser_grid = domain.GetCollection(It_Finest, 0);
    domain.RegRefineAll();
  }
  TCollection* finer_grid = domain.GetCollection(It_Finest, 0);


  std::vector<size_t> ansatz_orders = {0, 1, 2, 3};
  std::vector<BoundCondFunct2D*> bound_conds({BoundConditionDirichlet,
                                              BoundConditionNeumann,
                                              BoundConditionBoth});

  std::vector<std::vector<const TFESpace2D*>> original_fe_spaces(4, {nullptr, nullptr, nullptr});

  for(size_t ans =0 ; ans < ansatz_orders.size() ; ++ans )
  {
    for(size_t bc = 0; bc < bound_conds.size() ;++bc)
    {
      original_fe_spaces[ans][bc] =
          new TFESpace2D(finer_grid, (char*) "",(char*) "", bound_conds[bc], ansatz_orders[ans],nullptr);
    }
  }

  TFESpace2D* target_space =
            new TFESpace2D(coarser_grid, (char*) "target_space",
                           (char*) "space to project into",
                           BoundConditionDirichlet, 0, nullptr);

  // Set up the grid transfer tool, a reference solution and
  // function and values into which to write the interpolation.
  GridTransferTool trafo(target_space, GridTransferType::MidPointEvaluation);

  std::vector<double> reference_1_vals(target_space->GetN_DegreesOfFreedom());
  TFEFunction2D reference_1(target_space, "coarse-grid-reference-1",
                            " ", &reference_1_vals.at(0), reference_1_vals.size());
  reference_1.Interpolate(analytic_function_1);

  std::vector<double> reference_2_vals(target_space->GetN_DegreesOfFreedom());
  TFEFunction2D reference_2(target_space, "coarse-grid-reference-2",
                            " ", &reference_2_vals.at(0), reference_2_vals.size());
  reference_2.Interpolate(analytic_function_2);

  std::vector<double> into_vals(target_space->GetN_DegreesOfFreedom());
  TFEFunction2D into_function(target_space, "function-into-which-to-interpolate",
                              " ", &into_vals.at(0), into_vals.size());

  // Now go through all the spaces, set up fe representations of the
  // analytic functions, and try to transfer them to the target space.
  for(size_t ans =0 ; ans < ansatz_orders.size() ; ++ans )
  {
    for(size_t bc = 0; bc < bound_conds.size() ;++bc)
    {
      //set up the "original" fe function
      std::vector<double> original_vals(original_fe_spaces[ans][bc]->GetN_DegreesOfFreedom());
      TFEFunction2D original_function(original_fe_spaces[ans][bc], "original-function",
                                  " ", &original_vals.at(0), original_vals.size());

      // transfer the first function
      original_function.Interpolate(analytic_function_1);
      trafo.transfer(original_function, into_function, into_vals);
      //...and check if it worked!
      double tol = 1e-10;
      if(ans != 0 && euclidean_distance(into_vals,reference_1_vals) > tol)
        ErrThrow("Transfer tool failed: ", ans, " ",
                 bc, " ", euclidean_distance(into_vals,reference_1_vals));


      // transfer the second function
      original_function.Interpolate(analytic_function_2);
      trafo.transfer(original_function, into_function, into_vals);
      if(ans != 0 && euclidean_distance(into_vals,reference_2_vals) > tol)
        ErrThrow("Transfer tool failed: ", ans, " ",
                 bc, " ", euclidean_distance(into_vals,reference_2_vals));

    }
  }

  delete target_space;
  for(size_t ans =0 ; ans < ansatz_orders.size() ; ++ans )
    for(size_t bc = 0; bc < bound_conds.size() ;++bc)
      delete original_fe_spaces[ans][bc];
}
