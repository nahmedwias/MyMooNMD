/**
 * @file New declaration of a multigrid object, which holds the necessary
 * grid information for executing a multigrid iteration.
 *
 * TODO Can't do anything so far, this is only a dummy.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#include <BlockVector.h>
#include <GridTransfer.h>
#include <Multigrid.h>
#include <MooNMD_Io.h>
#include <ParameterDatabase.h>


Multigrid::Multigrid(const ParameterDatabase& db,
                     std::list<BlockFEMatrix*> matrices)
{
  for(auto mat : matrices)
  {
    double beta = 1; //TODO from db
    SmootherCode sm = SmootherCode::DIRECT_SOLVE; //TODO from db

    levels_.push_back(MultigridLevel(mat, beta, sm));
  }

  size_t n_levels = levels_.size();

  for(size_t i =0; i<n_levels ;++i)
    damp_factors_.push_back(1); //TODO from db

  cycle_ = MGCycle::F; //TODO from db

  set_cycle_control();

  //CB DEBUG
  print_cycle_control();
  Output::print("Dummy Multigrid object constructed.");
  //END DEBUG
}

void Multigrid::set_finest_sol(const BlockVector& bv)
{
  //TODO Eventually here a call to "RestrictToAllGrids" is necessary
  // - so that for step length control or Braess-Sarazin an initial
  // residual can be calculated on each level
  ;
}

void Multigrid::set_finest_rhs(const BlockVector& bv)
{
  ;
}

void Multigrid::cycle()
{
  ;
}

const BlockVector& Multigrid::get_finest_sol()
{
  ;
}

void Multigrid::set_cycle_control()
{
  size_t n_levels = levels_.size();

  std::vector<int> mg_recursions(n_levels);
  // coarsest grid
  if (cycle_ == MGCycle::V)
  {
    std::fill(mg_recursions.begin(), mg_recursions.end(), 1);
  }
  else if (cycle_ == MGCycle::W)
  {
    std::fill(mg_recursions.begin(), mg_recursions.end(), 2);
  }
  else if (cycle_ == MGCycle::F)
  {
    std::fill(mg_recursions.begin(), mg_recursions.end(), 2);
  }
  mg_recursions[n_levels-1] = 1;

  fill_recursively(mg_recursions, n_levels-1);

}

void Multigrid::fill_recursively(std::vector<int>& mg_recursions, int level)
{
  Output::print("Level ", level);

  int coarsest_level = 0;
  int finest_level = levels_.size() - 1;

  if(level == coarsest_level)
  {
    cycle_control_.push_back(MGDirection::Up);
    return;
  }

  for(int j=0;j<mg_recursions[level];j++)
  {
    cycle_control_.push_back(MGDirection::Down);
    fill_recursively(mg_recursions, level-1);
  }

  if (cycle_ == MGCycle::F)
    mg_recursions[level] = 1;

  if(level == finest_level) //finest level
  {
    cycle_control_.push_back(MGDirection::End);
  }
  else //other levels
  {
    cycle_control_.push_back(MGDirection::Up);
  }

}

void Multigrid::print_cycle_control() const
{
  Output::print("Look at that cycle_control_ !");
  for(auto cc : cycle_control_)
  {
    if(cc == MGDirection::Up)
      Output::print("Up");
    if(cc == MGDirection::Down)
      Output::print("Down");
    if(cc == MGDirection::End)
      Output::print("End");
  }
}


ParameterDatabase Multigrid::default_multigrid_database()
{
  Output::print<3>("creating a default multigrid parameter database");
  ParameterDatabase db("default multigrid database");

  db.add<size_t>("multigrid_n_levels", 2,
         "Determine how many levels the multigrid cycle consists of.", 0, 5);

  db.add("multigrid_cycle_type", std::string("V"),
         "The recursion type how to traverse the multigrid levels."
         "So far the three standard cycle V, W and F are implemented."
         , {"V", "W", "F"});

  db.add("multigrid_smoother", std::string("direct_solve"),
         "The smoother to use on all but the coarsest level."
         "You should take care, that the smoother you chose fits"
         "your problem type, e.g. Vanka smoothers are best fitted"
         "for saddle point problems.",
         {"direct_solve", "nodal_vanka", "cell_vanka", "batch_vanka"});

  db.add("multigrid_smoother_coarse", std::string("direct_solve"),
         "The smoother to use on the coarsest level."
         "You should take care, that the smoother you chose fits"
         "your problem type, e.g. Vanka smoothers are best fitted"
         "for saddle point problems.",
         {"direct_solve", "nodal_vanka", "cell_vanka", "batch_vanka"});

  db.add("multigrid_correction_damp_factor", 1.0,
         "The damping factor which is used when applying the coarse"
         "grid correction to a finer grid. A factor of 1.0 means:"
         "no damping, a factor of 0.0: nothing changes.",
         0.0,1.0);

  db.add<size_t>("multigrid_n_pre_smooth", 1,
                 "The number of smoothing steps to apply per level"
                 "before going down to the next coarsest level.",
                 1, 10);

  db.add<size_t>("multigrid_n_post_smooth", 1,
                 "The number of smoothing steps to apply per level"
                 "after coming up from the next coarsest level.",
                 1, 10);

  db.add("multigrid_coarse_residual", 1.0e-1,
         "The target residual on the coarsest grid."
         "When this residual is reached on the coarsest grid by solving"
         "or smoothing, the coarse level will return and the process"
         "continue on the next finest level.",
         1.0e-20, 1.0);


  return db;
}














