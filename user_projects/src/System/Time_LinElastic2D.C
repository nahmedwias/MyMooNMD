#include <Time_LinElastic2D.h>
#include <Database.h>
#include <DirectSolver.h>
#include <MainUtilities.h>


/**************************************************************************** */
ParameterDatabase get_default_Time_LinElastic2D_parameters()
{
  Output::print<5>("creating a default Time_LinElastic2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TLinElastic2D parameter database");

  //set up a nonlinit_database and merge
  ParameterDatabase nl_db = ParameterDatabase::default_nonlinit_database();
  db.merge(nl_db,true);

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  // a default time database
  ParameterDatabase time_db = ParameterDatabase::default_time_database();
  db.merge(time_db,true);

  return db;
}


/**************************************************************************** */
Time_LinElastic2D::System_per_grid::System_per_grid(const Example_TimeLinElastic2D& example,
                                                    TCollection& coll)
: fe_space_(&coll, (char*)"space", (char*)"time_linelastic2d space", example.get_bc(0),
            TDatabase::ParamDB->ANSATZ_ORDER, nullptr),
  stiffness_matrix_({&fe_space_,&fe_space_}),
  mass_matrix_({&fe_space_,&fe_space_}),
  rhs_(this->stiffness_matrix_,true),
  solution_(this->stiffness_matrix_,false),
  u_(&fe_space_,(char*)"u", (char*)"u", solution_.block(0),
     solution_.length(0), 2)
{
  // K has 4 blocks, M has 2 diagonal blocks
//  stiffness_matrix_ = BlockFEMatrix::LinElastic2D(fe_space_);
//  mass_matrix_ = BlockFEMatrix::Mass_LinElastic2D(fe_space_);

  cout << "CONSTRUCTOR OF SYSTEM PER GRID OK!!" << endl;
}


/**************************************************************************** */
Time_LinElastic2D::Time_LinElastic2D(const TDomain& domain,
                                     const ParameterDatabase& param_db,
                                     int reference_id)
: Time_LinElastic2D(domain, param_db, Example_TimeLinElastic2D(param_db), reference_id)
{
  cout << "CONSTRUCTOR1 OK!" << endl;
}


/**************************************************************************** */
Time_LinElastic2D::Time_LinElastic2D(const TDomain& domain,
                                     const ParameterDatabase& param_db,
                                     const Example_TimeLinElastic2D& ex,
                                     int reference_id)
: db_(get_default_Time_LinElastic2D_parameters()),
  example_(ex), solver_(param_db), outputwriter_(param_db),
  systems_()
{
  db_.merge(param_db,false);

  /* Construction of the systems per Grid */
  bool usingMultigrid = this->solver_.is_using_multigrid();
  if(!usingMultigrid)
  {
    // create the collection of cells from the domain (finest grid)
    TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
    // create finite element space and function, a matrix, rhs, and solution
    this->systems_.emplace_back(this->example_, *coll);

  }
  else
  {
    ErrThrow("MULTIGRID NOT IMPLEMENTED YET!!!");
  }


  cout << "CONSTRUCTOR2 OK!" << endl;
}
