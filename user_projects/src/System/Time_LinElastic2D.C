#include <Time_LinElastic2D.h>
#include <Database.h>




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
: systems_()
{
  cout << "CONSTRUCTOR2 OK!" << endl;
}
