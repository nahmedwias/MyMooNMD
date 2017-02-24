#include "RungeKuttaTable.h"
#include <MooNMD_Io.h>
#include <Database.h>

unsigned int get_n_stages_from_database(std::string time_disc)
{
  if(time_disc.compare("crank_nicolson") == 0)
    return 2;
  else
    return 1;
}

RungeKuttaTable::RungeKuttaTable(std::string time_disc)
: n_stages(get_n_stages_from_database(time_disc)), order_p(0), order_q(0), type(0),
   a(n_stages, std::vector<double>(n_stages, 0.)), b(n_stages, 0.), 
   c(n_stages, 0.), bh(n_stages, 0.)
{
  if(time_disc.compare("forward_euler") == 0)
  {
    type = 0;
    order_p = 1;
    c[0] = 0.; 
    a[0][0] = 0.;
    b[0] = 1.;
  }
  else if(time_disc.compare("backward_euler") == 0)
  {
    type = 1;
    // stability function 1/(1-z)
    order_p = 1;
    c[0] = 1.;
    a[0][0] = 1.;
    b[0] = 1.;
    
    // for internal usage set the variable in the database
    TDatabase::TimeDB->THETA1 = 1.;
    TDatabase::TimeDB->THETA2 = 0.;
    TDatabase::TimeDB->THETA3 = 0.;
    TDatabase::TimeDB->THETA4 = 1.;
  }
  else if(time_disc.compare("crank_nicolson") == 0)
  {
    type = 1;
    // stability function (1+z/2) / (1-z/2)
    order_p = 2;
    c[0] = 0;  /*|*/ a[0][0]=0.; a[0][1]=0.;
    c[1] = 1.; /*|*/ a[1][0]=.5; a[1][1]=.5;
               /*|*/ b[0]=.5;    b[1]=.5;
               
    TDatabase::TimeDB->THETA1 = 0.5;
    TDatabase::TimeDB->THETA2 = 0.5;
    TDatabase::TimeDB->THETA3 = 0.5;
    TDatabase::TimeDB->THETA4 = 0.5;
  }
  else
  {
    ErrThrow("time_disc = ", time_disc , " is not implemented yet" );
  }
}
