#include "Particles.h"
#include "Database.h"

// Drag coefficients according to 
// An investigation of particle trajectories in two-phase flow systems
// S. A. Morsi and A. J. Alexander, 1972
void get_coeffs_drag(double &b1, double &b2, double &b3, double reynold_pr)
{
  if(reynold_pr < 0.1)
  {
    b1 = 0.; b2 = 24.; b3=0.;    
  }
  else if(reynold_pr < 1.)
  {
    b1 = 3.69; b2 = 22.73; b3=0.0903;
  }
  else if(reynold_pr < 10.)
  {
    b1 = 1.222; b2 = 29.1667; b3 = -3.8889;
  }
  else if(reynold_pr < 100.)
  {
    b1 = 0.6167; b2 = 46.5; b3 = -116.67;
  }
  else if(reynold_pr < 1000.)
  {
    b1 = 0.3644; b2 = 98.33; b3 = -2778.0;
  }
  else if(reynold_pr < 5000.)
  {
    b1 = 0.357; b2 = 148.62; b3 = -4.75e4;
  }
  else if(reynold_pr < 10000.)
  {
    b1 = 0.46; b2 = -490.546; b3 = 57.87e4;
  }
  else if(reynold_pr < 50000.)
  {
    b1 = 0.5191; b2 = -1662.5; b3 = 5.4167e6;
  }
  else
  {
    b1 = b2 = b3 =0.;
  }
}

ParameterDatabase get_default_Particle_parameters()
{
  Output::print<5>("creating a default Particle parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NSE3D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("Particle parameter database");

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  // a default time database
  ParameterDatabase time_db = ParameterDatabase::default_time_database();
  db.merge(time_db,true);

  return db;
}


Particles::Particles(const ParameterDatabase& param_db)
: db_(get_default_Particle_parameters()), velocity(2)
{
  db_.merge(param_db);
  diameters.push_back(0.1);
  positions.push_back(Point(0, 2.));
}


void Particles::update_particle(const std::vector< double >& fluid_velocity, 
                                double fluid_density)
{
  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double viscos = db_["viscosity"];
  double part_slip = db_["particle_slip"];
  double current_diam = diameters.back();
  double reynold_pr = fluid_density/viscos * part_slip * current_diam;
  
  double b1, b2, b3; 
  get_coeffs_drag(b1, b2, b3, reynold_pr);
  
  double c_drag_coeff = b1 + b2/reynold_pr + b3/(reynold_pr*reynold_pr);
  
  double hdynamic_force_coeff = 3.* reynold_pr * c_drag_coeff * viscos;
  hdynamic_force_coeff /= 4. * current_diam * current_diam * density;
  
  double gz = -981.0;
  
  // preparing the right hand side 
  std::vector<double> rhs_part(2);
  for(size_t i=0; i<fluid_velocity.size(); ++i)
    rhs_part[i] = hdynamic_force_coeff * fluid_velocity[i];
  
  rhs_part[1] += gz * (1-fluid_density/density);
  
  double t2 = TDatabase::TimeDB->THETA2;
  for(size_t i=0; i< velocity.size(); ++i)
  {
    rhs_part[i] += -t2*hdynamic_force_coeff*velocity[i];
  
    rhs_part[i] *= dt;
    rhs_part[i] += velocity[i];
  }
  
  double t1 = TDatabase::TimeDB->THETA1;
  double factor = 1 + t1*hdynamic_force_coeff*dt;
  
  for(size_t i=0; i< velocity.size(); ++i)
    velocity[i] = rhs_part[i]/factor;
  
  // update the position of the particles
  std::vector<double> new_position(fluid_velocity.size());
  
  new_position[0] = positions.back().x() + dt*velocity[0];
  new_position[1] = positions.back().y() + dt*velocity[1];
  if(fluid_velocity.size() == 3)
  {
    new_position[2] = positions.back().z() + dt*velocity[2];    
  }
  positions.push_back(Point(new_position));
}


