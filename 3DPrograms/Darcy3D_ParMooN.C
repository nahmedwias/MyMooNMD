// =======================================================================
//
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <Darcy.h>

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  double t_start = GetTime();
  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase3D FEDatabase;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.read(std::string(argv[1]));
  
  //open OUTFILE, this is where all output is written to (addionally to console)
  Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain domain(parmoon_db, argv[1]);
 
  // write all Parameters to the OUTFILE (not to console) for later reference
  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  
  // refine grid
  size_t n_ref = domain.get_n_initial_refinement_steps();
  for(size_t i = 0; i < n_ref; i++)
    domain.RegRefineAll();
  
  // choose example according to the value of db["example"]
  Example_Darcy3D example(parmoon_db);
  
  //=========================================================================
  Darcy<3> darcy3d(domain, parmoon_db, example);
  darcy3d.assemble();
  darcy3d.solve();
//   double a = 1./3.;
//   darcy3d.get_solution()[12] = -a;
//   darcy3d.get_solution()[13] = a;
//   darcy3d.get_solution()[14] = a;
//   darcy3d.get_solution()[15]= -a;
// 
//   darcy3d.get_solution()[24] = a;
//   darcy3d.get_solution()[25] = -a;
//   darcy3d.get_solution()[26] = -a;
//   darcy3d.get_solution()[27] = a;
//   auto& space = darcy3d.get_velocity_space();
//   double v[3];
//   darcy3d.get_velocity().FindValueLocal(space.GetCollection()->GetCell(0), 
//                                         0, 1., 0.5, 0.5, v);
//   Output::print("values: ", v[0], ", ", v[1], ", ", v[2]);
//   darcy3d.get_solution().print("sol");
  darcy3d.output();
  //=========================================================================
  
  Output::print<1>("used time: ", GetTime() - t_start, " seconds");
  Output::close_file();
  return 0;
} // end main
