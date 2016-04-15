#include <ParameterDatabase.h>
#include <MooNMD_Io.h>

int main(int argc, char* argv[])
{
  Output::print("starting parameter test");
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  ParameterDatabase db2 = ParameterDatabase::parmoon_default_database();
  db.set_name("some test database");
  
  db.add("double_parameter", 1.2, "a dummy parameter", -2.0, 2.0);
  db.add("int_parameter", 0, "test parameter", {-4, -2, 0, 2, 4});
  
  try
  {
    // add a parameter with a name which is already in the database
    // this should throw an exception
    db.add("int_parameter", 1, "wrong parameter");
    db.info();
    Output::print("I was able to add two parameter with the same name");
    db.info();
    return 1;
  }
  catch(...)
  {
    // fine, caught the expected exception
  }
  
  // get a parameter
  Parameter & p = db["double_parameter"];
  p.set_range(-3.0, 12.64789);
  p = 10.0;
  if(!p.is(10.0))
  {
    Output::print("parameter is not 10.0 as it should be");
    return 1;
  }
  
  if(db.contains("not_existing_parameter"))
  {
    Output::print("ParameterDatabase::contains gives true where it shouldn't");
    return 1;
  }
  
  std::stringstream os;
  db.write(os);
  db2.read(os);
  if(db.get_n_parameters() != db2.get_n_parameters())
  {
    Output::print("ParameterDatabase write->read does not work");
    return 1;
  }
  
  Output::print("successful test");
  return 0;
}
