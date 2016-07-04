#include <ParameterDatabase.h>
#include <MooNMD_Io.h>

bool test_for_parameter_class()
{
  Parameter p("parameter", 1, "dummy description");
  Parameter q("unsigned_parameter", (size_t)1, "dummy description");
  try
  {
    p.impose(q);
    Output::print("I was able to impose one parameter to another one, even "
                  "though their names are different");
    return false;
  }
  catch(...) { }
  
  Parameter p2("parameter", 2, "some description");
  p.set_range(-2, 4);
  p2.set_range<int>({-2, 2, 5});
  // impose a parameter with the same type but different type of range
  p2.impose(p);
  p2.info();
  
  Parameter p3("parameter", (size_t)6, "some description");
  // impose a parameter with a different type
  p.impose(p3);
  p.info();
  
  Parameter p4("parameter", 2.0, "double parameter description");
  p4.impose(p3);
  p4.info();
  p4.impose(p);
  p4.info();
  return true;
}

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
  
  // test piping a Parameter into a stream
  os << p;
  
  ParameterDatabase timedb = ParameterDatabase::default_time_database();

  if(!test_for_parameter_class())
    return 1;
  
  Output::print("successful test");
  return 0;
}
