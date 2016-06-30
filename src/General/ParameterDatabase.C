#include <ParameterDatabase.h>
#include <iostream>
#include <algorithm>
#include <ctime>
#include <memory>
#include <MooNMD_Io.h>

// helper functions
std::list<Parameter>::const_iterator 
find_parameter(std::string name, const std::list<Parameter>& parameters)
{
  auto search_lambda = [=](const Parameter& p){ return p.get_name() == name; };
  return std::find_if(parameters.begin(), parameters.end(), search_lambda);
}
std::list<Parameter>::iterator 
find_parameter(std::string name, std::list<Parameter>& parameters)
{
  auto search_lambda = [=](const Parameter& p){ return p.get_name() == name; };
  return std::find_if(parameters.begin(), parameters.end(), search_lambda);
}

bool ParameterDatabase::contains(std::string param_name) const
{
  if(find_parameter(param_name, this->parameters) == this->parameters.end())
    return false;
  else
    return true;
}

/* ************************************************************************** */
ParameterDatabase::ParameterDatabase(std::string name) 
 : name(name), parameters()
{
  
}

/* ************************************************************************** */
ParameterDatabase::ParameterDatabase(const ParameterDatabase& db)
 : name(db.name), parameters()
{
  // iterate through the parameters in db
  for(const auto& p : db.parameters)
  {
    this->parameters.emplace_back(Parameter(p));
  }
}

/* ************************************************************************** */
ParameterDatabase& ParameterDatabase::operator=(const ParameterDatabase& db)
{
  this->name = db.get_name();
  // remove parameters previously stored
  this->parameters.clear();
  // add parmeter copies from db to this
  // iterate through the parameters in db
  for(const auto& p : db.parameters)
  {
    this->parameters.emplace_back(Parameter(p));
  }
  return *this;
}

/* ************************************************************************** */
template <typename T>
void ParameterDatabase::add(std::string name, T value, std::string description)
{
  Parameter p(name, value, description);
  this->add(std::move(p));
}
template void ParameterDatabase::add(std::string n, bool v,      std::string d);
template void ParameterDatabase::add(std::string n, int v,       std::string d);
template void ParameterDatabase::add(std::string n, size_t v,    std::string d);
template void ParameterDatabase::add(std::string n, double v,    std::string d);
template void ParameterDatabase::add(std::string n,std::string v,std::string d);
template<> void ParameterDatabase::add(std::string n,const char* v,
                                       std::string d)
{
  this->add(n, std::string(v), d);
}

/* ************************************************************************** */
template <typename T>
void ParameterDatabase::add(std::string name, T value, std::string description,
                            T min, T max)
{
  Parameter p(name, value, description);
  p.set_range(min, max);
  this->add(std::move(p));
}
template void ParameterDatabase::add(std::string n, int v,    std::string d,
                                     int min, int max);
template void ParameterDatabase::add(std::string n, size_t v, std::string d,
                                     size_t min, size_t max);
template void ParameterDatabase::add(std::string n, double v, std::string d,
                                     double min, double max);

/* ************************************************************************** */
template <typename T>
void ParameterDatabase::add(std::string name, T value, std::string description, 
                            std::set<T> range)
{
  Parameter p(name, value, description);
  p.set_range(range);
  this->add(std::move(p));
}
template void ParameterDatabase::add(std::string n, bool v,      std::string d,
                                     std::set<bool> range);
template void ParameterDatabase::add(std::string n, int v,       std::string d,
                                     std::set<int> range);
template void ParameterDatabase::add(std::string n, size_t v,    std::string d,
                                     std::set<size_t> range);
template void ParameterDatabase::add(std::string n,std::string v,std::string d,
                                     std::set<std::string> range);
template<> void ParameterDatabase::add(std::string n,const char* v,
                                       std::string d, 
                                       std::set<const char*> range)
{
  std::set<std::string> range_of_strings;
  for(const auto& r : range) range_of_strings.insert(std::string(r));
  this->add(n, std::string(v), d, range_of_strings);
}

/* ************************************************************************** */
void ParameterDatabase::add(Parameter&& p)
{
  // check if such a parameter already exists in this database
  if(!this->contains(p.get_name()))
    this->parameters.emplace_back(std::move(p));
  else
    ErrThrow("parameter with this name already exists");
    // what happens to p now?
}

/* ************************************************************************** */
const Parameter& ParameterDatabase::operator[](std::string parameter_name) const
{
  auto it = find_parameter(parameter_name, this->parameters);
  if(it == this->parameters.end())
    ErrThrow("unknown parameter ", parameter_name);
  return *it;
}
Parameter& ParameterDatabase::operator[](std::string parameter_name)
{
  auto it = find_parameter(parameter_name, this->parameters);
  if(it == this->parameters.end())
    ErrThrow("unknown parameter ", parameter_name);
  return *it;
}

/* ************************************************************************** */
std::string ParameterDatabase::get_name() const
{
  return this->name;
}

/* ************************************************************************** */
void ParameterDatabase::set_name(std::string new_name)
{
  this->name = new_name;
}

/* ************************************************************************** */
size_t ParameterDatabase::get_n_parameters() const
{
  return this->parameters.size();
}

/* ************************************************************************** */
void ParameterDatabase::write(std::ostream& os, size_t verbose) const
{
  if(!os.good())
  {
    Output::print("Error in ParameterDatabase::write. stream not good");
    return;
  }
  if(verbose > 2) // make sure verbose is either 0,1, or 2
  {
    Output::print("Error in ParameterDatabase::write. unknown verbosity level");
    return;
  }
  
  if(verbose == 2)
  {
    time_t rawtime = time(NULL);
    os << "# current date and time: " << ctime(&rawtime);
    os << "# Writing a ParMooN parameter database with "
       << this->get_n_parameters() << " parameters\n";
    os << "\n";
  
    os << "# The name of the database. This is usually not of any importance\n";
    os << "[ " << this->name << " ]\n\n";
  }
  for(const auto& p : this->parameters)
  {
    if(verbose == 1 || verbose == 2) // verbose != 0
    {
      p.print_description(os, "## ", 60, "");
    }
    os << p.get_name() << ": " << p << "   " 
       << p.range_as_string() << "\n\n";
  }
}

/* ************************************************************************** */
// helper function to remove whitespace at the beginning of a string
void remove_leading_whitespace(std::string& s)
{
  // position of first non white space character
  auto pos = s.find_first_not_of(" \t");
  if(pos != std::string::npos)
    s = s.substr(pos);
  else 
    s = std::string(); // empty string
}
// helper function to remove whitespace at the end of a string
void remove_trailing_whitespace(std::string& s)
{
  // position of first non white space character
  auto pos = s.find_last_not_of(" \t");
  if(pos != std::string::npos)
    s = s.substr(0, pos+1);
  else 
    s = std::string(); // empty string
}
// helper function to create a list of all lines for one database from a stream
// additionally get the name of the database
std::list<std::string> get_lines_of_database(std::istream& is,
                                             std::string& name)
{
  std::list<std::string> lines_read; // read all lines into local variable
  std::string line; // read each line into this string
  bool found_name = false; // indicate if a database name has been found
  auto position = is.tellg(); // remember position of stream before 'getline'
  // loop over all lines until a new database is found. This is indicated by a
  // pair of brackets surrounding the name of the database 
  // '[_name_of_database_]'
  while(std::getline(is, line))
  {
    if(!line.empty() && line.at(0) == '[')
    {
      // starting position of the new name in string 'line', is right after the
      // opening bracket. The name ends at the closing bracket
      auto end_position_for_name = line.find("]");
      if(end_position_for_name != std::string::npos) // found closing bracket
      {
        if(found_name) // we found a second database here
        {
          // reset the stream position such that the next line to be read is
          // 'line', so that the second database can be read as well
          is.seekg(position);
          break; // no more lines are read into the list
        }
        found_name = true;
        name = line.substr(1, end_position_for_name-1);
        remove_leading_whitespace(name);
        remove_trailing_whitespace(name);
      }
      // else // we do nothing, maybe give a warning?
    }
    else
      lines_read.push_back(line);
    position = is.tellg(); // update to new position
  }
  if(!is)
    Output::print<3>("read last line of input stream");
  
  if(!found_name)
  {
    // no parameter database name found
    name = "custom parameter database without a given name";
  }
  return lines_read;
}

// helper function to identify the type of a variable in `value_string`. Exactly
// one of the other parameters are updated accordingly.
Parameter::types get_parameter_and_type(std::string value_string,
                                        bool & bool_ret,
                                        int & int_ret, size_t & size_t_ret, 
                                        double & double_ret,
                                        std::string & string_ret)
{
  if(value_string == "True" || value_string == "true")
  {
    bool_ret = true;
    return Parameter::types::_bool;
  }
  else if(value_string == "False" || value_string == "false")
  {
    bool_ret = false;
    return Parameter::types::_bool;
  }
  else if(value_string.find('.') != std::string::npos 
         || value_string.find(std::string("e-")) != std::string::npos)
  {
    try // this could be a double
    {
      double_ret = stod(value_string);
      return Parameter::types::_double;
    }
    catch(...)
    {
      // this could be a filename or a path
    }
  }
  else 
  {
    try
    {
      int value = stoi(value_string);
      if(value < 0) // this is an int not a size_t
      {
        int_ret = value;
        double_ret = (double)value;
        return Parameter::types::_int;
      }
      // this is an unsigned parameter
      size_t_ret = (size_t)value; // safe, because value >= 0
      // we set int_ret and double_ret as well because maybe this parameter 
      // turns out to be of those types really, then we need these numbers
      int_ret = value;
      double_ret = (double)value;
      return Parameter::types::_size_t;
    }
    catch(...)
    {}
  }
  
  // then it must be a string
  string_ret = value_string;
  return Parameter::types::_string;
}

// if the parameter range consists of values of a different type than `type`,
// that type needs to be adjusted. For example if the type is size_t but the 
// range has negative values, then the type is set to int. Also if the range 
// has a double value in it and the type is an integer type, we change the type
// to double.
void adjust_type(Parameter::types& type,
                 const std::pair<bool,std::set<std::string>>& range_list)
{
  // we need the following values to call get_parameter_and_type
  bool bool_val;
  int int_val;
  size_t size_t_val;
  double double_val;
  std::string string_val;
  
  for(const std::string & s : range_list.second)
  {
    Parameter::types t = get_parameter_and_type(s, bool_val, int_val, 
                                                size_t_val, double_val, 
                                                string_val);
    if(type == t)
      continue;
    // type and t are different
    if(type == Parameter::types::_bool || t == Parameter::types::_bool)
      ErrThrow("value and range have different types, one is bool");
    if(type == Parameter::types::_string || t == Parameter::types::_string)
      ErrThrow("value and range have different types, one is string");
    if(type == Parameter::types::_size_t && t == Parameter::types::_int)
      type = Parameter::types::_int;
    if((type == Parameter::types::_size_t || type == Parameter::types::_int)
       && t == Parameter::types::_double)
      type = Parameter::types::_double;
    // in other cases we don't change the type 
  }
}

// the bool is true in case of an interval and false otherwise
std::pair<bool,std::set<std::string>> get_range_list(std::string range)
{
  std::pair<bool, std::set<std::string>> ret; // to be returned
  if(range.empty())
    return ret; // nothing to do here
  auto npos = std::string::npos;
  // expected format: [ min, max ] or { a, b, c, d, e, ... }
  auto begin = range.find("[");
  if(begin != npos) // interval range
  {
    begin += 1;
    auto middle = range.find(",", begin);
    auto end    = range.find("]", begin);
    if(middle == npos || end == npos)
    {
      ErrThrow("wrong range format for interval.");
    }
    if(range.find(",", middle+1) < end)
    {
      // another comma between '[' and ']'
      Output::print("WARNING a parameter interval should only have two values. "
                    "Did you mean a set? Then use '{' and '}'.");
    }
    std::string min_string = range.substr(begin, middle-begin);
    std::string max_string = range.substr(middle+1, end-middle-1);
    remove_leading_whitespace(min_string);
    remove_leading_whitespace(max_string);
    remove_trailing_whitespace(min_string);
    remove_trailing_whitespace(max_string);
    ret.first = true;
    ret.second.insert(min_string);
    ret.second.insert(max_string);
  }
  else // set range
  {
    auto begin = range.find("{");
    if(begin == npos)
      return ret; // no meaningful range given
    begin += 1; // just after the "{"
    auto middle = range.find(",", begin); // find first comma
    auto end    = range.find("}", begin);
    ret.first = false;
    while(middle < end && middle != npos)
    {
      std::string value = range.substr(begin, middle-begin);
      remove_trailing_whitespace(value);
      remove_leading_whitespace(value);
      ret.second.insert(value);
      begin = middle+1; // just after the comma
      middle = range.find(",", begin); // next comma
    }
    // find last entry
    std::string value = range.substr(begin, end-begin);
    remove_trailing_whitespace(value);
    remove_leading_whitespace(value);
    ret.second.insert(value);
  }
  return ret;
}

// helper function to call set_range on the parameter with the right type
void add_range_to_parameter(Parameter& p, 
                            std::pair<bool,std::set<std::string>> range_list)
{
  if(range_list.second.empty())
    return; // nothing we can do here
  switch(p.get_type())
  {
    case Parameter::types::_bool:
    {
      if(range_list.first)
        ErrThrow("unable to set an interval range for a boolean parameter");
      auto rl = range_list.second;
      if(rl.count("true") == 1 || rl.count("True") == 1)
      {
        if(rl.count("false") == 1 || rl.count("False") == 1)
          p.set_range(std::set<bool>{true, false});
        else
          p.set_range(std::set<bool>{ true });
      }
      else if(rl.count("false") == 1 || rl.count("False") == 1)
        p.set_range(std::set<bool>{ false });
      // else not a valid range for a bool
      break;
    }
    case Parameter::types::_int:
    case Parameter::types::_size_t:
    {
      if(range_list.first)
      {
        // interval range
        std::string min_string = *range_list.second.begin();
        std::string max_string;
        if(range_list.second.size() == 1)
          max_string = min_string;
        else
          max_string = *range_list.second.rbegin();
        try
        {
          long min = stol(min_string);
          long max = stol(max_string);
          if(p.get_type() == Parameter::types::_size_t)
          {
            if(min < 0 || max < 0)
              ErrThrow("reading negative value for size_t??");
            p.set_range((size_t)min, (size_t)max);
          }
          else
            p.set_range((int)min, (int)max);
        }
        catch(...)
        {
          ErrThrow("could not read interval range for int parameter");
        }
      }
      else
      {
        // set range (not necessarily interval)
        std::set<long> long_range;
        try
        {
          for(const std::string& s : range_list.second)
          {
            long_range.insert(stol(s));
          }
        }
        catch(...)
        {
          ErrThrow("could not read range for int parameter");
        }
        // convert to the right `set` and set_range
        if(p.get_type() == Parameter::types::_size_t)
        {
          std::set<size_t> size_t_range;
          for(auto t : long_range)
          {
            if(t < 0)
              ErrThrow("reading negative value for size_t??");
            size_t_range.insert(t);
          }
          p.set_range(size_t_range);
        }
        else
        {
          std::set<int> int_range;
          for(auto t : long_range) int_range.insert(t);
          p.set_range(int_range);
        }
      }
      break;
    }
    case Parameter::types::_double:
    {
      if(!range_list.first)
        ErrThrow("unable to set a set range for a double parameter");
      std::string min_string = *range_list.second.begin();
      std::string max_string;
        if(range_list.second.size() == 1)
          max_string = min_string;
        else
          max_string = *range_list.second.rbegin();
      try
      {
        double min = stod(min_string);
        double max = stod(max_string);
        p.set_range(min, max);
      }
      catch(...)
      {
        ErrThrow("could not read range for double parameter");
      }
      break;
    }
    case Parameter::types::_string:
    {
      if(range_list.first)
        ErrThrow("unable to set an interval range for a string parameter");
      p.set_range(range_list.second);
      break;
    }
    default:
      ErrThrow("unknown parameter type");
      break;
  }
}

// helper function to create a new parameter where the type is to be determined
Parameter create_parameter(std::string name, std::string value_string, 
                           std::string description, std::string range_string)
{
  bool bool_val;
  int int_val;
  size_t size_t_val;
  double double_val;
  std::string string_val;
  
  // p will be returned at the end. We use a pointer because the Parameter class
  // has no default constructor.
  std::unique_ptr<Parameter> p;
  // find out what type this is
  Parameter::types type = get_parameter_and_type(value_string, bool_val,
                                                 int_val, size_t_val, 
                                                 double_val, string_val);
  
  auto range_list = get_range_list(range_string);
  
  adjust_type(type, range_list);
  
  // create a new parameter according to `type`
  switch(type)
  {
    case Parameter::types::_bool:
      p.reset(new Parameter(name, bool_val, description));
      break;
    case Parameter::types::_int:
      p.reset(new Parameter(name, int_val, description));
      break;
    case Parameter::types::_size_t:
      p.reset(new Parameter(name, size_t_val, description));
      break;
    case Parameter::types::_double:
      p.reset(new Parameter(name, double_val, description));
      break;
    case Parameter::types::_string:
      p.reset(new Parameter(name, string_val, description));
      break;
    default:
      ErrThrow("unknown type");
      break;
  }
  add_range_to_parameter(*p, range_list);
  return std::move(*p);
}

// helper function to read a single parameter from one line and return strings 
// for name, value and range
bool read_parameter(const std::string& line, std::string& name, 
                    std::string& value, std::string& range_string)
{
  auto npos = std::string::npos;
  // find ':'. Whatever is before that is the parameter name
  auto position_for_name = line.find(":");
  if(position_for_name == npos)
    return false;
  
  name = line.substr(0, position_for_name);
  remove_leading_whitespace(name);
  remove_trailing_whitespace(name);
  
  // whatever is after ':' is the value
  position_for_name += 1;
  // what remains is the value and possibly the range
  std::string remainder = line.substr(position_for_name);
  // remove leading whitespace
  remove_leading_whitespace(remainder);
  
  if(remainder.size() == 0)
    return false; // no value given, maybe write a warning
  
  auto position_for_value = remainder.find(" "); // may be npos
  value = remainder.substr(0, position_for_value);
  
  range_string.clear();
  if(position_for_value != npos)
  {
    // the now remaining part is the range which is not always given
    range_string = remainder.substr(position_for_value);
    remove_leading_whitespace(range_string);
    remove_trailing_whitespace(range_string);
  }
  return true;
}

/* ************************************************************************** */
void ParameterDatabase::read(std::istream& is)
{
  if(!is.good())
  {
    Output::print("Error in ParameterDatabase::read. stream not good");
    return;
  }
  
  ParameterDatabase tmp("");
  
  Output::print<2>("\nReading database from stream");
  
  // get all lines of the input stream as a list of strings
  std::list<std::string> lines_read = get_lines_of_database(is, tmp.name);
  
  std::string description; // accumulates from several lines
  // loop over all lines
  for(std::string & line : lines_read)
  {
    remove_leading_whitespace(line);
    // check if this line is empty
    if(line.length() == 0)
    {
      // no parameter on this line, reset the description
      description.clear();
      continue; // go to next line
    }
    if(line.at(0) == '#')
    {
      // documentation for a parameter on this line
      std::string des = line.substr(1);
      remove_leading_whitespace(des);
      // if there were two '#' at the beginning, then interpret this as a 
      // comment which should not become part of the description
      if(des.length() != 0 && des.at(0) == '#')
        description.clear();
      else
        // add this line to the descripion
        description += des + " ";
      continue; // go to next line
    }
    
    // check for a parameter
    std::string param_name, value_string, range_string;
    bool found_parameter = read_parameter(line, param_name, value_string,
                                          range_string);
    Output::print<4>("Reading parameter from stream: ", param_name, " \t ",
                     value_string, " \t ->", range_string, "<-     ",
                     found_parameter);
    if(!found_parameter)
    {
      description.clear();
      continue;
    }
    
    tmp.add(create_parameter(param_name, value_string, description,
                             range_string));
    description.clear();
  }
  Output::print<2>("Done reading database from stream. Read ",
                   this->parameters.size(), " parameters");
  
  this->name = tmp.get_name();
  this->merge(tmp);
}

/* ************************************************************************** */
void ParameterDatabase::merge(const ParameterDatabase &other,
                              bool create_new_parameters)
{
  for(const Parameter & p : other.parameters)
  {
    if(this->contains(p.get_name()))
    {
      this->operator[](p.get_name()).impose(p);
    }
    else if(create_new_parameters)
      this->add(Parameter(p)); // add a copy of the parameter p
  }
}

/* ************************************************************************** */
void ParameterDatabase::info(bool only_names) const
{
  Output::print("Parameter database: ", this->name);
  Output::print("  number of parameters: ", this->parameters.size());
  for(const auto& p : this->parameters)
  {
    if(only_names)
      Output::print("    ", p.get_name(), ": ", p.value_as_string());
    else
      p.info();
  }
}

/* ************************************************************************** */
// set default parameters
ParameterDatabase ParameterDatabase::parmoon_default_database()
{
  ParameterDatabase db("default ParMooN parameter database");
  
  // add parameters which are needed by all ParMooN programs and don't belong
  // anywhere else.
  
  db.add("outfile", "default_parmoon_outfile.out",
         "This is the file where all output of ParMooN is (usually) written "
         "to. In general ParMooN produces text output on the console as well "
         "as in this file. For this to properly work, you should call "
         "`Output::set_outfile(db[\"outfile\"]);` in your main program.");

  db.add("mesh_file", "__nofile__",
         "This files describes the computational mesh in .mesh format. "
         "Set this to the path of your desired mesh file.");

  db.add("problem_type", (size_t)0, 
         "Determine which kind of problem you want to solve. A value of 0 "
         "means not set. Other values have the following meanings: "
         "1: stationary convection-diffusion,  2: time-dependent "
         "convection-diffusion,  3: stationary Stokes,  4: time-dependent "
         "Stokes,  5: stationary Navier-Stokes,  6: time-dependent "
         "Navier-Stokes.",
         (size_t)0, (size_t)6);

  db.add("output_write_ps", false,
	 "Draw a postscript file of the domain. This only works in two space "
	 "dimensions. Usually this is used in the main program.",
	 {true,false});

  db.add("verbosity", (size_t)1,
         "Set the verbosity of ParMooN. The higher the number, the more will "
         "output you will get. Such output will be written to console and the "
         "'outfile'.", (size_t)1, (size_t)5);
  
  db.add("example", 0,
         "Choose which example to run. \nNote that depending on the type of "
         "problem you want to solve, different values are meaningful here. See "
         "the class 'Example' and its derived classes.", -5, 200);

  return db;
}

/* ************************************************************************** */
// set default parameters for time discretization
ParameterDatabase ParameterDatabase::default_time_database()
{
  ParameterDatabase db("default ParMooN time parameters database");
  
  db.add("time_start", 0.,
         "This is the start time of the simulation. The total simulated time "
         "interval is [start_time, end_time].", -1000.,1000.);
  
  db.add("time_end", 1.,
         "This is the end time of the simulation. The total simulated time "
         "interval is [start_time, end_time].", -1000.,1000.);
  
  db.add("time_step_length", 0.05,
         "This is the time step length. Without time adaptivity, it "
         "is the same constant value for all time steps, whereas for "
         "time adaptivity, it only corresponds to the initial value.",
         0., 0.5);
  
  db.add("time_discretization", (size_t)2,
         "This is the time discretization scheme. The following values are "
         "implemented :"
         "0 -> Forward Euler, "
         "1 -> Backward Euler, "
         "2 -> Crank-Nicholson, "
         "3 -> Fractional step."
         "4 -> Extrapolated Crank_Nicholson (or IMplicit-EXplicit, IMEX)",
         (size_t)0 , (size_t)4 );
  
  return db;
}

ParameterDatabase ParameterDatabase::default_nonlinit_database()
{
  ParameterDatabase db("default ParMooN nonlinear iteration parameters database");

  //TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE, TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR;
  db.add("nonlinloop_maxit", (size_t) 100,
         "The maximum number of iterations to perform in a non-linear loop.",
         (size_t) 0, size_t (1000));

  // TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE, TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR,
  db.add("nonlinloop_epsilon", 1e-10,
         "At which absolute residual to break the nonlinear loop.",
         0. , 1. );
  
  db.add("nonlinloop_damping_factor", 1.0,
         "Damping factor 'w' for the nonlinear iteration. The solution of the "
         "k-th iterate will be scaled by 'w'. Then The previous solution, "
         "scaled by '1-w', will be added. Setting to it to zero makes no "
         "sense.", 
         0., 1.);

  //TDatabase::ParamDB->SC_NONLIN_DIV_FACTOR
  db.add("nonlinloop_slowfactor", (double) 1e10,
         "Determines at which reduction rate over x iterations"
         "(usually x = 10, see system classes) a convergence is interpreted"
         "as too slow and therefore the iteration is stopped.",
         0., (double) 1e10);

  //TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE, TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SCALAR
  db.add("nonlinloop_scale_epsilon_with_size", false,
         "Whether or not to scale the absolute residual breaking criterion"
         "with the square root of the problem size.",
         {true,false});

  return db;
}

ParameterDatabase ParameterDatabase::default_output_database()
{
	  ParameterDatabase db("default ParMooN output control parameters database");

	  db.add("output_write_vtk", false,
			  "This parameter can control, whether an output method"
			  "of a system class will produce VTK output or not.",
			  {true,false});

	  db.add("output_write_case", false,
			  "This parameter can control, whether an output method"
			  "of a system class will produce CASE output or not.",
			  {true,false});

	  db.add("output_compute_errors", true,
	         "Do or do not compute errors after computing a solution. This makes "
	         "much sense if an analytical solution is known. If not then it is "
	         "often simply set to zero and computing errors then means computing "
	         "norms, e.g. the L^2-norm of the solution.",
			 {true,false});

	  db.add("output_vtk_directory", ".",
	         "This directory is where the VTK output is written. This "
	         "directory will be created, if it does not exist already. Files in "
	         "this directory will be overwritten without any warning.");

	  db.add("output_basename", "parmoon",
	         "This string is prepended to most files written by ParMooN. "
	         "This includes also vtk- and case-files");

	  db.add("steps_per_output", 1,
	         "This integer specifies how many (time) steps are performed "
		 "before writing the results ");

	  return db;
}

ParameterDatabase ParameterDatabase::default_tetgen_database()
{
  ParameterDatabase db("default ParMooN mesh generation using TetGen "
                        "parameters database");

  // maximum can be infnity
  db.add("tetgen_quality", 1.4, " This value is for the aspect ratio", 
         1.0, 1000.);
  
  db.add("tetgen_steiner", (size_t) 0,
         "this parameter 'preserve the mesh on the exterior boundary' ", 
         (size_t) 0, (size_t) 1);
  
  // maximum can be anything depending on geometry
  db.add("tetgen_volume", 1.0, "this parameter is for the maximum volume "
          "that depends on the geometry", 0.0, 1000.);
  
  db.add("tetgen_merge_colplaner", (size_t) 0, 
         "This parameter is for the coplanar facets to be merge "
         "or very close vertices", (size_t) 0, (size_t) 1);
  
  db.add("tetgen_quiet", (size_t) 1, "Quiet: No terminal output except errors ",
         (size_t) 0, (size_t) 1);
  return db;
}
