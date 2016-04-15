#include <Parameter.h>
#include <iostream>
#include <algorithm>
#include <iomanip>      // std::setw
#include <sstream>
#include <iterator>
#include <MooNMD_Io.h>

/// @brief Check that name is allowed (no spaces or dots), it returns the name
/// without any changes, or throws
std::string check_name(std::string name)
{
  // Space and punctuation not allowed in names
  for (std::size_t i = 0; i < name.size(); i++)
  {
    if (name[i] == ' ' || name[i] == '.' || name[i] == ':')
    {
      ErrThrow("illegal name for a Parameter. No spaces or dots allowed.");
    }
  }
  return name;
}

/// @brief check if the given Parameter type is the corresponding template type
template <typename T> bool check_type(Parameter::types t) { return false; }
template<> bool check_type<bool>(Parameter::types t)
{ return Parameter::types::_bool == t; }
template<> bool check_type<int>(Parameter::types t)
{ return Parameter::types::_int == t; }
template<> bool check_type<size_t>(Parameter::types t)
{ return Parameter::types::_size_t == t; }
template<> bool check_type<double>(Parameter::types t)
{ return Parameter::types::_double == t; }
template<> bool check_type<std::string>(Parameter::types t)
{ return Parameter::types::_string == t; }

template <typename T>
std::string convert_to_string(T x)
{
  // std::to_string does not work for std::string and more importantly for 
  // T=double it truncates to 6 digits, which is not what we want here
  std::ostringstream os;
  os << std::setprecision(std::numeric_limits< double >::digits10+1) << x;
  return os.str();
}

template <typename T>
bool check_value_in_range(T new_value, const std::set<T>& range)
{
  return (std::find(range.begin(), range.end(), new_value) != range.end());
}

std::string type_as_string(const Parameter::types t)
{
  switch(t)
  {
    case Parameter::types::_bool:   return std::string("bool");   break;
    case Parameter::types::_int:    return std::string("int");    break;
    case Parameter::types::_size_t: return std::string("size_t"); break;
    case Parameter::types::_double: return std::string("double"); break;
    case Parameter::types::_string: return std::string("string"); break;
    default: return std::string("unknown Parameter type"); break;
  }
}


/* ************************************************************************** */
Parameter::Parameter(std::string key, bool new_value, std::string description)
 : type(types::_bool), bool_value(new_value), int_value(), unsigned_value(), 
   double_value(), string_value(), access_count(0), change_count(0),
   name(check_name(key)), description(description), 
   not_bool_value_allowed(true), int_range(), unsigned_range(), min(), max(),
   string_range()
{
}

//template <typename T>
Parameter::Parameter(std::string key, int new_value, std::string description)
 : type(types::_int), bool_value(), int_value(new_value), unsigned_value(), 
   double_value(), string_value(), access_count(0), change_count(0),
   name(check_name(key)), description(description), not_bool_value_allowed(),
   int_range({new_value}), unsigned_range(), min(), max(), string_range()
{
}

Parameter::Parameter(std::string key, size_t new_value, std::string description)
 : type(types::_size_t), bool_value(), int_value(), unsigned_value(new_value), 
   double_value(), string_value(), access_count(0), change_count(0),
   name(check_name(key)), description(description), not_bool_value_allowed(),
   int_range(), unsigned_range({new_value}), min(), max(), string_range()
{
}

Parameter::Parameter(std::string key, double new_value, std::string description)
 : type(types::_double), bool_value(), int_value(), unsigned_value(), 
   double_value(new_value), string_value(), access_count(0), change_count(0),
   name(check_name(key)), description(description), not_bool_value_allowed(),
   int_range(), unsigned_range(), min(new_value), max(new_value), string_range()
{
}

Parameter::Parameter(std::string key, std::string new_value, std::string description)
 : type(types::_string), bool_value(), int_value(), unsigned_value(), 
   double_value(), string_value(new_value), access_count(0),
   change_count(0), name(check_name(key)), description(description),
   not_bool_value_allowed(), int_range(), unsigned_range(), min(), max(),
   string_range({new_value})
{
}

Parameter::Parameter(const Parameter& p)
 : type(p.type), bool_value(p.bool_value), int_value(p.int_value),
   unsigned_value(p.unsigned_value), double_value(p.double_value),
   string_value(p.string_value), access_count(0), change_count(0), name(p.name),
   description(p.description), not_bool_value_allowed(p.not_bool_value_allowed),
   int_range(p.int_range), unsigned_range(p.unsigned_range), min(p.min),
   max(p.max), string_range(p.string_range)
{
  Output::print<3>("Parameter(const Parameter& p)");
}


/* ************************************************************************** */
void Parameter::impose(const Parameter& p)
{
  if(this->name != p.get_name())
  {
    ErrThrow("cannot impose another parameter to this one because it has a "
             "different name: ", this->name, " ", p.get_name());
  }
  if(this->type != p.get_type())
  {
    ErrThrow("cannot impose another parameter to this one because it has a "
             "different type: ", type_as_string(this->type), " ", 
             type_as_string(p.get_type()));
  }
  // reset change_count and access_count
  this->access_count = 0;
  this->change_count = 0;
  // handle the description
  if(this->description.empty())
  {
    this->description = p.get_description();
  }
  else if(!p.get_description().empty()
          && this->description != p.get_description())
  {
    // concatenate the descriptions
    this->description += std::string(". ") + p.get_description();
  }
  // else: we just leave the description as it is, because the description of p
  // is either empty or equal.
  
  // handle the range (it will be at least as large as that of p)
  this->not_bool_value_allowed |= p.not_bool_value_allowed; // logical or
  this->int_range.insert(p.int_range.begin(), p.int_range.end());
  this->unsigned_range.insert(p.unsigned_range.begin(), p.unsigned_range.end());
  this->min = std::min(this->min, p.min);
  this->max = std::max(this->max, p.max);
  this->string_range.insert(p.string_range.begin(), p.string_range.end());
  
  // handle the value
  // the value is included in the range, because it was included in the range of
  // p and this->range is a superset of p.range. So no additional check needed 
  // here.
  this->bool_value = p.bool_value;
  this->int_value = p.int_value;
  this->unsigned_value = p.unsigned_value;
  this->double_value = p.double_value;
  this->string_value = p.string_value;
}

/* ************************************************************************** */
std::string Parameter::get_name() const
{
  return name;
}

/* ************************************************************************** */
std::string Parameter::get_description() const
{
  return description;
}

/* ************************************************************************** */
std::size_t Parameter::get_access_count() const
{
  return access_count;
}

/* ************************************************************************** */
std::size_t Parameter::get_change_count() const
{
  return change_count;
}

/* ************************************************************************** */
std::string Parameter::range_as_string() const
{
  std::string ret; // this will be returned
  switch(this->type)
  {
    case types::_bool:
      if(this->not_bool_value_allowed) // range consists of true and false
        return std::string("{ true, false }");
      else if(this->bool_value) // range consists only of the current value
        return std::string("{ true }");
      else
        return std::string("{ false }");
      break;
    case types::_int: // almost equal to the size_t case
    {
      int n_param = this->int_range.size();
      auto minmax = std::minmax_element(this->int_range.begin(),
                                        this->int_range.end());
      bool range_is_interval = (n_param == (*minmax.second - *minmax.first + 1));
      if(range_is_interval)
      {
        ret += "[ " + std::to_string(*minmax.first) + ", ";
        ret += std::to_string(*minmax.second) + " ]";
      }
      else
      {
        auto it = this->int_range.begin();
        ret += "{ " + std::to_string(*it);
        ++it; // it points to second entry in this->int_range
        for(; it != this->int_range.end(); ++it)
        {
          ret += ", " + std::to_string(*it);
        }
        ret += " }";
      }
      break;
    }
    case types::_size_t: // almost equal to the int case
    {
      size_t n_param = this->unsigned_range.size();
      auto minmax = std::minmax_element(this->unsigned_range.begin(),
                                        this->unsigned_range.end());
      bool range_is_interval = (n_param == (*minmax.second - *minmax.first + 1));
      if(range_is_interval)
      {
        ret += "[ " + std::to_string(*minmax.first) + ", ";
        ret += std::to_string(*minmax.second) + " ]";
      }
      else
      {
        auto it = this->unsigned_range.begin();
        ret += "{ " + std::to_string(*it);
        ++it; // it points to second entry in this->unsigned_range
        for(; it != this->unsigned_range.end(); ++it)
        {
          ret += ", " + std::to_string(*it);
        }
        ret += " }";
      }
      break;
    }
    case types::_double:
      ret += "[ " + convert_to_string(this->min) + ", ";
      ret += convert_to_string(this->max) + " ]";
      break;
    case types::_string:
    {
      auto it = this->string_range.begin();
      ret += "{ " + *it;
      ++it; // it points to second entry in this->string_range
      for(; it != this->string_range.end(); ++it)
      {
        ret += ", " + *it;
      }
      ret += " }";
      break;
    }
    default:
      ErrThrow("unknown parameter type");
      break;
  }
  return ret;
}

/* ************************************************************************** */
Parameter::types Parameter::get_type() const
{
  return this->type;
}

/* ************************************************************************** */
template<>
bool Parameter::get<bool>() const
{
  if(!check_type<bool>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != bool");
  access_count++;
  return this->bool_value;
}
template<>
int Parameter::get<int>() const
{
  if(!check_type<int>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != int");
  access_count++;
  return this->int_value;
}
template<>
size_t Parameter::get<size_t>() const
{
  if(!check_type<size_t>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != size_t");
  access_count++;
  return this->unsigned_value;
}
template<>
double Parameter::get<double>() const
{
  if(!check_type<double>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != double");
  access_count++;
  return this->double_value;
}
template<>
std::string Parameter::get<std::string>() const
{
  if(!check_type<std::string>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != std::string");
  access_count++;
  return this->string_value;
}

/* ************************************************************************** */
template <typename T>
bool Parameter::is(T value) const
{
  if(!check_type<T>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != typename T");
  return this->get<T>() == value;
}
template bool Parameter::is(bool) const;
template<> bool Parameter::is(int value) const
{
  // we implement this explicitly, because int is not implicitly converted to 
  // size_t, so Parameter::is<int>(int) gives wrong type errors if the type of
  // this parameter is size_t.
  if(check_type<int>(this->type))
    return this->int_value == value;
  else if(check_type<size_t>(this->type) && value >= 0)
    return this->unsigned_value == (size_t)value;
  else
    ErrThrow("wrong type: ", type_as_string(this->type), " != int or size_t");
}
template bool Parameter::is(size_t) const;
template bool Parameter::is(double) const;
template bool Parameter::is(std::string) const;
template<> bool Parameter::is(const char* v) const
{ return this->is<std::string>(std::string(v)); }


/* ************************************************************************** */
Parameter::operator bool() const
{
  return this->get<bool>();
}
Parameter::operator int() const
{
  return this->get<int>();
}
Parameter::operator size_t() const
{
  return this->get<size_t>();
}
Parameter::operator unsigned int() const
{
  return this->get<size_t>();
}
Parameter::operator double() const
{
  return this->get<double>();
}
Parameter::operator std::string() const
{
  return this->get<std::string>();
}
Parameter::operator const char*() const
{
  return this->get<std::string>().c_str();
}

/* ************************************************************************** */
std::string Parameter::value_as_string() const
{
  switch(this->type)
  {
    case Parameter::types::_bool:
      return this->bool_value ? "true" : "false";
      break;
    case Parameter::types::_int:
      return convert_to_string(this->int_value);
      break;
    case Parameter::types::_size_t:
      return convert_to_string(this->unsigned_value);
      break;
    case Parameter::types::_double:
      return convert_to_string(this->double_value);
      break;
    case Parameter::types::_string:
      return this->string_value;
      break;
    default:
      ErrThrow("unknown type: ", type_as_string(this->type));
  }
}

/* ************************************************************************** */
template<> void Parameter::get_range(int& min_value, int& max_value) const
{
  int n_param = this->int_range.size();
  auto minmax = std::minmax_element(this->int_range.begin(),
                                    this->int_range.end());
  bool range_is_interval = (n_param == (*minmax.second - *minmax.first + 1));
  if(range_is_interval)
  {
    min_value = *minmax.first;
    max_value = *minmax.second;
  }
  else
  {
    // range is not an interval and this function makes no sense
    ErrThrow("integral range is not an interval");
  }
}
template<> void Parameter::get_range(size_t& min_value, size_t& max_value) const
{
  size_t n_param = this->unsigned_range.size();
  auto minmax = std::minmax_element(this->unsigned_range.begin(),
                                    this->unsigned_range.end());
  bool range_is_interval = (n_param == (*minmax.second - *minmax.first + 1));
  if(range_is_interval)
  {
    min_value = *minmax.first;
    max_value = *minmax.second;
  }
  else
  {
    // range is not an interval and this function makes no sense
    ErrThrow("integral range is not an interval");
  }
}
template<> void Parameter::get_range(double& min_value, double& max_value) const
{
  if(!check_type<double>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != double");
  min_value = this->min;
  max_value = this->max;
}

/* ************************************************************************** */
template<> void Parameter::get_range(std::set<int>& range) const
{
  range = this->int_range;
}
template<> void Parameter::get_range(std::set<size_t>& range) const
{
  range = this->unsigned_range;
}
template<> void Parameter::get_range(std::set<bool>& range) const
{
  if(this->not_bool_value_allowed)
    range = std::set<bool>( { true, false } );
  else if(this->bool_value)
    range = std::set<bool>({ true });
  else
    range = std::set<bool>({ false });
}
template<> void Parameter::get_range(std::set<double>& range) const
{
  range = std::set<double>({this->min, this->max});
}
template<> void Parameter::get_range(std::set<std::string>& range) const
{
  range = this->string_range;
}

/* ************************************************************************** */
template<> void Parameter::set(bool new_value)
{
  if(!check_type<bool>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != bool");
  // check if new_value is in the specified parameter range
  if(!not_bool_value_allowed && new_value != this->bool_value)
    ErrThrow("new parameter value out of range");
  
  this->change_count++;
  this->bool_value = new_value;
}
template<> void Parameter::set(size_t new_value)
{
  if(!check_type<size_t>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != size_t");
  // check if new_value is in the specified parameter range
  if(!check_value_in_range(new_value, this->unsigned_range))
    ErrThrow("new parameter value out of range");
  
  this->change_count++;
  this->unsigned_value = new_value;
}
template<> void Parameter::set(int new_value)
{
  if(!check_type<int>(this->type) && !check_type<size_t>(this->type))
    throw(std::runtime_error("wrong type"));
  if(check_type<size_t>(this->type))
  {
    if(new_value < 0)
      throw(std::runtime_error("wrong type"));
    // setting a positive integer to a size_t parameter
    this->set<size_t>(new_value);
    return;
  }
  // check if new_value is in the specified parameter range
  if(!check_value_in_range(new_value, this->int_range))
    ErrThrow("new parameter value out of range");
  
  this->change_count++;
  this->int_value = new_value;
}
template<> void Parameter::set(double new_value)
{
  if(!check_type<double>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != double");
  // check if new_value is in the specified parameter range
  if(new_value < this->min || new_value > this->max)
    ErrThrow("new parameter value out of range");
  
  this->change_count++;
  this->double_value = new_value;
}
template<> void Parameter::set(std::string new_value)
{
  if(!check_type<std::string>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != std::string");
  
  // check if new_value is in the specified parameter range
  if(!check_value_in_range(new_value, this->string_range))
    ErrThrow("new parameter value out of range");
  
  this->change_count++;
  this->string_value = new_value;
}
template<> void Parameter::set(const char* new_value)
{
  this->set<std::string>(std::string(new_value));
}

/* ************************************************************************** */
Parameter& Parameter::operator=(bool new_value)
{
  this->set<bool>(new_value);
  return *this;
}
Parameter& Parameter::operator=(int new_value)
{
  this->set<int>(new_value);
  return *this;
}
Parameter& Parameter::operator=(size_t new_value)
{
  this->set<size_t>(new_value);
  return *this;
}
Parameter& Parameter::operator=(double new_value)
{
  this->set<double>(new_value);
  return *this;
}
Parameter& Parameter::operator=(std::string new_value)
{
  this->set<std::string>(new_value);
  return *this;
}
Parameter& Parameter::operator=(const char* new_value)
{
  this->set<const char*>(new_value);
  return *this;
}

/* ************************************************************************** */
template<> void Parameter::set_range(int min_value, int max_value)
{
  if(!check_type<int>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != int");
  if(this->int_value < min_value || this->int_value > max_value)
    ErrThrow("current value not in specified range");
  if(max_value - min_value > 1000)
    ErrThrow("more than 1000 parameters in range");
  this->int_range.clear();
  for(int i = min_value; i <= max_value; i++)
      this->int_range.insert(i);
}
template<> void Parameter::set_range(size_t min_value, size_t max_value)
{
  if(!check_type<size_t>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != size_t");
  if(this->unsigned_value < min_value || this->unsigned_value > max_value)
    ErrThrow("current value not in specified range");
  if(max_value - min_value > 1000)
    ErrThrow("more than 1000 parameters in range");
  this->unsigned_range.clear();
  for(size_t i = min_value; i <= max_value; i++)
      this->unsigned_range.insert(i);
}
template<> void Parameter::set_range(double min_value, double max_value)
{
  if(!check_type<double>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != double");
  if(this->double_value < min_value || this->double_value > max_value)
    ErrThrow("current value not in specified range");
  this->min = min_value;
  this->max = max_value;
}
template<> void Parameter::set_range(bool min_value, bool max_value)
{ ErrThrow("Parameter::set_range with two bool arguments makes no sense"); }
template<> void Parameter::set_range(std::string min_value,
                                     std::string max_value)
{ ErrThrow("Parameter::set_range with two string arguments makes no sense"); }

/* ************************************************************************** */
template<> void Parameter::set_range(std::set<bool> new_range)
{
  if(!check_type<bool>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != bool");
  if(new_range.size() >= 2)
    this->not_bool_value_allowed = true;
  // check if (single element) range covers current value
  else if(*new_range.begin() == this->bool_value)
  {
    // only allow either true or false
    this->not_bool_value_allowed = false;
  }
  else
  {
    ErrThrow("setting a range which does not cover the current value");
  }
}
template<> void Parameter::set_range(std::set<int> new_range)
{
  if(!check_type<int>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != int");
  // check if current value is in new range
  auto it = std::find(new_range.begin(), new_range.end(), this->int_value);
  if(it == new_range.end())
    ErrThrow("setting a range which does not cover the current value");
  this->int_range = new_range;
}
template<> void Parameter::set_range(std::set<size_t> new_range)
{
  if(!check_type<size_t>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != size_ts");
  // check if current value is in new range
  auto it = std::find(new_range.begin(), new_range.end(), this->unsigned_value);
  if(it == new_range.end())
    ErrThrow("setting a range which does not cover the current value");
  this->unsigned_range = new_range;
}
template<> void Parameter::set_range(std::set<std::string> new_range)
{
  if(!check_type<std::string>(this->type))
    ErrThrow("wrong type: ", type_as_string(this->type), " != std::string");
  // check if current value is in new range
  auto it = std::find(new_range.begin(), new_range.end(), this->string_value);
  if(it == new_range.end())
    ErrThrow("setting a range which does not cover the current value");
  this->string_range = new_range;
}
template<> void Parameter::set_range(std::set<double>)
{ ErrThrow("Parameter::set_range with a set of double does not make sense"); }

/* ************************************************************************** */
void Parameter::info() const
{
  std::ostringstream os;
  os << "  Parameter: " << this->name << "\n";
  this->print_description(os, "", 60, "    description    ");
  std::string t = type_as_string(this->type);
  os << "    value(" << t << ")" << std::string(8-t.length(), ' ') 
     << this->value_as_string() << "\n";
  os << "    accessed       " << this->access_count << " times\n";
  os << "    changed        " << this->change_count << " times\n";
  os << "    range          " << this->range_as_string();
  Output::print(os.str());
}

/* ************************************************************************** */
void Parameter::print_description(std::ostream& os, std::string prepend, 
                                  size_t max_width, std::string key) const
{
  // words in description as a vector of strings
  // see http://stackoverflow.com/a/19137669
  // unfortunately this loses explicit line breaks (`\n`) within the description
  std::vector<std::string> words;
  std::istringstream iss(this->description);
  std::copy(std::istream_iterator<std::string>(iss),
            std::istream_iterator<std::string>(),
            std::back_inserter(words));
  // print the words and break the line after some fixed width
  size_t length = 0; // length of printed line (not counting indentation)
  os << prepend << key; // write the string key first
  
  size_t indent = key.length() + prepend.length();
  if(indent > max_width)
    max_width = indent; // this will print every word on a new line
  for(std::string& s : words) // loop over vector of words in description
  {
    if(length + s.length() + indent >= max_width)
    {
      // the next word has to go on a new line, add new line and indent
      os << "\n" << prepend << std::string(key.length(), ' ');
      length = 0; // reset
    }
    os << s << " ";
    length += s.length();
  }
  os << "\n";
 }