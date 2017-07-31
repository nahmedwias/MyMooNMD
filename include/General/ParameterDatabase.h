#ifndef __PARAMETERDATABASE__
#define __PARAMETERDATABASE__

#include <Parameter.h>
#include <list>
#include <iostream>


/** @brief store a number of Parameter objects
 * 
 * This class stores `Parameter` objects where no two Parameters have the same
 * name. One can add parameters via some `add` methods and access them via
 * `operator[]` given the name of the desired parameter as a string.
 * 
 * Furthermore this class can write all its contents to a output stream (e.g. a
 * file stream) and also read from an input stream, see methods `read` and 
 * `write` for a description of the format.
 * 
 * If one has two ParameterDatabase objects, one can merge the one into the 
 * other, see the method `merge` for details.
 * 
 * To see all stored Parameters call `ParameterDatabase::info()`.
 * 
 * Additionally to any local ParameterDatabase there is a global one which is
 * accessible from anywhere in ParMooN (if this header is included). It is 
 * called `parmoon_db`.
 */
class ParameterDatabase
{
  public:
    /// @brief construct an empty parameter database with a given name
    ParameterDatabase(std::string name);
    
    /// @brief construct a database filled with parameters of general interest
    ///
    /// These parameters include "outfile", "boundary_file", "mesh_file",
    /// "problem_type", "base_name", ...
    static ParameterDatabase parmoon_default_database();
    
    /// @brief construct a database filled with parameters for time
    /// discretization
    ///
    static ParameterDatabase default_time_database();

    /// @brief construct a database filled with parameters which holds
    /// controls and stopping criteria for a nonlinear iteration loop
    static ParameterDatabase default_nonlinit_database();

    /// @brief construct a database which holds all those parameters
    /// that are typically needed in the output methods of the system
    /// classes
    /// TODO CB The more of these default databases I add here, the more I am
    /// convinced: they are a symptom for classes which we ought to have but don't.
    static ParameterDatabase default_output_database();
    
    /// A database to control a feature that is typically used in time dependent
    /// problems. A solution can be read from a binary file and used as initial
    /// solution and the computed solution can be written to a binary file at
    /// regular intervals. TODO CB Build in the parameters needed to control
    /// the same feature in MPI.
    static ParameterDatabase default_solution_in_out_database();

    /// construct a database filled with parameters for 
    /// mesh generation using TetGen
    static ParameterDatabase default_tetgen_database();


    /// @brief delete all parameters from this database
    ~ParameterDatabase() = default;
    
    /// @brief copy constructor, this is a deep copy
    ParameterDatabase(const ParameterDatabase&);
    
    /// @brief default move constructor
    ParameterDatabase(ParameterDatabase&&) = default;
    
    /// @brief copy assignment, this is a deep copy
    ParameterDatabase& operator=(const ParameterDatabase&);
    
    /// @brief default move assignemt
    ParameterDatabase& operator=(ParameterDatabase&&) = default;
    
    
    /// @brief construct and add one parameter without given range
    ///
    /// T can be bool, int, size_t, double, or std::string.
    template <typename T>
    void add(std::string name, T value, std::string description);
    
    /// @brief construct and add one parameter with a given interval as range
    ///
    /// T can be int, size_t, or double. The range is `[min, max]`.
    template <typename T>
    void add(std::string name, T value, std::string description, T min, T max);
    
    /// @brief construct and add one parameter with a given set as range
    ///
    /// T can be bool, int, size_t, or string.
    template <typename T>
    void add(std::string name, T value, std::string description, 
             std::set<T> range);
    
    /// @brief add an already existing parameter
    void add(Parameter&& p);
    
    /// @brief return parameter with a given name
    const Parameter& operator[](std::string parameter_name) const;
    Parameter& operator[](std::string parameter_name);
    
    /// @brief return the name of this database
    std::string get_name() const;
    
    /// @brief change the name of this database
    void set_name(std::string);
    
    /// @brief get the number of parameters in this database
    size_t get_n_parameters() const;
    
    /// @brief find out if a parameter with a given name exists in this database
    /// 
    /// better name? one would write e.g.: if(db.contains("param_name"))
    bool contains(std::string name) const;
    
    /// @brief write database to a stream, which can be read again
    ///
    /// The idea is to use this function to write the database into a file which
    /// then can be used to run the program again.
    ///
    /// Set the parameter verbose to 
    /// - false, if you only want the parameters,
    /// - true, if you want the parameters, their description (documentation), 
    ///   and range.
    ///
    /// Additionally to the above the date, the database name, some hg revision 
    /// information and the host name is printed.
    void write(std::ostream& stream, bool verbose = false) const;
    
    /// @brief read parameters from a stream
    ///
    /// There are formatting restrictions. In general try to use 
    /// ParameterDatabase::write to get a conforming stream. 
    ///
    /// The name of the database is surrounded by '[' and ']' where the first
    /// character of that line must be '['. In one stream there can be multiple
    /// databases to be read. They must then be seperated by lines indicating 
    /// names, e.g.
    /// \code
    /// [name of first database]
    /// parameter_name_1: value
    /// parameter_name_2: value
    /// 
    /// [name of second database]
    /// parameter_name_3: value
    /// parameter_name_4: value
    /// \endcode
    /// The parameters for the database are all between the first occurence of 
    /// '[_some_name_]' and the second (or the end of the stream respectively).
    /// That means any parameter before the first '[_some_name_]' is ignored!
    /// If no name is found at all, then a default name is given and all lines 
    /// are considered to be read.
    ///
    /// Each parameter must be on one line. Empty lines and lines without a 
    /// colon (':') are ignored. Parameters must not have spaces in their names.
    /// The general form is one of the following
    ///     parameter_name: value
    ///     parameter_name: value [ range_min, range_max ]
    ///     parameter_name: value { range_1, range_2, range_3 }
    ///
    /// There must be a colon (':') after the name. For booleans the value can 
    /// be 'true' or 'false'.
    ///
    /// All parameters can have a range. That means the parameter will never
    /// be outside some specified values (e.g. an interval or some set of 
    /// values). Ranges for boolean parameters should be determined using braces 
    /// (e.g. '{true}' or '{ true, false }'). For parameters of type int or 
    /// size_t either braces ('{}') or brackets ('[]') are ok, one indicating a
    /// set of individual entries and the other an interval. For double 
    /// parameters only brackets (e.g. '[-2.5, 5.7]') are possible, indicating
    /// an interval. There has to be a space (' ') between the value and the 
    /// range.
    ///
    /// It is possible to read documentation for each parameter which will then
    /// become its Parameter::description. All lines directly before the line 
    /// with the parameter which start with a '#' are considered. The line 
    /// before the documentation should be either empty, some other parameter, 
    /// or consist of anything but a colon (':'). If you want to document the
    /// Parameters in your input stream, but don't want that to be read, use 
    /// two '#' instead of just one. This is how the default documentation in
    /// ParMooN can be kept inside the Parameter objects.
    ///
    /// If the input stream `is` is read until the end (this happens if there is
    /// no second name found), the `operator bool()` on `is` will return false.
    /// If it does return true, another database name has been found and a call
    /// to `std::getline(is, line)` will return a line such as 
    /// '[name_of_another_database_]'.
    void read(std::istream& is);
    
    /// @brief merge another database into this one
    ///
    /// Parameters which exist in this one get the values from the other
    /// database. All parameters in `other` which do not exist in this database,
    /// are copied into this one only if `create_new_parameters` is true.
    void merge(const ParameterDatabase &other,
               bool create_new_parameters = true);
    
    /// @brief out some information on the parameters
    ///
    /// Setting `only_names_and_values` to true prints only little information
    /// on each parameter. Setting it to false prints all details.
    void info(bool only_names_and_values = true) const;
  
  private:
    /// @brief name of this parameter database
    std::string name;
    
    // We need a container of Parameter objects with the following properties:
    // - easy insertion (easy deleting is not so much of importance)
    // - unique elements, no two Parameter objects shall have the same name
    // - full access to the Parameters (not just const)
    // A std::set<Parameter> with a custom compare function seems to be a good 
    // candidate to store all the Parameter objects. However a set will not let
    // you access the items as non-const references, as this might break the 
    // ordering. Another idea would be std::map<std::string, Parameter>, but
    // then the name of each parameter would be safed as the key and in the 
    // Parameter class, not so nice. These would have to be kept in sync. So it
    // seems std::list<Parameter> is ok as long as one makes sure no two entries
    // have the same name.
    // Note that finding a parameter in a list is slower than in a set::set or
    // std::map. However we expect this list to be rather short, so there 
    // should be no performance problem.
    
    /// @brief the set of all parameters stored in this parameter database
    std::list<Parameter> parameters;
    
};

#endif // __PARAMETERDATABASE__
