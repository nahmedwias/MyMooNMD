#ifndef __PARAMETER__
#define __PARAMETER__

#include <set>
#include <string>

/** @brief store a single parameter
 * 
 * A Parameter has a value of type bool, int, size_t, double, or string which 
 * can be modified with a `set` method and returned by a `get` method. There is
 * a separate constructor for each Parameter type. Instead of the `get` method
 * casts are available as well which usually reduce the amount of code.
 * 
 * Additionally each Parameter stores a range describing a set in which valid
 * values are. It is assured that the value always is in this range. This range
 * can be either a set of discrete values, e.g. `["value1", "value2", "value3"]`
 * for a string Parameter, or an interval. While int and size_t Parameters 
 * support both kinds of ranges, double Parameters can only use intervals and
 * string and bool Parameters can only use sets of discrete values. There are 
 * `set_range` and `get_range` methods to change and inspect the range.
 * 
 * Each Parameter keeps track of the number of changes and the accesses to its
 * value, see `get_access_count` and `get_change_count`. This is mainly for 
 * debugging purposes.
 * 
 * As an identification each Parameter has a name which must not consist of 
 * spaces, dots, or colons. There is a `get_name` method.
 * 
 * In order to allow each Parameter to be documented, there is a member variable
 * `description` in the Parameter class. You can get a copy of the description
 * through the method `get_description`. 
 * 
 * To only print a summary of this Parameter use the `info` method.
 * 
 * Some remarks on the design of this class: 
 * - We did not use a template class where the template parameter is the type 
 *   of the value. The reason is that we use this Parameter class in a database
 *   which stores a set of such Parameters. However you can not create a 
 *   container of Parameters where each of the entries has a different template
 *   argument.
 * - We did not use a base class and a number of derived classes to overcome the
 *   above difficulty. This would work in principle but produced a lot of 
 *   duplicate code which we really want to avoid.
 * - We use a single class with template `get`/`set` methods to get the benefits
 *   of templates (less code) and still use many Parameters of different types
 *   in one standard container.
 * - Disadvantage: Its not clear how to store the parameter value. One idea is
 *   to store a std::string and convert it whenever accessed/changed. This is 
 *   not very nice, but allows a nice interface. A better approach could be to 
 *   use unions, however they don't support std::string easily, but only plain
 *   old data. Because the memory requirements for a single parameter is small
 *   anyway and hopefully we do not need too many of them, we decided to simply
 *   store one member variable for each possible type of parameter. This also is
 *   not nice but is easier in the implementation and maintenance. 
 */
class Parameter
{
  public:
    // the types which are supported by this class
    enum class types { _bool, _int, _size_t, _double, _string };
    

    /// @brief Create bool parameter with given name, value and description
    explicit Parameter(std::string name, bool value, std::string description);
    
    /// @brief Create int parameter with given name, value and description
    explicit Parameter(std::string name, int value, std::string description);
    
    /// @brief Create size_t parameter with given name, value and description
    explicit Parameter(std::string name, size_t value, std::string description);
    
    /// @brief Create double parameter with given name, value and description
    explicit Parameter(std::string name, double value, std::string description);
    
    /// @brief Create string parameter with given name, value and description
    explicit Parameter(std::string name, std::string value,
                       std::string description);
    

    /// @brief destructor
    ~Parameter() = default;
    
    /// @brief copy constructor
    ///
    /// This resets the access_count and the change_count to zero. These are 
    /// not tracked among copies of the same Parameter.
    ///
    /// This constructor is set as explicit to avoid accidental copies
    explicit Parameter(const Parameter&);
    
    /// @brief move constructor
    Parameter(Parameter&&) = default;
    
    /// @brief copy assignment
    Parameter& operator=(const Parameter&) = delete;
    
    /// @brief move assignment
    Parameter& operator=(Parameter&&) = default;


    /// @brief take the value and range from a given other parameter
    ///
    /// The other parameter `p` must have the same name and type. The 
    /// descriptions are concatenated (unless equal), the range is merged with 
    /// that of `p`. We think of this as if a new Parameter is defined, 
    /// therefore both change_count and access_count are reset to 0.
    void impose(const Parameter& p);
    
    /// @brief return parameter name
    std::string get_name() const;

    /// @brief return parameter description
    std::string get_description() const;

    /// @brief return access count (number of times parameter has been accessed)
    std::size_t get_access_count() const;

    /// @brief return change count (number of times parameter has been changed)
    std::size_t get_change_count() const;
    
    /// @brief range of parameter as a string
    std::string range_as_string() const;
    
    /// @brief return the value as a string
    std::string value_as_string() const;
    
    /// @brief return the type of this parameter
    Parameter::types get_type() const;

    /// @brief set range as an interval for parameter of types int, size_t, 
    /// double
    template <typename T>
    void set_range(T min_value, T max_value);

    /// @brief set range for parameter of types bool, int, size_t, std::string
    template <typename T>
    void set_range(std::set<T> range);

    /// @brief get range (interval) of parameter if T is int, size_t, or double
    ///
    /// Note that if T is int or size_t, the range of this parameter may not be
    /// be an interval, then this function throws an exception.
    template <typename T>
    void get_range(T& min_value, T& max_value) const;

    /// @brief get range (if T is floating point type, this is a lower and 
    /// upper bound)
    template <typename T>
    void get_range(std::set<T>& range) const;

    /// @brief assignment of new value
    template <typename T>
    void set(T value);
    
    /// @brief return the value of this parameter
    template <typename T>
    T get() const;
    
    /// @brief check if the value of this parameter is `value`
    template <typename T>
    bool is(T value) const;
    
    /// @name casts for easier access
    /// @brief explicit casts, throws if wrong type
    //@{
    operator bool() const;
    operator int() const;
    operator size_t() const;
    operator double() const;
    operator std::string() const;
    operator const char*() const;
    //@}
    
    /// @brief print some information on this parameter
    void info() const;
    
    /// @brief write the description of this parameter to a given stream
    ///
    /// The stream output is formatted: First the string `key` then the 
    /// description is printed where line breaks are inserted whenever the 
    /// total line width exceeds max_width. The new line is prepended by the 
    /// string `prepend` and indented by the width of the string `key`.
    ///
    /// A typical use case is `print_description(os, "# ", 60, "")` to print the
    /// description in files as a comment. See also the implementation of 
    /// Parameter::info() for another use case.
    void print_description(std::ostream& os, std::string prepend = "", 
                           size_t max_width = 80, std::string key = "") const;
  protected:
    
    /// @brief the type of this parameter (not changeable after construction)
    const types type;
    
    /// @name storing the value of this parameter
    /// @brief the value of this parameter, exactly one is actually used
    /// 
    /// Which one is used depends on Parameter::type. The others are simply 
    /// not used.
    //@{
    bool bool_value;
    int int_value;
    size_t unsigned_value;
    double double_value;
    std::string string_value;
    //@}
    
    /// @brief Access count
    mutable std::size_t access_count;

    /// @brief Change count
    std::size_t change_count;

    /// @brief Parameter key
    std::string name;

    /// @brief Parameter description
    std::string description;
    
    /// @name describing the range (valid values) of this parameter
    /// @brief Parameter range describing all valid parameter values
    ///
    /// Exactly one is actually used (int case of a Parameter of type `double`,
    /// it is `min` and `max`). The others are meaningless and are not set or 
    /// used.
    //@{
    /// @brief if false, the range is {bool_value}, otherwise it is {true,false}
    bool not_bool_value_allowed;
    /// @brief set of all possible integers for this parameter
    ///
    /// Even if the range is an interval, here all numbers are listed
    std::set<int> int_range;
    /// @brief set of all possible `size_t`s for this parameter
    ///
    /// Even if the range is an interval, here all numbers are listed
    std::set<size_t> unsigned_range;
    /// @brief the range is the interval [min, max]
    double min;
    /// @brief the range is the interval [min, max]
    double max;
    /// @brief set of all possible strings for this parameter
    std::set<std::string> string_range;
    //@}
};

#endif // __PARAMETER__
