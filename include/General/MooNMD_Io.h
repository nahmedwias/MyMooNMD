// =======================================================================
// 
// Purpose:     Collection of routine and methods for IO
//
// Author:      Gunar Matthies
//
// Date:        2004/02/19
// 
// =======================================================================

#ifndef __MOONMD_IO__
#define __MOONMD_IO__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>

using std::cout;
using std::cerr;
using std::endl;
using std::ends;
using std::setw;
using std::setprecision;


#define ErrMsg(x) { cerr <<"Error in file "<<__FILE__<<", line "<<__LINE__<<":\n"<<x<<endl; OutFile<<"Error in file "<<__FILE__<<", line "<<__LINE__<<":\n"<<x<<endl;}
#define Error(x) {ErrMsg(x)}
[[ noreturn ]] void throw_with_message(const std::string& x, std::string file,
                                       int line);
#define ErrThrow(x) {throw_with_message(std::string() + x, __FILE__, __LINE__);}

extern std::ofstream& OutFile;



/**
@brief everything related to writing to console and files
 
Especially the methods print and print_file are defined in this namespace. The
proposed usage within ParMooN is as follows:

At the beginning of any ParMooN program make a call to 

    Output::set_outfile("test_outfile.txt");
This will make sure that all calls to print will write its arguments to
`std::cout` as well as to the file 'test_outfile.txt'. Before calling 
`Output::set_outfile` the print methods only write to `std::cout`.

Use print with any number of arguments. An additional template parameter 
(`unsigned int`) is recommended to indicate at which verbosity level you 
really want to print. Higher numbers mean smaller likelihood of being 
printed. There is a threshold value stored in this namespace (Output.c++) 
which determines whether or not to print and is in the range of 1 and 5 
(maxVerbosity). For example

    print<2>("print only if verbosity level is 2 or larger")
    print<3>("print only if verbosity level is 3 or larger")
    print<4>("print only if verbosity level is 4 or larger")

Calling `Output::print` without a template argument, is equivalent to calling 
`Output::print<1>`. This means it will be printed for any verbosity level.
There is only one exception: If `Output::suppressAll()` is called, then nothing
is ever written to `std::cout` or the outfile. However in that case it is still
possible to write to a second outfile through calling redirect.


In summary the `Output::print` method writes to

- only `std::cout` if `Output::set_outfile` has not yet been called
  (discouraged)
- `std::cout` and the outfile if `Output::set_outfile` has been called
  (recommended)
- only to an extra outfile if `Output::redirect` has been called
- nowhere if `Output::suppressAll` has been called

where only in the first two cases the verbosity matters at all.

You can control the verbosity level using the methods 
`Output::increaseVerbosity()` and `Output::decreaseVerbosity()`.

At the end of the program call `Output::close_file()`.
 */
namespace Output
{
  /// @brief increase the verbosity by \p i
  void increaseVerbosity(unsigned int i = 1);
  /// @brief decrease the verbosity by \p i
  void decreaseVerbosity(unsigned int i = 1);
  /// @brief set the verbosity to a given level (must be between 1 and 5)
  void setVerbosity(unsigned int i);
  /// @brief suppress all output to std::cout and the outfile
  ///
  /// Only calls to Output::print are suppressed. Direct calls to std::cout are
  /// still possible. This lasts until Output::increaseVerbosity is called 
  /// again.
  void suppressAll();
  
  /// @brief set the outfile
  ///
  /// The print methods print to both std::cout and the outfile. This method 
  /// opens the file and it should be closed again calling Output::close_file().
  void set_outfile(std::string filename);
  
  /// @brief close files which have been previously opened
  ///
  /// Files can be opened through Output::set_outfile and Output::redirect.
  void close_file();
  
  /// \brief redirect all output to this file, nowhere else
  ///
  /// This method opens the file and all subsequent calls to print and std::cout
  /// will write into this file (no matter what the verbosity template parameter
  /// is set to). To return to the normal behavior of the print methods call
  /// Output::resetOutfile(). Typically this is used to write all output of a
  /// block of code into a separate file:
  ///     Output::redirect("file_for_output_of_f");
  ///     f(); // call some function which calls print a lot
  ///     Output::resetOutfile();
  ///
  /// \note Using `printf` will still give you output on the console
  void redirect(std::string filename);
  /// @brief close the file after a call to Output::redirect.
  ///
  /// This restores the normal behavior of the print methods.
  void resetOutfile();
  
  
  /// @brief Write something to std::cout and the outfile 
  ///
  /// This is printed only if the template parameter verbosity is not too small.
  /// That means verbosity 1 (the default) will always be printed, while with 
  /// larger numbers it is only printed if the verbosity level is high enough.
  /// Note that verbosity must be within 1 and 5 (maxVerbosity), other values
  /// (including 0) give a compile error. 
  /// 
  /// The second template parameter T can be any type which can be inserted into
  /// a std::ostream via '<<'.
  template<unsigned int verbosity = 1, typename T>
  void print(T const& t);
  
  /// @brief variadic template to enable variable number of arguments to print
  /// 
  /// This method recursively calls the itself until only two template 
  /// parameters are left, then it calls the other print method. For more 
  /// documentation, see there.
  template<unsigned int verbosity = 1, typename First, typename ... Rest>
  void print(First const& first, Rest const& ... rest);
  
  
  /// @brief print only to the outfile depending on verbosity
  ///
  /// This method behaves like the corresponding print, but does not print to
  /// std::cout, only to the outfile.
  template<unsigned int verbosity = 1, typename T>
  void printToFile(T const& t);
  /// @brief variadic template to enable variable number of arguments to print
  /// only to the outfile depending on verbosity
  ///
  /// This method behaves like the corresponding print, but does not print to
  /// std::cout, only to the outfile.
  template<unsigned int verbosity = 1, typename First, typename ... Rest>
  void printToFile(First const& first, Rest const& ... rest);
};


// this is an attempt to separate declaration and implementation. The template
// methods unfortunately can not be put in the source file, so they are 
// implemented here together with some helping methods. But in general you only
// have to worry about the above declarations.
namespace Output
{
  // maximum allowed verbosity
  constexpr unsigned int maxVerbosity = 5;
  // v should be between 1 and maxVerbosity
  bool verbose(unsigned int v);
  
  bool writeOnlyToFile();
  std::ofstream& get_outfile();
  
  // usual print to print to cout and file depending on verbosity
  template<unsigned int verbosity, typename T>
  void print(T const& t)
  {
    static_assert(verbosity > 0,
                  "It makes no sense to call print with 0 verbosity. "
                  "Use a value greater than 0");
    static_assert(verbosity <= maxVerbosity, 
                  "calling Output::print with verbosity too large");
    
    if(!writeOnlyToFile())
    {
      if(verbose(verbosity))
      {
        std::cout << t << std::endl;
        get_outfile() << t << std::endl;
      }
    }
    else
    {
      get_outfile() << t << std::endl;
    }
  }
  
  template<unsigned int verbosity, typename First, typename ... Rest>
  void print(First const& first, Rest const& ... rest)
  {
    static_assert(verbosity > 0,
                  "It makes no sense to call print with 0 verbosity"
                  "Use a value greater than 0");
    static_assert(verbosity <= maxVerbosity, 
                  "calling Output::print with verbosity too large");
    
    if(!writeOnlyToFile())
    {
      if(verbose(verbosity))
      {
        std::cout << first;
        get_outfile() << first;
        Output::print<verbosity>(rest ...);
      }
    }
    else
    {
      get_outfile() << first;
      Output::print<verbosity>(rest ...);
    }
  }
  
  
  // print only to file depending on verbosity
  template<unsigned int verbosity, typename T>
  void printToFile(T const& t)
  {
    static_assert(verbosity > 0,
                  "It makes no sense to call printToFile with 0 verbosity"
                  "Use a value greater than 0");
    static_assert(verbosity <= maxVerbosity, 
                  "calling Output::printToFile with verbosity too large");
    
    if(verbose(verbosity))
    {
      get_outfile() << t << std::endl;
    }
  }
  
  template<unsigned int verbosity, typename First, typename ... Rest>
  void printToFile(First const& first, Rest const& ... rest)
  {
    static_assert(verbosity > 0,
                  "It makes no sense to call printToFile with 0 verbosity"
                  "Use a value greater than 0");
    static_assert(verbosity <= maxVerbosity, 
                  "calling Output::printToFile with verbosity too large");
    
    if(verbose(verbosity))
    {
      get_outfile() << first;
      Output::printToFile<verbosity>(rest ...);      
    }
  }  
};

// old macro to print things into a file and t std::cout, deprecated, use 
// instead the Output::print methods below.
#define OutPut(x) {std::stringstream temporary; temporary << x; Output::print(temporary.str());}




#endif
