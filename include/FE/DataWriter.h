/** ***************************************************************************
 *
 * @class      DataWriter
 * @brief      store given data and write output
 *             Supported formats: VTK, CASE
 *
 * @author     Claudia Wohlgemuth, Ulrich Wilbrandt
 * @date       18.05.18
 *
 ******************************************************************************/
#ifndef __DATAWRITER__
#define __DATAWRITER__

#include <ParameterDatabase.h>
#include <vector>
#ifdef _MPI
#include "mpi.h"
#endif

class TVertex;
class TCollection;
class TFEFunction2D;
class TFEFunction3D;
class TFEVectFunct2D;
class TFEVectFunct3D;

template <int d>
struct Function_names
{
};
template <>
struct Function_names<2>
{
  typedef TFEFunction2D FEFunction;
  typedef TFEVectFunct2D FEVectFunct;
};
template <>
struct Function_names<3>
{
  typedef TFEFunction3D FEFunction;
  typedef TFEVectFunct3D FEVectFunct;
};

template <int d>
class DataWriter
{
  public:
  ///@brief default constructor: parameter are copied from Database
  DataWriter(const ParameterDatabase& param_db);

  using FEFunction = typename Function_names<d>::FEFunction;
  using FEVectFunct = typename Function_names<d>::FEVectFunct;

  /// @brief add a FEFunction into this output object
  void add_fe_function(const FEFunction* fefunction);

  /// @brief add a FEVectFunct into this output object
  void add_fe_vector_function(const FEVectFunct* fevectfunct);

  /**
     @brief write data to files during time dependent problems.
     You need to have added at least one fe function or fe vector function.
     @param[in] t the physical time
  */
  void write(double t);
  /**
     @brief write data a file (mainly for stationary problems).
     You need to have added at least one fe function or fe vector function.
  */
  void write();

  /**
   * @brief Write mesh in VTK format
   */
  void writeMesh(std::string name);

  // -------------------------------------------------------------------------
  protected:
  ///@brief output directory
  std::string testcaseDir;
  ///@brief base name of output files
  std::string testcaseName;

  /**
     @brief file formats option
     ParMooN currently supports output in:
     -VTK (one file per time step, containing the geometry and all the variables
     -CASE (one file per variable/time step, separate geometry file)
     CASE output is suitable for time dependent (esp. large) simulations, as the
     geometry does not need to be rewritten.
     Both formats are supported by Paraview
  */
  bool writeVTK;
  bool writeCASE;

  /// @brief number of iterations between two outputs
  size_t n_steps_per_output;
  /// @brief counting the steps until next output is written
  mutable int n_steps_until_next_output;

  ///@brief array of stored FE functions
  std::vector<const FEFunction*> FEFunctionArray;
  ///@brief array of stored vector-valued variables
  std::vector<const FEVectFunct*> FEVectFunctArray;

  ///@brief geometry (collection of cells)
  TCollection* Coll;

  /**
     @brief vector containing the (physical) output times
     This array allows to display the physical time when looking at results
  */
  std::vector<double> timeValues;

  // -------------------------------------------------------------------------

  /** @brief writes the coordinates of the vertices (used in Case and VTK)
   */
  void writeCoord(std::ofstream& f);

  /** @brief method computes the value in of a given function in each vertex
 * Input FEFunction2D or FEVectFunction2D function,
 *empty double vector, just declared int dimension
 * Returns vector of size (dimension)*(#nodes), where dimension is the number of
 * components of the values of the function, i.e. GetBaseVectDim() or
 *GetN_Components()
 */
  template <class T>
  void computeNodeValues(const T* function, std::vector<double>& solutionAtNode,
                         unsigned int& dimension);

  /** @brief writes function values in an suitable output for case */
  void writeVectCase(std::ofstream& dat, unsigned int N_Vertices,
                     unsigned int N_Comp, std::vector<double> solutionAtNode);

  /** @brief writes x-component for all vertices and then y- component for all
   * vertices*/
  void printVectCompwise(std::ofstream& dat, std::string name,
                         unsigned int N_Vertices, unsigned int N_Comp,
                         std::vector<double> solutionAtNode);

  /** @brief writes the vector (x- and y-component) for all vertices*/
  void printVectPointwise(std::ofstream& dat, std::string name,
                          unsigned int N_Vertices, unsigned int N_Comp,
                          std::vector<double> solutionAtNode);

  /** @brief writes the absolute value for all vertices*/
  void printVectAbsValue(std::ofstream& dat, std::string name,
                         unsigned int N_Vertices, unsigned int N_Comp,
                         std::vector<double> solutionAtNode);


  /// @brief VTK output
  void writeVtk(std::string name);

  ///@brief writes an extra vtk-file for (scalar) discontinuous functions
  void writeVtkDiscontinuous(std::string fileName, int N_LocVertices,
                             std::vector<const TVertex*> Vertices);


/** write stored PARALLEL data into a pvtu and vtu files (XML files for
 * paraview) */
#ifdef _MPI
  void Write_ParVTK(MPI_Comm comm, int img, char* subID,
                    std::string directory = std::string("."),
                    std::string basename = std::string("parmoon_solution"));
#endif


  // -------------------------------------------------------------------------
  /** @brief: .case output, suitable for time-depending problems.
      This format write the output onto multiple files:
      - The geometry is stored on a basename.[timestep].geo file
      - The results are stored on basename_fct.timestep.{scl,vct} files,
      where fct identify the written fied (e.g. temp, conc, pres, vel), scl
      stays for 'scalar', vct for 'vector'
      - The 'controller' file is a .case (text) file, which is used by Paraview
      to display the results.

      The case file will be rewritten every time this function is called. That
      enables us to view this file during a time consuming simulation.
   */
  void writeCaseFile();
  void writeCaseGeo();
  void writeCaseVars(int iter);
};

typedef DataWriter<3> DataWriter3D;
typedef DataWriter<2> DataWriter2D;

#endif // __DATAWRITER__
