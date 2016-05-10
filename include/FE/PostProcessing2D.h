/** ***************************************************************************
 *
 * @class      PostProcessing2D
 * @brief      store given data and write output
 *             Supported formats: VTK, CASE, GNUPLOT
 *             (Note: temporary class to eventually replace Output2D)
 *
 * @author     Alfonso Caiazzo
 * @date       01.03.16
 *
 ******************************************************************************/
#ifndef __POSTPROC2D__
#define __POSTPROC2D__

#include <FEVectFunct2D.h>
#include <Vertex.h>
#include <vector>

#include <ParameterDatabase.h>

class PostProcessing2D
{
 protected:

  ///@brief output directory
  std::string testcaseDir;
  ///@brief base name of output files
  std::string testcaseName;

  ///@brief file formats option
  bool writeVTK;
  bool writeCASE;

  ///@brief array of stored FE functions
  std::vector< const TFEFunction2D* > FEFunctionArray;

  ///@brief array of stored vector-valued variables
  std::vector< const TFEVectFunct2D*> FEVectFunctArray;

  ///@brief geometry (collection of cells)
  TCollection* Coll;
  
  ///@brief start time
  double t0;

  ///@brief number of iterations between two outputs 
  int period;

  ///@brief time step
  double dt;


  ///@brief vector containing the (physical) output times
  std::vector<double> timeValues;


  ///@brief sort a list of vertices according to nodes id
  ///@todo replace the array by a vector
  void sort(TVertex **Array, int length);
  
  ///@brief get the index of a vertex in the whole list
  ///@todo this should be done using vectors
  int getIndex(TVertex **Array, int Length, TVertex *Element);
  
 public:
  ///@brief default constructor: parameter are copied from Database
  PostProcessing2D();

  /** @todo constructor that allows to specify the basename
      This is useful to create to output objects in the
      same simulation (e.g. for Stokes-Darcy problems)
   */
  //PostProcessing2D(std::string basename);
  
  ~PostProcessing2D(){};
  
  ///@brief initialize the class parameters from the Database
  void init(const ParameterDatabase& param_db); 

  ///@brief set a different testcase name (e.g. to use multiple ouput sets)
  void setTestcaseName(std::string basename);

  /// @brief add a FEFunction into this output object
  void addFEFunction(const TFEFunction2D *fefunction);

  /// @brief add a FEVectFunct into this output object
  void addFEVectFunct(const TFEVectFunct2D *fevectfunct);
  
  ///@brief write data
  void write(const char *name, int i=1, double t=0.);
  void write(std::string basename, int i=1, double t=0);
  void write(int i=1, double t=0);
  
  /** GNUPLOT */
  void writeGnuplot(const char *name);

  /** VTK (PARAVIEW) */
  void writeVtk(const char *name);


  ///@brief writes an extra vtk-file for (scalar) discontinuous functions
  void writeVtkDiscontinuous(const char *fileName, 
			     int N_LocVertices, 
			     TVertex **Vertices);

  /** @brief .case output
      Suitable for time-depending problems: The geometry can be written only once 
      (in a .geo file) while the solution (scalar and/or vector) files are written 
      at each time iteration. The .case file takes care of combining the outputs
   */
  void writeCaseFile();
  void writeCaseGeo();
  void writeCaseVars(int iter);
  

};

#endif
