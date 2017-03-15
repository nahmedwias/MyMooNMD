/** ***************************************************************************
 *
 * @class      PostProcessing2D
 * @brief      store given data and write output
 *             Supported formats: VTK, CASE
 *             (Note: temporary class to eventually replace Output2D)
 *
 * @author     Alfonso Caiazzo
 * @date       01.03.16
 *
 ******************************************************************************/
#ifndef __POSTPROCESSING2D__
#define __POSTPROCESSING2D__

#include <vector>

#include <FEVectFunct2D.h>
#include <Vertex.h>
#include <ParameterDatabase.h>

class PostProcessing2D
{
 protected:

  // -------------------------------------------------------------------------
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

  ///@brief number of iterations between two outputs 
  size_t period;

  // -------------------------------------------------------------------------

  ///@brief array of stored FE functions
  std::vector< const TFEFunction2D* > FEFunctionArray;

  ///@brief array of stored vector-valued variables
  std::vector< const TFEVectFunct2D*> FEVectFunctArray;

  ///@brief geometry (collection of cells)
  TCollection* Coll;


  /**
     @brief vector containing the (physical) output times
     This array allows to display the physical time when looking at results
  */
  std::vector<double> timeValues;

  /**
   * This pointer points to scalar values for all cells set by
   * @a setCellValues .
   * 
   * Currently there can be just one value.
   */
  double* cellValues;
  
  /**
   * A descriptive name for @a cellValues .
   */
  std::string cellValuesName;

  ///@brief sort a list of vertices according to nodes id
  ///@todo replace the array by a vector
  void sort(TVertex **Array, int length);
  
  ///@brief get the index of a vertex in the whole list
  ///@todo this should be done using vectors
  int getIndex(TVertex **Array, int Length, TVertex *Element);
  
 public:
  ///@brief default constructor: parameter are copied from Database
  PostProcessing2D(const ParameterDatabase& param_db);
  
  /// @brief add a FEFunction into this output object
  void add_fe_function(const TFEFunction2D *fefunction);

  /// @brief add a FEVectFunct into this output object
  void add_fe_vector_function(const TFEVectFunct2D *fevectfunct);

  /**
   * @brief add values for each cell
   *
   * The input vector is assumed to be order according to the cells in the 
   * collection. The size is GetN_Cells().
   * 
   * @param in_cellValues vector containing the values for each cell
   *
   * @todo care about const correctness
   */
  void setCellValues(double *in_cellValues);
  
  void setCellValuesName(std::string &name);
  
  /**
     @brief write data to files during time dependent problems.
     You need to have added at least one fe function or fe vector function.
     @param[in] t: the physical time
  */
  void write(double t);
  /**
     @brief write data a file (mainly for stationary problems).
     You need to have added at least one fe function or fe vector function.
     If i is negative, it will not be used for the filename.
     
     If you want to write more than one case file with this function, I am not 
     sure if it works correctly.
     @param[in] i: an integer used for the naming of the file
  */
  void write(int i=-1);
  
  // -------------------------------------------------------------------------
 protected:
  /// @brief VTK output
  void writeVtk(std::string name);
  
  ///@brief writes an extra vtk-file for (scalar) discontinuous functions
  void writeVtkDiscontinuous(std::string fileName, int N_LocVertices,
                             TVertex** Vertices);

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

#endif
