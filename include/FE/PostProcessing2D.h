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
  bool listCreated;

  /// @brief number of iterations between two outputs 
  size_t n_steps_per_output;
  /// @brief counting the steps until next output is written
  mutable int n_steps_until_next_output;

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

  /** @brief writes the coordinates of the vertices (used in Case and VTK)
   */
  void writeCoord(std::ofstream & f, int dimension);
  
  /** @brief method computes the value in of a given function in each vertex
 * Input FEFunction2D or FEVectFunction2D function, 
 *empty double vector, just declared int dimension 
 * Returns vector of size (dimension)*(#nodes), where dimension is the number of
 * components of the values of the function, i.e. GetBaseVectDim() or GetN_Components()
 */
  template <class T>
  void computeNodeValues(const T* function,
			 std::vector<double>& solutionAtNode, unsigned int &dimension);
  
  /** @brief auxiliary function for writeVectCase */
  void printEntry(std::ofstream & dat, double value, int counter);
   
  /** @brief writes function values in an suitable output for case */
  void writeVectCase(std::ofstream & dat,
		       unsigned int N_Vertices, unsigned int N_Comp,
		       std::vector<double> solutionAtNode);
   
  /** @brief writes x-component for all vertices and then y- component for all vertices*/
  void printVectCompwise(std::ofstream & dat, std::string name,
		        unsigned int N_Vertices, unsigned int N_Comp,
			 std::vector<double> solutionAtNode);
  
  /** @brief writes the vector (x- and y-component) for all vertices*/
  void printVectPointwise(std::ofstream & dat, std::string name,
		        unsigned int N_Vertices, unsigned int N_Comp,
			  std::vector<double> solutionAtNode);
    
  /** @brief writes the absolute value for all vertices*/
  void printVectAbsValue(std::ofstream & dat, std::string name,
		        unsigned int N_Vertices, unsigned int N_Comp,
			 std::vector<double> solutionAtNode);
  

 public:
  ///@brief default constructor: parameter are copied from Database
  PostProcessing2D(const ParameterDatabase& param_db);
  
  /// @brief add a FEFunction into this output object
  void add_fe_function(const TFEFunction2D *fefunction);

  /// @brief add a FEVectFunct into this output object
  void add_fe_vector_function(const TFEVectFunct2D *fevectfunct);
  
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
  
  // -------------------------------------------------------------------------
 protected:
  /// @brief VTK output
  void writeVtk(std::string name);
  
  ///@brief writes an extra vtk-file for (scalar) discontinuous functions
  void writeVtkDiscontinuous(std::string fileName, int N_LocVertices,
                             std::vector<TVertex *> Vertices);

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
