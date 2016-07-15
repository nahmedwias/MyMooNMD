// =======================================================================
// %W% %G%
// 
// Class:       TOutput3D
// Purpose:     store given data and realize output
//
// Author:      Gunar Matthies (24.07.2000)
//
// History:     start of implementation 24.07.2000 (Gunar Matthies)
//
// =======================================================================

#ifndef __OUTPUT3D__
#define __OUTPUT3D__

#ifdef _MPI
#  include "mpi.h"
#endif
#include <string>

class TFEVectFunct3D;
class TDomain;

/** store given data and realize output */
class TOutput3D
{
  protected:
    /** collection for all spaces and functions */
    TCollection *Coll;

    /** number of stored FESpace */
    int N_FESpaces;

    /** maximal storage for FESpaces */
    int MaxN_FESpaces;

    /** array of stored FESpaces */
    const TFESpace3D **FESpaceArray;

    /** number of stored scalar variables = TFEFunction */
    int N_ScalarVar;

    /** maximal storage for FEFunction */
    int MaxN_ScalarVar;

    /** array of stored scalar variables */
    TFEFunction3D **FEFunctionArray;

    /** number of stored vector-valued variables = TFEFunction */
    int N_VectorVar;

    /** maximal storage for FEVectFunct */
    int MaxN_VectorVar;

    /** array of stored vector-valued variables */
    TFEVectFunct3D **FEVectFunctArray;

    /** number of stored paramters */
    int N_Parameters;

    /** maximal storage for parameters */
    int MaxN_Parameters;

    /** values of parameters */
    double *ParameterValues;

    /** description for parameters */
    const char **ParameterDescription;

    /** corresponding domain */
    TDomain *Domain;

    /** add a FESpace into this output object (internal use) */
    int AddFESpace(const TFESpace3D *fespace);

     /** internal data storage **/
    struct TOutputData
    {
      int N_Nodes;
      int N_Data;
      TVertex **Nodes;
      int *ConList;

      double **FEFuncValues;

      enum CELLTYPE {TETRAHEDRON=4, BRICK} Type;

      TOutputData() : N_Nodes(0), N_Data(0), Nodes(0), ConList(0),
		       FEFuncValues(0) {};
      ~TOutputData();
    };

    TOutputData *Data;

    char *Name;
  
  protected:
    void ComputeOutputData();
    void ComputeFEValues();

  public:
    /** constructor: maximum number of these things */
    TOutput3D(int maxn_fespaces, int maxn_scalar, int maxn_vect, 
              int maxn_parameters, TDomain *domain, TCollection *coll=NULL,
	      const char *name=NULL);

    /** destructor: freeing all allocated space */
    ~TOutput3D();

    /** add a FEFunction into this output object */
    int AddFEFunction(TFEFunction3D *fefunction);

    /** add a FEVectFunct into this output object */
    int AddFEVectFunct(TFEVectFunct3D *fevectfunct);

    /** add parameter into this output object */
    int AddParameter(double value, const char *descr);
    
    /** write stored data into an vtk file */
    int WriteVtk(const char *name);
    
     /** write a discontinuous function into a VTK file */
    void WriteVtkDiscontinuous(const char *fileName, 
                               int N_LocVertices, TVertex **Vertices);

    /** write stored PARALLEL data into a pvtu and vtu files (XML files for paraview) */
    int Write_ParVTK(
#ifdef _MPI
                     MPI_Comm comm,
#endif
                     int img, char *subID,
                     std::string directory = std::string("."),
                     std::string basename = std::string("parmoon_solution"));
};

#endif
