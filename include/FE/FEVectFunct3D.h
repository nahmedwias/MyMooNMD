// =======================================================================
// @(#)FEVectFunct3D.h        1.2 07/20/99
// 
// Class:       TFEVectFunct3D
// Purpose:     a function from a finite element space in 3D
//
// Author:      Gunar Matthies (13.07.2000)
//
// History:     start of implementation 13.07.2000 (Gunar Matthies)
//
//              WriteSol/ReadSol    13.12.10 (Sashikumaar Ganesan)
// =======================================================================

#ifndef __FEVECTFUNCT3D__
#define __FEVECTFUNCT3D__

#include <FEFunction3D.h>

/** a function from a finite element space */
class TFEVectFunct3D : public TFEFunction3D
{
  protected:
    /** number of components */
    int N_Components;

  public:

    /// Default constructor. Constructs an empty object.
    TFEVectFunct3D();

    /** constructor with vector initialization */
    TFEVectFunct3D(TFESpace3D *fespace3D, std::string name, std::string description,
                  double *values, int length, int n_components);

    /// Copy assignment operator. Shallow copy, as the
    /// FEFunction does not take any memory responsibility.
    TFEVectFunct3D& operator=( const TFEVectFunct3D & );

    /** return number of components */
    int GetN_Components()
    { return N_Components; }

    /** return i-th component as FEFunction3D */
    TFEFunction3D *GetComponent(int i) const
    {
      // the name of the component will include the index i
      std::string fct_name(Name);
      fct_name += std::to_string(i);
      return new TFEFunction3D(FESpace3D, fct_name, Description,
                               Values+i*Length, Length);
    }

    /** convert current grid to vector-values FE function */
    void GridToData();

    /** use current data for grid replacement */
    void DataToGrid();

    /** calculate errors to given vector function */
    void GetDeformationTensorErrors( 
        DoubleFunct3D *Exact, DoubleFunct3D *Exact1,
        DoubleFunct3D *Exact2,
        int N_Derivatives,
        MultiIndex3D *NeededDerivatives,
        int N_Errors, ErrorMethod3D *ErrorMeth, 
        CoeffFct3D Coeff, TAuxParam3D *Aux,
        int n_fespaces, TFESpace3D **fespaces,
        double *errors);
        
    /** write the solution into a data file **/
    void WriteSol(double t,
    		std::string directory=std::string("."),
    		std::string basename=std::string("parmoon_solution"));

    /** Read the solution from a given data file **/
    void ReadSol(std::string BaseName);
   
};

#endif
