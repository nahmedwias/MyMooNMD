// =======================================================================
// %W% %G%
// 
// Class:       TFEFunction3D
// Purpose:     a function from a finite element space in 3D
//
// Author:      Gunar Matthies (17.01.98)
//
// History:     start of implementation 17.01.98 (Gunar Matthies)
//
//              start of reimplementation 06.08.1998 (GM)
//
// =======================================================================

#ifndef __FEFUNCTION3D__
#define __FEFUNCTION3D__

#include <AllClasses.h>
#include <FESpace3D.h>
#include <AuxParam3D.h>
#include <Constants.h>

/** a function from a finite element space */
class TFEFunction3D
{
  protected:
    /** name of the function */
    char *Name;

    /** some more words describing the function */
    char *Description;

    /** space to which this function belongs to */
    TFESpace3D *FESpace3D;

    /** double vector according to FE isomorphism */
    double *Values;

    /** length of vector */
    int Length;

  public:
    /** constructor with vector initialization */
    TFEFunction3D(TFESpace3D *fespace3D, char *name, char *description,
                  double *values, int length);

    /** destructor */
    ~TFEFunction3D();

    /** return name */
    char *GetName()
    { return Name; }

    /** return description */
    char *GetDescription()
    { return Description; }

    /** return fe space */
    TFESpace3D *GetFESpace3D()
    { return FESpace3D; }

    /** return length */
    int GetLength()
    { return Length; }

    /** return vector of data */
    double *GetValues()
    { return Values; }

    /** calculate errors to given function */
    void GetErrors(DoubleFunct3D *Exact, int N_Derivatives,
                   MultiIndex3D *NeededDerivatives,
                   int N_Errors, ErrorMethod3D *ErrorMeth, 
                   CoeffFct3D *Coeff, TAuxParam3D *Aux,
                   int n_fespaces, TFESpace3D **fespaces,
                   double *errors);
    void GetErrorsAdapt(DoubleFunct3D *Exact, int N_Derivatives,
			MultiIndex3D *NeededDerivatives,
			int N_Errors, ErrorMethod3D *ErrorMeth, 
			CoeffFct3D *Coeff, 
			TAuxParam3D *Aux,
			int n_fespaces, TFESpace3D **fespaces,
			double *errors);
    
    /** calculate errors to given function taylored to OPTPDE */
    void GetErrorsOPTPDE(DoubleFunct3D *Exact, int N_Derivatives,
		   MultiIndex3D *NeededDerivatives,
		   int N_Errors, ErrorMethod3D *ErrorMeth, 
		   CoeffFct3D *Coeff, TAuxParam3D *Aux,
		   int n_fespaces, TFESpace3D **fespaces,
		   double radius, double upper, double lower, double *errors);
    

    /** calculate errors to given function */
    void GetMeshCellParams(DoubleFunct3D *Exact, int N_Derivatives,
                   MultiIndex3D *NeededDerivatives,
                   int N_Errors, ErrorMethod3D *ErrorMeth, 
                   CoeffFct3D *Coeff, TAuxParam3D *Aux,
                   int n_fespaces, TFESpace3D **fespaces,
                   double *errors, double *cell_parameters);

    /** determine the value of function and its first derivatives at
        the given point */
    void FindGradient(double x, double y, double z, double *values);

    /** determine the value of function and its first derivatives at
        the given point */
    void FindGradientLocal(TBaseCell *cell, int cell_no, 
                           double x, double y, double z, 
                           double *values);

    /** determine the value of function
        the given point */
    void FindValueLocal(TBaseCell *cell, int cell_no, 
                           double x, double y, double z, 
                           double *values);

    /** calculate the interpolation of an exact function */
    void Interpolate(DoubleFunct3D *Exact);

   /** calculate the super-convergence interpolation of an exact function */
    void InterpolateSuper(DoubleFunct3D *Exact);
    
    /**  interpolation of an exact function with give FeFunction values */
    void Interpolate(int N_Coord, double *Coords, int N_AuxFeFcts,  TFEFunction3D **AuxFeFcts, DoubleFunctND *Exact);
    
    
};

#endif
