// =======================================================================
// @(#)FEFunction2D.h        1.3 04/13/00
// 
// Class:       TFEFunction2D
// Purpose:     a function from a finite element space in 2D
//
// Author:      Gunar Matthies (17.01.98)
//
// History:     start of implementation 17.01.98 (Gunar Matthies)
//
//              start of reimplementation 06.08.1998 (GM)
//
// =======================================================================

#ifndef __FEFUNCTION2D__
#define __FEFUNCTION2D__

#include <AllClasses.h>
#include <FESpace2D.h>
#include <AuxParam2D.h>
#include <Constants.h>

/** a function from a finite element space */
class TFEFunction2D
{
  protected:
    /** name of the function */
    char *Name;

    /** some more words describing the function */
    char *Description;

    /** space to which this function belongs to */
    TFESpace2D *FESpace2D;

    /** double vector according to FE isomorphism */
    double *Values;

    /** length of vector */
    int Length;

  public:
    /** constructor with vector initialization */
    TFEFunction2D(TFESpace2D *fespace2D, char *name, char *description,
                  double *values, int length);

    /** destructor */
    ~TFEFunction2D();

    /** return name */
    char *GetName()
    { return Name; }

    /** return description */
    char *GetDescription()
    { return Description; }

    /** return fe space */
    TFESpace2D *GetFESpace2D()
    { return FESpace2D; }

    /** return length */
    int GetLength()
    { return Length; }

    /** return vector of data */
    double *GetValues()
    { return Values; }

    /** calculate errors to given function */
    void GetErrors(DoubleFunct2D *Exact, int N_Derivatives,
                   MultiIndex2D *NeededDerivatives,
                   int N_Errors, ErrorMethod2D *ErrorMeth, 
                   CoeffFct2D *Coeff, TAuxParam2D *Aux,
                   int n_fespaces, TFESpace2D **fespaces,
                   double *errors);
    
    void GetErrorsAdapt(DoubleFunct2D *Exact, int N_Derivatives,
		   MultiIndex2D *NeededDerivatives,
		   int N_Errors, ErrorMethod2D *ErrorMeth, 
		   CoeffFct2D *Coeff, TAuxParam2D *Aux,
		   int n_fespaces, TFESpace2D **fespaces,
		   double *errors);
    
    void GetErrorsOPTPDE(DoubleFunct2D *Exact, int N_Derivatives,
		   MultiIndex2D *NeededDerivatives,
		   int N_Errors, ErrorMethod2D *ErrorMeth, 
		   CoeffFct2D *Coeff, TAuxParam2D *Aux,
		   int n_fespaces, TFESpace2D **fespaces,
		   int& kink, double upper, double lower, double *errors);
    
    void GetErrorsAdaptOPTPDE(DoubleFunct2D *Exact, int N_Derivatives,
			MultiIndex2D *NeededDerivatives,
			int N_Errors, ErrorMethod2D *ErrorMeth, 
			CoeffFct2D *Coeff, TAuxParam2D *Aux,
			int n_fespaces, TFESpace2D **fespaces,
			double radius, double upper, double lower,double *errors);

    /** determine the value of function and its first derivatives at
        the given point */
    void FindGradient(double x, double y, double *values);

    /** determine the value of function and its first derivatives at
        the given point */
    void FindGradientLocal(TBaseCell *cell, int cell_no, double x, double y, double *values);

    /** determine the value of function at
        the given point */
    void FindValueLocal(TBaseCell *cell, int cell_no, double x, double y, double *values);

    /** calculate the interpolation of an exact function */
    void Interpolate(DoubleFunct2D *Exact);

    /** calculate parameters which are connected to a mesh cell */
    void GetMeshCellParams(DoubleFunct2D *Exact, int N_Derivatives,
                           MultiIndex2D *NeededDerivatives,
                           int N_Errors, ErrorMethod2D *ErrorMeth, 
                           CoeffFct2D *Coeff, 
                           TAuxParam2D *Aux,
                           int n_fespaces, TFESpace2D **fespaces,
                           double *errors, double *parameters);

    /** calculate the super-convergence interpolation of an exact function */
    void InterpolateSuper(DoubleFunct2D *Exact);
    
   /** write the solution into a data file - written by Sashi **/
   void WriteSol();

   /** Read the solution from a given data file - written by Sashi **/
   void ReadSol(char *BaseName);
   
   /** interpolate the old mesh fe function values to the new fe function **/
   void Interpolate(TFEFunction2D *OldFeFunction);

   /** sol will be correct to conserve the Old_Mass (remessing, temp, surfact, psd, etc) - added by sashi */   
   void CorrectMass(double OldMass);

   /** Retun the mass, domain volume and mean values of the function - added by sashi */
   void GetMassAndMean(double *OutVal);
};
#endif
