// =======================================================================
// @(#)QuadBilinear.h        1.5 04/13/00
//
// Class:      TQuadBilinear
//
// Purpose:    Bilinear reference transformations for parallelogram
//
// Author:     Gunar Matthies
//
// History:    08.07.97 start implementation
// 
// =======================================================================

#ifndef __QUADBilinear__
#define __QUADBilinear__

#include <RefTrans2D.h>

/** reference transformations for triangle */
class TQuadBilinear : public TRefTrans2D
{
  protected:
    /** x coordinate */
    double x0, x1, x2, x3;

    /** y coordinate */
    double y0, y1, y2, y3;

    /** x parameters for reference transformation */
    double xc0, xc1, xc2, xc3;

    /** y parameters for reference transformation */
    double yc0, yc1, yc2, yc3;

  public:
    /** constuctor */
    TQuadBilinear();

    /** transfer form reference element to original element */
    void GetOrigFromRef(double eta, double xi, double &x, double &y);

    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, double *eta, double *xi, 
                        double *x, double *y, double *absdetjk);

    /** transfer form reference element to original element */
    void GetOrigFromRef(double *ref, double *orig);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double &eta, double &xi);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double *orig, double *ref);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(BaseFunct2D BaseFunct,
                       int N_Points, double *xi, double *eta, 
                       int N_Functs, QuadFormula2D QuadFormula);

    /** calculate functions and derivatives from reference element
        to original element, for all given elements */
    void GetOrigValues(int N_Sets, BaseFunct2D *BaseFunct,
                       int N_Points, double *xi, double *eta,
                       QuadFormula2D QuadFormula,
                       bool *Needs2ndDer);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(double xi, double eta, int N_BaseFunct,
                       double *uref, double *uxiref, double *uetaref,
                       double *uorig, double *uxorig, double *uyorig,
                       int _BaseVectDim = 1);
    
    void GetOrigValues(int joint, double zeta, int N_BaseFunct,
                       double *uref, double *uxiref, double *uetaref,
                       double *uorig, double *uxorig, double *uyorig,
                       int _BaseVectDim = 1);

    /** set element to cell */
    void SetCell(const TBaseCell * cell);

    /** return outer normal vector */
    void GetOuterNormal(int j, double zeta,
                                double &n1, double &n2);

    /** return tangent normal vector */
    void GetTangent(int j, double zeta,
                                double &t1, double &t2);
    /** return volume of cell according to reference transformation */
    double GetVolume();
    
    /** @brief Piola transformation for vector valued basis functions */
    void PiolaMapOrigFromRefNotAffine(int N_Functs, double *refD00, 
                                      double *origD00, double xi, double eta);
    /** @brief Piola transformation for the derivatives of vector valued basis 
     *         functions */
    void PiolaMapOrigFromRefNotAffine(int N_Functs, double *refD00, 
                                      double *refD10, double *refD01, 
                                      double *origD10, double *origD01,
                                      double xi, double eta);
};

#endif
