// =======================================================================
// @(#)QuadAffin.h        1.5 04/13/00
//
// Class:      TQuadAffin
//
// Purpose:    affin reference transformations for parallelogram
//
// Author:     Gunar Matthies
//
// History:    08.07.97 start implementation
// 
// =======================================================================

#ifndef __QUADAFFIN__
#define __QUADAFFIN__

#include <RefTrans2D.h>
#include <Enumerations.h>

/** reference transformations for triangle */
class TQuadAffin : public TRefTrans2D
{
  protected:
    /** x coordinate */
    double x0, x1, x2, x3;

    /** y coordinate */
    double y0, y1, y2, y3;

    /** x parameters for reference transformation */
    double xc0, xc1, xc2;

    /** y parameters for reference transformation */
    double yc0, yc1, yc2;

    /** detjk */
    double detjk;

    /** 1/detjk */
    double rec_detjk;

  public:
    /** constuctor */
    TQuadAffin();

    /** transfer form reference element to original element */
    void GetOrigFromRef(double eta, double xi, double &x, double &y) override;

    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, const double *eta, const double *xi, 
                        double *x, double *y, double *absdetjk) override;

    /** transfer form reference element to original element */
    void GetOrigFromRef(const double *ref, double *orig) override;

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double &eta, double &xi) override;

    /** transfer from original element to reference element */
    void GetRefFromOrig(const double *orig, double *ref) override;

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(BaseFunct2D BaseFunct,
                       int N_Points, const double *xi, const double *eta,
                       int N_Functs, QuadFormula2D QuadFormula);

    /** calculate functions and derivatives from reference element
        to original element, for all given elements */
    void GetOrigValues(int N_Sets, BaseFunct2D *BaseFunct,
                       int N_Points, const double *xi, const double *eta,
                       QuadFormula2D QuadFormula,
                       bool *Needs2ndDer);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(double xi, double eta, int N_BaseFunct,
                       const double *uref, const double *uxiref,
                       const double *uetaref,
                       double *uorig, double *uxorig, double *uyorig,
                       int _BaseVectDim = 1);
    
    void GetOrigValues(int joint, double zeta, int N_BaseFunct,
                       const double *uref, const double *uxiref,
                       const double *uetaref,
                       double *uorig, double *uxorig, double *uyorig,
                       int _BaseVectDim = 1);

    /** set element to cell */
    void SetCell(const TBaseCell * cell) override;

    /** return outer normal vector */
    void GetOuterNormal(int j, double zeta, double &n1, double &n2) const
      override;

    /** return tangent */
    void GetTangent(int j, double zeta, double &t1, double &t2) const override;

    /** return volume of cell according to reference transformation */
    double GetVolume() const override;
    
    
    /** transfer a  set of boundary points from reference to original element */
    void GetOrigBoundFromRef(int joint, int N_Points, double *zeta, double *X, double *Y);   
    
    /** @brief Piola transformation for vector valued basis functions */
    void PiolaMapOrigFromRef(double xi, double eta, int N_Functs,
                             const double *refD00, double *origD00) override;
    /** @brief Piola transformation for the derivatives of vector valued basis 
     *         functions */
    void PiolaMapOrigFromRef(double xi, double eta, int N_Functs,
                             const double *refD00, const double *refD10,
                             const double *refD01, double *origD10,
                             double *origD01) override;
};

#endif // __QUADAFFIN__
