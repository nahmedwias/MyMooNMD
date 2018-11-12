// =======================================================================
// @(#)QuadFormula.h        1.3 05/04/99
//
// Class:    TQuadFormula
//
// Purpose:  base class for quadrature formulas
// Author:   Gunar Matthies
//
// History:  29.08.1997 start implementation
//           07.04.1999 add Accuracy
// 
// =======================================================================

#ifndef __QUAD_FORMULA__
#define __QUAD_FORMULA__

#include <Constants.h>
#include <MooNMD_Io.h>

/** base class for quadrature formulas */
class TQuadFormula
{
  protected:
    /** number of quadrature points */
    int N_QuadPoints;
    /** weights for the formula */
    double *Weights;

    /** accuracy of this formula */
    int Accuracy;

  protected:
    /** constructor */
    TQuadFormula();

  public:
    /** return number of quadrature points */
    int GetN_QuadPoints() const
    { return N_QuadPoints; }

    /** return weights of the formula */
    const double *GetWeights() const
    { return Weights; }
    /** return coordinates of the formula */
    virtual const double *GetCoords(int i);

    /** print information on this formula */
    friend std::ostream & operator << (std::ostream &s, const TQuadFormula *qf);
};

#endif
