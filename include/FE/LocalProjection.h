// =======================================================================
// LocalProjection.h
//
// Purpose:   routines for local projection stabilization
//
// Author:    Gunar Matthies  2007/03/06
//
// =======================================================================

#ifndef __LOCAL_PROJECTION__
#define __LOCAL_PROJECTION__

#include <LinAlg.h>

#ifdef __2D__
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *C,
        double *x, double *b, double *r);

void Defect_NSE2C(TSquareMatrix **A, TMatrix **B, double *x,
                  double *b, double *r);

void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *C,
        double *x, double *y);

void MatVect_NSE2C(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void UltraLocalProjection(void* A, boolean ForPressure, CoeffFct2D *Coeff);
void UltraLocalProjection(void* A, boolean ForPressure);
void UltraLocalProjectionSD(void* A, boolean ForPressure);

double UltraLocalError(TFEFunction2D *uh, DoubleFunct2D *ExactU,
        double lpcoeff, double lpexponent, int OrderDiff);

double UltraLocalErrorDivergence(TFEFunction2D *uh1, TFEFunction2D *uh2,
                       DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2,
                       double lpcoeff, double lpexponent, int OrderDiff);

double UltraLocalErrorStreamline(TFEFunction2D *uh, DoubleFunct2D *ExactU,
                       TFEFunction2D *b1, TFEFunction2D *b2,
                       double lpcoeff, double lpexponent, int OrderDiff);

double UltraLocalErrorStreamlinePWConst(TFEFunction2D *uh, DoubleFunct2D *ExactU,
                       TFEFunction2D *b1, TFEFunction2D *b2,
                       double lpcoeff, double lpexponent, int OrderDiff);

void AddStreamlineTerm(TSquareMatrix2D* A, TFEFunction2D *uh1,
                       TFEFunction2D *uh2,
                       double lpcoeff, double lpexponent, int OrderDiff);

void AddStreamlineTermPWConst(TSquareMatrix2D* A, TFEFunction2D *uh1,
                              TFEFunction2D *uh2,
                              double lpcoeff, double lpexponent, int OrderDiff);

void AddDivergenceTerm(TSquareMatrix2D *A11,TSquareMatrix2D *A12,
                       TSquareMatrix2D *A21,TSquareMatrix2D *A22,
                       double lpcoeff, double lpexponent, int OrderDiff);

void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                    TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                    TMatrix *B1T, TMatrix *B2T,
                    TMatrix *C,
                    double *x, double *y);

void MatVect_NSE4C(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                   TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                   TMatrix *B1T, TMatrix *B2T,
                   TMatrix *C,
                   double *x, double *b, double *r);

void Defect_NSE4C(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

void AdaptivePostProcess(TFEFunction2D *FeFunction, double *PostSol, boolean DirichletBC);


void AddALEStreamlineLPS(TSquareMatrix2D* A, int N_FeFunct, TFEFunction2D **FeFunct,
                         double lpcoeff, double lpexponent, int OrderDiff);
#else 

void AddStreamlineTerm(TSquareMatrix3D* A, TFEFunction3D *uh1,
                       TFEFunction3D *uh2, TFEFunction3D *uh3,
                       double lpcoeff, double lpexponent, int OrderDiff); 

void UltraLocalProjection(TSquareMatrix3D* A, 
                          double lpcoeff, double lpexponent, int OrderDiff);

#endif                           
#endif
