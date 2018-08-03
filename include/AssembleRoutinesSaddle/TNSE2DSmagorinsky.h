#ifndef TNSE2DSMAGORINSKY_H
#define TNSE2DSMAGORINSKY_H


// ======================================================================
// Type 1, Smagorinsky
// ======================================================================
void TimeNSType1Smagorinsky(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 2, Smagorinsky
// ======================================================================
void TimeNSType2Smagorinsky(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Smagorinsky, only nonlinear part
// Type 2, Smagorinsky, only nonlinear part
// ======================================================================
void TimeNSType1_2NLSmagorinsky(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 3, Standard Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType3Smagorinsky(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType4Smagorinsky(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 3, Standard Smagorinsky, (grad u, grad v), only nonlinear part
// Type 4, Standard Smagorinsky, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLSmagorinsky(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType3SmagorinskyDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType4SmagorinskyDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLSmagorinskyDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

#endif // TNSE2DSMAGORINSKY_H
