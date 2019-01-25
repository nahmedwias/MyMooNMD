#ifndef INCLUDE_ASSEMBLEROUTINES_CD_LOCAL_ASSEMBLING_ROUTINES_H
#define INCLUDE_ASSEMBLEROUTINES_CD_LOCAL_ASSEMBLING_ROUTINES_H

///////////////////////////////////////////////////////////////////////////////
// standard terms (needed for Galerkin)

template <int d>
void TCDStiff(double Mult, double *coeff, double *param, double hK, double**OrigValues, 
              int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template <int d>
void TCDMass(double Mult, double *coeff, double *param, double hK, double**OrigValues, 
             int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

template <int d>
void TCDMassPOD(double Mult, double *coeff, double *param, double hK, double**OrigValues, 
             int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

template <int d>
void TCDRhs(double Mult, double *coeff, double *param, double hK, double**OrigValues, 
             int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template <int d>
void TCDStiffSUPG(double Mult, double *coeff, double *param, double hK, double**OrigValues, 
              int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template <int d>
void TCDMassSUPG(double Mult, double *coeff, double *param, double hK, double**OrigValues, 
             int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template <int d>
void TCDRhsSUPG(double Mult, double *coeff, double *param, double hK, double**OrigValues, 
             int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

#endif // INCLUDE_ASSEMBLEROUTINES_CD_LOCAL_ASSEMBLING_ROUTINES_H
