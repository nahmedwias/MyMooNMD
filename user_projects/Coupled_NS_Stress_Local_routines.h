#ifndef COUPLED_NS_STRESS_LOCAL_ROUTINES_H
#define COUPLED_NS_STRESS_LOCAL_ROUTINES_H

template<int d>
void Stress_Stress(double Mult, double *coeff, double *param, double hK, double**OrigValues, 
              int *N_BaseFuncts, double ***LocMatrices, double **LocRhs, double eta);

template<int d>
void Stress_velocity(double Mult, double *coeff, double *param, double hK, double**OrigValues, 
              int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

template<int d>
void Velocity_stress(double Mult, double *coeff, double *param, double hK, double**OrigValues, 
                     int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
#endif // COUPLED_NS_STRESS_LOCAL_ROUTINES_H
