#ifndef __TCD2D__POD_LocalAssemble__
#define __TCD2D__POD_LocalAssemble__

#include <Enumerations.h>


/** @brief Mass matrix (p,q) for scalar problems
*/
void mat_p_q(double Mult, double *coeff, double *param, double hK, double **OrigValues, 
             int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

/** @brief Mass matrix (p,q) + SUPG term for scalar problems
*/
void mat_p_q_supg(double Mult, double *coeff, double *param, double hK, double **OrigValues, 
                  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

/** @brief Matrix (grad p, grad q) for scalar problems
*/
void mat_gradp_gradq(double Mult, double *coeff, double *param, double hK, double **OrigValues, 
                     int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

/** @brief Source term (f, q) for scalar problems
*/
void rhs_f_q(double Mult, double *coeff, double *param, double hK, double **OrigValues, 
             int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

/** @brief Source term (f, q) + SUPG terms for scalar problems
*/
void rhs_f_q_supg(double Mult, double *coeff, double *param, double hK, double **OrigValues, 
                  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

/** @brief Convection-diffusion-reaction matrix for scalar problems
*/
void mat_cdr(double Mult, double *coeff, double *param, double hK, double **OrigValues, 
             int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

/** @brief SUPG convection-diffusion-reaction matrix for scalar problems
*/
void mat_cdr_supg(double Mult, double *coeff, double *param, double hK, double **OrigValues, 
                  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

#endif
