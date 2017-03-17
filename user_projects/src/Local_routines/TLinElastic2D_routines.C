/*
 * TLinElastic2D_routines.C
 *
 *  Created on: Mar 17, 2017
 *      Author: alia
 */


#include <Database.h>
#include <Convolution.h>
#include <MooNMD_Io.h>

#include <stdlib.h>


// ======================================================================
// Stiffness, mass matrix and rhs
// ======================================================================
void TimeLinearElasticityWholeSystem(double Mult, double *coeff,
                                     double *param, double hK,
                                     double **OrigValues, int *N_BaseFuncts,
                                     double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22, **MatrixM11, **MatrixM22;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow11, *MatrixRow12, *MatrixRow21, *MatrixRow22,
  *MatrixMRow1, *MatrixMRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2; // corresponds to number of terms
  int i,j,N_U;
  double c2, c3;         // c0, c1;      // coeff from get_coeff
  double lambda, mu, rho;//, u1, u2;     // parameters

  MatrixA11  = LocMatrices[0];
  MatrixA12  = LocMatrices[1];
  MatrixA21  = LocMatrices[2];
  MatrixA22  = LocMatrices[3];
  MatrixM11  = LocMatrices[4];
  MatrixM22  = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];   // u,x  derivative
  Orig1 = OrigValues[1];   // u,y
  Orig2 = OrigValues[2];   // u

  /* note: the lambda and mu, here, are not taken from the example with
   * get_coeff, they are rather used with param.
   */
//  c0 = coeff[0];           // lambda
//  c1 = coeff[1];           // mu
  c2 = coeff[2];           // f1
  c3 = coeff[3];           // f2

//  u1     = param[0];       // u1old not needed in most cases
//  u2     = param[1];       // u2old not needed in most cases
  lambda = param[2];       // lame coefficient
  mu     = param[3];       // lame coefficient
  rho    = param[4];       // material density

  for(i=0;i<N_U;i++)
  {
    MatrixRow11 = MatrixA11[i];
    MatrixRow12 = MatrixA12[i];
    MatrixRow21 = MatrixA21[i];
    MatrixRow22 = MatrixA22[i];
    MatrixMRow1 = MatrixM11[i];
    MatrixMRow2 = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c2;
    Rhs2[i] += Mult*test00*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = mu*(2*test10*ansatz10+test01*ansatz01);
      val += lambda*(test10+test01)*(ansatz10+ansatz01);
      MatrixRow11[j] += Mult * val;

      val  = mu*(test01*ansatz10);
      MatrixRow12[j] += Mult * val;

      val  = mu*(test10*ansatz01);
      MatrixRow21[j] += Mult * val;

      val  = mu*(test10*ansatz10+2*test01*ansatz01);
      val += lambda*(test10+test01)*(ansatz10+ansatz01);
      MatrixRow22[j] += Mult * val;

      val = rho*Mult*(ansatz00*test00);
      MatrixMRow1[j] += val;
      MatrixMRow2[j] += val;
    } // endfor j
//    cout << "lambda = " << lambda << " mu =  " << mu << " rho = " << rho << endl;
  } // endfor i
}



// ======================================================================
// Stiffness matrix
// ======================================================================
void TimeLinearElasticityStiffness(double Mult, double *coeff,
                                   double *param, double hK,
                                   double **OrigValues, int *N_BaseFuncts,
                                   double ***LocMatrices, double **LocRhs)
{

}

// ======================================================================
// Mass matrix
// ======================================================================
void TimeLinearElasticityMass(double Mult, double *coeff,
                                   double *param, double hK,
                                   double **OrigValues, int *N_BaseFuncts,
                                   double ***LocMatrices, double **LocRhs)
{

}

// ======================================================================
// Rhs
// ======================================================================
void TimeLinearElasticityRhs(double Mult, double *coeff,
                                   double *param, double hK,
                                   double **OrigValues, int *N_BaseFuncts,
                                   double ***LocMatrices, double **LocRhs)
{

}


// ======================================================================
// Parameter Fct
// ======================================================================
void ParameterFunction(double *in, double *out)
{
  // in[0] and in[1] are x and y coordinates
  out[0] = in[2];   //  u1
  out[1] = in[3];   //  u2
  out[2] = in[4];   //  lambda
  out[3] = in[5];   //  rho
  out[4] = in[6];   //  mu
}


