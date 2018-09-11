#include <Database.h>
#include <FEDatabase2D.h>
#include <FEFunction2D.h>
#include <LinAlg.h>
#include <MainUtilities.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <MooNMD_Io.h>
#include <ConvDiff.h>

/** @brief the local assembling routines */

/**************************************************************************** */
void mat_p_q(double Mult, double *coeff,
              double *param, double hK,
              double **OrigValues, int *N_BaseFuncts,
              double ***LocMatrices, double **LocRhs)
{
  double ansatz00, test00, val;
  double **MatrixC = LocMatrices[0];
  int N_P = N_BaseFuncts[0];
  double *Orig  = OrigValues[0];          // p

  for(int i=0;i<N_P;i++)
  {
    test00 = Orig[i];
    for(int j=0;j<N_P;j++)
    {
      ansatz00 = Orig[j];
      val = Mult * test00*ansatz00;
      MatrixC[i][j] += val;
    }
  }
}

void mat_p_q_supg(double Mult, double *coeff,
        double *param, double hK,
        double **OrigValues, int *N_BaseFuncts,
        double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c5;
  double delta, bgradv;

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];            // u
  Orig1 = OrigValues[1];            // u_x
  Orig2 = OrigValues[2];            // u_y

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c5 = coeff[5];                                  // \|b\|_infty

  //NOTE definition of SUPG parameter as for a stationary problem
  delta = 0;//Compute_SDFEM_delta<2>(hK, c0, c1, c2, c3, c5);

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test00 = Orig0[i];
    test10 = Orig1[i];
    test01 = Orig2[i];
    bgradv = (c1*test10+c2*test01) * delta;
    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig0[j];
      // (p,q)
      val  = ansatz00 * test00;
      // supg part
      val += ansatz00 * bgradv;
      val *= Mult;
      MatrixRow[j] += val;
    }                                             // endfor j
  }
}

/**************************************************************************** */
void mat_gradp_gradq(double Mult, double *coeff,
              double *param, double hK,
              double **OrigValues, int *N_BaseFuncts,
              double ***LocMatrices, double **LocRhs)
{
  double **MatrixC;
  double ansatz10, ansatz01;
  double test10, test01;
  double *Orig1, *Orig2;
  int N_P;

  // matrix for pressure/pressure term
  MatrixC = LocMatrices[0];
  N_P = N_BaseFuncts[0];


  Orig1 = OrigValues[1];         // p_x
  Orig2 = OrigValues[2];         // p_y


  for(int i=0;i<N_P;i++)
  {
    test10 = Orig1[i];
    test01 = Orig2[i];

    for(int j=0;j<N_P;j++)
    {
      ansatz10 = Orig1[j];
      ansatz01 = Orig2[j];

      MatrixC[i][j] += Mult * (test10*ansatz10 + test01*ansatz01);
    }
  }
}

/**************************************************************************** */
void rhs_f_q(double Mult, double *coeff,
             double *param, double hK,
             double **OrigValues, int *N_BaseFuncts,
             double ***LocMatrices, double **LocRhs)
{
  double *Orig0;
  double test00;
  double c1;
  int N_P;

  N_P = N_BaseFuncts[0];
  Orig0 = OrigValues[0]; // p

  c1 = coeff[4];         // f

  for(int i=0;i<N_P;i++)
  {
    test00 = Orig0[i];
    LocRhs[0][i] += Mult*test00*c1;
  }
}

/**************************************************************************** */
void rhs_f_q_supg(double Mult, double *coeff,
                  double *param, double hK,
                  double **OrigValues, int *N_BaseFuncts,
                  double ***LocMatrices, double **LocRhs)
{
  double *Orig0, *Orig1, *Orig2;
  double test00, test10, test01;
  double c0, c1, c2, c3, c4, c5;
  double bgradv, delta;
  int N_P;

  N_P = N_BaseFuncts[0];
  Orig0 = OrigValues[0];           // p
  Orig1 = OrigValues[1];           // p_x
  Orig2 = OrigValues[2];           // p_y

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  c5 = coeff[5];                                  // \|b\|_infty

  //NOTE definition of SUPG parameter as for a stationary problem
  delta = 0;//Compute_SDFEM_delta<2>(hK, c0, c1, c2, c3, c5);

  for(int i=0;i<N_P;i++)
  {
    test00 = Orig0[i];
    test10 = Orig1[i];
    test01 = Orig2[i];
    bgradv = (c1*test10+c2*test01) * delta;
    LocRhs[0][i] += Mult*c4*bgradv;
    LocRhs[0][i] += Mult*test00*c4;
  }
}

/**************************************************************************** */
void mat_cdr(double Mult, double *coeff, double *param, double hK,
             double **OrigValues, int *N_BaseFuncts,
             double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3;

  Matrix = LocMatrices[0];
  N_ = N_BaseFuncts[0];
  Orig0 = OrigValues[0];                         // p
  Orig1 = OrigValues[1];                         // p_x
  Orig2 = OrigValues[2];                         // p_y

  // coefficients
  c0 = coeff[0];                                  // eps
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c


  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    // test function
    test10 = Orig1[i];                            // xi derivative
    test01 = Orig2[i];                            // eta derivative
    test00 = Orig0[i];                            // function

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig1[j];                        // xi derivative
      ansatz01 = Orig2[j];                        // eta derivative
      ansatz00 = Orig0[j];                        // function

      // assemble viscous term
      // eps (test_x ansatz_x + test_y ansatz_y)
      val = c0*(test10*ansatz10+test01*ansatz01);
      // assemble convective term
      // (b_1 ansatz_x + b_2 ansatz_y) test
      val += (c1*ansatz10+c2*ansatz01)*test00;
      // assemble reactive term
      // c  ansatz test
      val += c3*ansatz00*test00;

      // quad weigth
      val *= Mult;

      // update matrix entry
      MatrixRow[j] += val;
    }                                             // endfor j
  }                                               // endfor i
}

/**************************************************************************** */
void mat_cdr_supg(double Mult, double *coeff, double *param,
                          double hK,
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c5;
  double delta, bgradv;
  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  //c4 = coeff[4];                                  // f
  c5 = coeff[5];                                  // \|b\|_infty

  //NOTE definition of SUPG parameter as for a stationary problem
  delta = 0;//Compute_SDFEM_delta<2>(hK, c0, c1, c2, c3, c5);

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test00 = Orig0[i];
    test10 = Orig1[i];
    test01 = Orig2[i];

    bgradv = (c1*test10+c2*test01) * delta;

    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig0[j];
      ansatz10 = Orig1[j];
      ansatz01 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      val = c0*(test10*ansatz10+test01*ansatz01); //diffusion
      val += (c1*ansatz10+c2*ansatz01)*test00; //convection
      val += c3*ansatz00*test00; //reaction

      val += (-c0*(ansatz20+ansatz02) //supg
          +c1*ansatz10+c2*ansatz01
          +c3*ansatz00) * bgradv;

      val *=Mult;

      MatrixRow[j] += val;

    }                                             // endfor j
  }                                               // endfor i
}
