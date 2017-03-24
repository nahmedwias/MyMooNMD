#include <assemble_routine_tnse2D_PBVMS.h>
#include <Database.h>
#include <CommonRoutineTNSE2D.h>

//================================================================================
// Local assemble routines for the Projection-Based VMS method:: It's implemented
// in the following way, the test nonlinearity is dealt using extrapolation
// of velocity and pressure, whereas the ansatz are deal with the fix-point
// iterations.
//================================================================================

//================================================================================
void TimeNSParamsVelo_GradVelo_LargeScale2D(double *in, double *out)
{
  out[0] = in[2]; // u1old
  out[1] = in[3]; // u2old

  out[2] = in[4]; // D1u1
  out[3] = in[5]; // D1u2
  out[4] = in[6]; // D2u1
  out[5] = in[7]; // D2u2

  out[6] = in[0]; // x - coordinate for van Driest damping
  out[7] = in[1]; // y - coordinate for van Driest damping

  out[8] = in[8]; // projection space label
}

//================================================================================
/* ======================================================================
 * Type 4, VMS_Projection, D(u):D(v)
 * ======================================================================
 */
void TimeNSType4VMS_ProjectionDD2D(double Mult, double *coeff,
                double *param, double hK,
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21;//**MatrixA13,
  double **MatrixA22;//, **MatrixA23, **MatrixA31, **MatrixA32;
//  double **MatrixA33;
  double **MatrixM11;
  double **MatrixB1, **MatrixB2;//,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T;//,  **MatrixB3T;
  double **MatrixL, **Matrix_tilde_G11, **Matrix_tilde_G22;//, **Matrix_tilde_G33;
  double **Matrix_G11, **Matrix_G22;//,   **Matrix_G33;
  double *Rhs1, *Rhs2, val, val1;//*Rhs3,
  double *Matrix11Row, *Matrix12Row, *Matrix21Row;//*Matrix13Row,
  double *Matrix22Row;//, *Matrix23Row, *Matrix31Row, *Matrix32Row;
//  double *Matrix33Row;
  double *MatrixM11Row;
  double *MatrixRow1, *MatrixRow2;//, *MatrixRow3;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;//, *Orig5;
  int i,j,N_U, N_P, N_L;
  double c0, c1, c2;//, c3;
  double u1, u2, mu, viscosity;//u3

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  //MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  //MatrixA23 = LocMatrices[5];
  //MatrixA31 = LocMatrices[6];
  //MatrixA32 = LocMatrices[7];
  //MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[4];
  MatrixL   = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];
  //MatrixB3  = LocMatrices[13];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];
  //MatrixB3T = LocMatrices[16];
  Matrix_tilde_G11  = LocMatrices[10];
  Matrix_tilde_G22  = LocMatrices[11];
  //Matrix_tilde_G33  = LocMatrices[19];
  Matrix_G11  = LocMatrices[12];
  Matrix_G22  = LocMatrices[13];
  //Matrix_G33  = LocMatrices[22];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  //Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  //Orig2 = OrigValues[2]; // u_z
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p
  Orig4 = OrigValues[4]; // l

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  //c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  //u3 = param[2]; // u3old


//  double *x = &param[12];
//  double *y = &param[13];
//  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[2];
  double *uConv = &param[0];
//  double *projection_space_label = &param[15];
  mu = turbulentviscosity(hK, gradu, u, uConv);

  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
//    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
//    Matrix23Row = MatrixA23[i];
//    Matrix31Row = MatrixA31[i];
//    Matrix32Row = MatrixA32[i];
//    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    val1 = Mult*test00;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
//    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;

      val  = viscosity*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

//      val  = viscosity*(test001*ansatz100);
//      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

//      val  = viscosity*(test001*ansatz010);
//      Matrix23Row[j] += Mult * val;
//      val  = viscosity*(test100*ansatz001);
//      Matrix31Row[j] += Mult * val;
//      val  = viscosity*(test010*ansatz001);
//      Matrix32Row[j] += Mult * val;
//      val  = viscosity*(test100*ansatz100+test010*ansatz010
//                   +2*test001*ansatz001);
//      val += val1;
//      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
//    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];
      val1 = Mult*ansatz00;

      val = -val1*test10;
      MatrixRow1[j] += val;
      val = -val1*test01;
      MatrixRow2[j] += val;
//      val = -val1*test001;
//      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
//    MatrixRow3 = MatrixB3[i];

    test00 = Orig3[i];
    val1 = Mult*test00;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
//      ansatz001 = Orig2[j];

      val = -val1*ansatz10;
      MatrixRow1[j] += val;

      val = -val1*ansatz01;
      MatrixRow2[j] += val;

//      val = -val1*ansatz001;
//      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
//     Matrix33Row = Matrix_tilde_G33[i];
     test10 = Orig0[i];
     test01 = Orig1[i];
//     test001 = Orig2[i];
     for(j=0;j<N_L;j++)
     {
        ansatz00 = Orig4[j];
        val =  Mult * 2*mu * ansatz00;
        Matrix11Row[j] -= val * test10;
        Matrix22Row[j] -= val * test01;
//        Matrix33Row[j] -= val * test001;
     }
  }

  for(i=0;i<N_L;i++)
  {
     Matrix11Row = Matrix_G11[i];
     Matrix22Row = Matrix_G22[i];
//     Matrix33Row = Matrix_G33[i];
     test00 = Orig4[i];
     val =  Mult * test00;

     for(j=0;j<N_U;j++)
     {
        ansatz10 = Orig0[j];
        ansatz01 = Orig1[j];
//        ansatz001 = Orig2[j];

        Matrix11Row[j] -= val * ansatz10;
        Matrix22Row[j] -= val * ansatz01;
//        Matrix33Row[j] -= val * ansatz001;
     }
  }

  for(i=0;i<N_L;i++)
  {
     test00 = Orig4[i];
     MatrixRow1 = MatrixL[i];
     for(j=0;j<N_L;j++)
     {
        ansatz00 = Orig4[j];
        MatrixRow1[j] += Mult * test00 * ansatz00;
     }
  }
}

//================================================================================
/* ======================================================================
 * Type 3, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
 * Type 4, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
 * ======================================================================*/
void TimeNSType3_4NLVMS_ProjectionDD2D(double Mult, double *coeff,
                double *param, double hK,
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21;//**MatrixA13,
  double **MatrixA22 ;//**MatrixA23, **MatrixA31, **MatrixA32;
//  double **MatrixA33;
  double val1,val2, val3, valu1, valu2;//, valu3;val4,
  double **Matrix_tilde_G11, **Matrix_tilde_G22;//, **Matrix_tilde_G33;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row;//*Matrix13Row,
  double *Matrix22Row;//, *Matrix23Row, *Matrix31Row, *Matrix32Row;
//  double *Matrix33Row;
  double ansatz10, ansatz01, ansatz00;//, ansatz001;  // double ansatz000;
  double test00, test10, test01;//, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;//, *Orig4;
  int i,j,N_U, N_L;  // int N_P;
  double c0, viscosity;
  double u1, u2, mu;//u3

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
//  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
//  MatrixA23 = LocMatrices[5];
//  MatrixA31 = LocMatrices[6];
//  MatrixA32 = LocMatrices[7];
//  MatrixA33 = LocMatrices[8];
  Matrix_tilde_G11  = LocMatrices[4];
  Matrix_tilde_G22  = LocMatrices[5];
//  Matrix_tilde_G33  = LocMatrices[11];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
//  Orig2 = OrigValues[2]; // u_z
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // l

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
//  u3 = param[2]; // u3old

//  double *x = &param[12];
//  double *y = &param[13];
//  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[2];
  double *uConv = &param[0];
//  double *projection_space_label = &param[15];

  mu = turbulentviscosity(hK, gradu, u, uConv);
  viscosity = Mult*(mu/2.0+c0);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
//    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
//    Matrix23Row = MatrixA23[i];
//    Matrix31Row = MatrixA31[i];
//    Matrix32Row = MatrixA32[i];
//    Matrix33Row = MatrixA33[i];
    test10 = viscosity*Orig0[i];
    test01 = viscosity*Orig1[i];
//    test001 = viscosity*Orig2[i];
    test00 = Mult*Orig2[i];
    valu1 = u1 * test00;
    valu2 = u2 * test00;
//    valu3 = u3 * test000;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
//      ansatz001 = Orig2[j];

      //val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val1 = valu1*ansatz10+valu2*ansatz01;//+valu3*ansatz001;

      val2 = test10*ansatz10;
      val3 = test01*ansatz01;
//      val4 = test001*ansatz001;
      val1 += val2+val3;//+val4;
      Matrix11Row[j] += val2+val1;
      Matrix12Row[j] += test01*ansatz10;
//      Matrix13Row[j] += test001*ansatz100;
      Matrix21Row[j] += test10*ansatz01;
      Matrix22Row[j] += val3+val1;
//      Matrix23Row[j] += test001*ansatz010;
//      Matrix31Row[j] += test100*ansatz001;
//      Matrix32Row[j] += test010*ansatz001;
//      Matrix33Row[j] += val4+val1;
    } // endfor j
  } // endfor i

  val2 = Mult * mu;
  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
//     Matrix33Row = Matrix_tilde_G33[i];
     test10 = Orig0[i];
     test01 = Orig1[i];
//     test001 = Orig2[i];

     for(j=0;j<N_L;j++)
     {
        ansatz00 = Orig3[j];
        val1 = val2 * ansatz00;
        Matrix11Row[j] -= val1 * test10;
        Matrix22Row[j] -= val1 * test01;
//        Matrix33Row[j] -= val1 * test001;
     }
  }
}
