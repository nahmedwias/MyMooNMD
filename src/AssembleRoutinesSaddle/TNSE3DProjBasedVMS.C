#include <Database.h>
#include <Convolution.h>
#include <MooNMD_Io.h>
#include <ConvDiff.h>
#include <stdlib.h>
#include <TNSE3D_Routines.h>
#include <MainUtilities.h>

void TimeNSType4VMS_ProjectionDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double **MatrixL, **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33;
  double **Matrix_G11, **Matrix_G22,   **Matrix_G33;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5;
  int i,j,N_U, N_P, N_L;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, viscosity;
  
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixL   = LocMatrices[12];
  MatrixB1  = LocMatrices[13];
  MatrixB2  = LocMatrices[14];
  MatrixB3  = LocMatrices[15];
  MatrixB1T = LocMatrices[16];
  MatrixB2T = LocMatrices[17];
  MatrixB3T = LocMatrices[18];
  Matrix_tilde_G11  = LocMatrices[19];
  Matrix_tilde_G22  = LocMatrices[20];
  Matrix_tilde_G33  = LocMatrices[21];
  Matrix_G11  = LocMatrices[22];
  Matrix_G22  = LocMatrices[23];
  Matrix_G33  = LocMatrices[24];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p
  Orig5 = OrigValues[5]; // l
  
  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old


  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  double *projection_space_label = &param[15];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, projection_space_label[0]);
  
  mu = mu/2.0;
  viscosity = c0+mu;
  
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      
      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += val1;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];
      val1 = Mult*ansatz000;
      
      val = -val1*test100;
      MatrixRow1[j] += val;
      val = -val1*test010;
      MatrixRow2[j] += val;
      val = -val1*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];
    val1 = Mult*test000;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixRow1[j] += val;

      val = -val1*ansatz010;
      MatrixRow2[j] += val;

      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
     Matrix33Row = Matrix_tilde_G33[i];
     test100 = Orig0[i];
     test010 = Orig1[i];
     test001 = Orig2[i];
     for(j=0;j<N_L;j++)
     {       
        ansatz000 = Orig5[j];
        val =  Mult * 2*mu * ansatz000;
        Matrix11Row[j] -= val * test100;
        Matrix22Row[j] -= val * test010;
        Matrix33Row[j] -= val * test001;
     }
  }   

  for(i=0;i<N_L;i++)
  {
     Matrix11Row = Matrix_G11[i];
     Matrix22Row = Matrix_G22[i];
     Matrix33Row = Matrix_G33[i];
     test000 = Orig5[i];     
     val =  Mult * test000;

     for(j=0;j<N_U;j++)
     {        
        ansatz100 = Orig0[j];
        ansatz010 = Orig1[j];
        ansatz001 = Orig2[j];

        Matrix11Row[j] -= val * ansatz100;
        Matrix22Row[j] -= val * ansatz010;
        Matrix33Row[j] -= val * ansatz001;
     }
  }   

  for(i=0;i<N_L;i++)
  {
     test000 = Orig5[i];     
     MatrixRow1 = MatrixL[i];
     for(j=0;j<N_L;j++)
     {
        ansatz000 = Orig5[j];
        MatrixRow1[j] += Mult * test000 * ansatz000;
     }
  }
}

// ======================================================================
// Type 3, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
// Type 4, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLVMS_ProjectionDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val1,val2, val3, val4, valu1, valu2, valu3;
  double **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33; 
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;  // double ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_L;  // int N_P;
  double c0, viscosity;
  double u1, u2, u3, mu;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  Matrix_tilde_G11  = LocMatrices[9];
  Matrix_tilde_G22  = LocMatrices[10];
  Matrix_tilde_G33  = LocMatrices[11];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // l

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  double *projection_space_label = &param[15];
  
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, projection_space_label[0]);
  viscosity = Mult*(mu/2.0+c0);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    test100 = viscosity*Orig0[i];
    test010 = viscosity*Orig1[i];
    test001 = viscosity*Orig2[i];
    test000 = Mult*Orig3[i];
    valu1 = u1 * test000;
    valu2 = u2 * test000;
    valu3 = u3 * test000;

    for(j=0;j<N_U;j++)
    {      
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val1 = valu1*ansatz100+valu2*ansatz010+valu3*ansatz001;
      
      val2 = test100*ansatz100;
      val3 = test010*ansatz010;
      val4 = test001*ansatz001;
      val1 += val2+val3+val4;
      Matrix11Row[j] += val2+val1;
      Matrix12Row[j] += test010*ansatz100; 
      Matrix13Row[j] += test001*ansatz100;
      Matrix21Row[j] += test100*ansatz010;
      Matrix22Row[j] += val3+val1;
      Matrix23Row[j] += test001*ansatz010;
      Matrix31Row[j] += test100*ansatz001;
      Matrix32Row[j] += test010*ansatz001;
      Matrix33Row[j] += val4+val1;
    } // endfor j
  } // endfor i

  val2 = Mult * mu;
  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
     Matrix33Row = Matrix_tilde_G33[i];
     test100 = Orig0[i];
     test010 = Orig1[i];
     test001 = Orig2[i];

     for(j=0;j<N_L;j++)
     {       
	 //ansatz000 = Orig4[j];
        val1 = val2 * Orig4[j];
        Matrix11Row[j] -= val1 * test100;
        Matrix22Row[j] -= val1 * test010;
        Matrix33Row[j] -= val1 * test001;
     }
  }   
}

// ========================================================================
// parameters: u1old, u2old, G^H
// used for : projection-based VMS
// ========================================================================
void TimeNSParamsVelo_GradVelo_LargeScale3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3

  out[12] = in[0]; // x - coordinate for van Driest damping
  out[13] = in[1]; // y - coordinate for van Driest damping
  out[14] = in[2]; // z - coordinate for van Driest damping
  
  out[15] = in[15]; // projection space label
}