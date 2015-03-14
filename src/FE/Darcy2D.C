#include <Convolution.h>
#include <Database.h>
#include <TNSE2D_Routines.h>
#include <stdlib.h>


// ======================================================================
// Type 1, Standard Galerkin with RT elements
// Note: this assembling is not used anymore in Darcy_2D_DirectSolver
// ======================================================================
void DarcyType1Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB;
  double *Rhs0, *Rhs1;
  double val;
  double *MatrixRow, *MatrixRow1;

  double ansatz_x_00, ansatz_y_00;
  double test00, test_x_00, test_y_00;
  double ansatz_div;
  double ansatz_x_10, ansatz_y_01;

  double q00,q01,q10,p00,p10,p01,u00,u01,u10,v00,v10,v01;

  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0;

  MatrixA = LocMatrices[0];
  MatrixB = LocMatrices[1];

  Rhs0 = LocRhs[0];
  Rhs1 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];                          // u
  Orig1 = OrigValues[1];                          // u_x
  Orig2 = OrigValues[2];                          // u_y
  Orig3 = OrigValues[3];                          // p

  c0 = coeff[0];                                  // sigma

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];

    // HERE CHECK WHETHER (ROW-)SIGN SHOULD BE INVERSED
    int newsign_i;
    if (N_U == 4)
    {
      newsign_i = TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[i];

    }
    else if (N_U == 3)
    {
      newsign_i = TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[i];
    }
    test_x_00 = newsign_i*Orig0[i];
    test_y_00 = newsign_i*Orig0[N_U+i];

    for(j=0;j<N_U;j++)
    {
      // HERE CHECK WHETHER (COLUMN-)SIGN SHOULD BE INVERSED
      int newsign_j;
      if (N_U == 4)
      {
        newsign_j = TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[j];
      }
      else if (N_U == 3)
      {
        newsign_j = TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[j];
      }
      ansatz_x_00 = newsign_j*Orig0[j];
      ansatz_y_00 = newsign_j*Orig0[N_U+j];

      val  = c0*(test_x_00*ansatz_x_00 + test_y_00*ansatz_y_00);
      MatrixRow[j] += Mult * val;

    }                                             // endfor j
  }                                               // endfor i

  for(i=0;i<N_P;i++)
  {

    MatrixRow1 = MatrixB[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {

      int newsign_j;
      if (N_U == 4)
      {
        newsign_j = TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[j];
      }
      else if (N_U == 3)
      {
        newsign_j = TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[j];
      }
      ansatz_x_10 = newsign_j*Orig1[j];
      ansatz_y_01 = newsign_j*Orig2[N_U+j];

      ansatz_div = ansatz_x_10 + ansatz_y_01;

      val = -Mult*test00*ansatz_div;
      MatrixRow1[j] += val;
    }                                             // endfor j

  }                                               // endfor i
}

/* ======================================================================
   depending on the global orientation of the normals at inner edges the sign 
   of the basis functions need to be inversed. This is only important for 
   Raviart-Thomas finite elements. (called from DarcyRaviartThomas(...))
   
   This information is stored in the FE-Descriptor. However the FE-Descriptor 
   is not accessible in the assembling routine DarcyRaviartThomas. Therefore we
   use global parameters in TDatabase. These need to be set for every cell.
   
   Here it is assumed that Raviart-Thomas elements of order 0,1,2, or 3 on
   triangles or quadrilaterals are used
*/
int GetSignOfThisDOF(int N_DOF, int DOF)
{
  switch (N_DOF)
  {
  case 3:// Raviart-Thomas zeroth order, triangles
    return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[DOF];
    break;
    // note that RT2 on quadrilaterals and RT3 on triangles have 24 basis functions
  case 4:// Raviart-Thomas zeroth order, quadrilaterals
    return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[DOF];
    break;
  case 6: // BDM first order, triangles
    // degree of freedom on an edge, no inner degree of freedom
    return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%2)/2];
    break;
  case 8:
    // Raviart-Thomas first order, triangles
    if(TDatabase::ParamDB->VELOCITY_SPACE == 1001)
    {
      if(DOF<6) // degree of freedom on an edge
        return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%2)/2];
      else // inner degree of freedom
        return 1;
    }
    else if (TDatabase::ParamDB->VELOCITY_SPACE == 1011) //BDM first order, quadrilaterals
    {
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%2)/2];
    }
    break;
  case 12:
    // Raviart-Thomas first order, quadrilaterals
    if(TDatabase::ParamDB->VELOCITY_SPACE == 1001)
    {
    if(DOF<8) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%2)/2];
    else // inner degree of freedom
      return 1;
    }
    else if (TDatabase::ParamDB->VELOCITY_SPACE == 1012) //BDM second order, triangles
    {
      if(DOF<9) // degree of freedom on an edge
        return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%3)/3];
      else // inner degree of freedom
        return 1;
    }
    break;
  case 14://BDM second order, quadrilaterals
    if(DOF<12) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%3)/3];
    else // inner degree of freedom
      return 1;
    break;
  case 15:// Raviart-Thomas second order, triangles
    if(DOF<9) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%3)/3];
    else // inner degree of freedom
      return 1;
    break;
  case 20:// BDM third order. triangles
    if(DOF<12) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%4)/4];
    else // inner degree of freedom
      return 1;
    break;
  case 22://BDM third order, quadrilaterals
    if(DOF<16) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%4)/4];
    else // inner degree of freedom
      return 1;
    break;
  case 24:
    if(TDatabase::ParamDB->VELOCITY_SPACE == 1002)
    {
      // Raviart-Thomas second order, quadrilaterals
      if(DOF<12) // degree of freedom on an edge
        return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%3)/3];
      else // inner degree of freedom
        return 1;
    }
    else if(TDatabase::ParamDB->VELOCITY_SPACE == 1003)
    {
      // Raviart-Thomas third order, triangles
      if(DOF<12) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%4)/4];
      else // inner degree of freedom
        return 1;
    }
    else
    {
      Error("VELOCITY_SPACE has to be set to either 1002 or 1003\n");
      exit(0);
    }
    break;
  case 40:// Raviart-Thomas third order, quadrilaterals
    if(DOF<16) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%4)/4];
    else
      return 1;
    break;
  default:
       OutPut("N_DOF" << N_DOF << endl);
       OutPut("WARNING: Unknown Raviart-Thomas or BDM element !" << endl);
       return 1;
       break;
  }
}

// ======================================================================
// (DarcyType 1)
// Standard Galerkin with Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM)
// elements
// ======================================================================
void DarcyRaviartThomas(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB;
  double *Rhs0, *Rhs1;
  double val;
  double *MatrixRow, *MatrixRow1;

  double ansatz_x_00, ansatz_y_00;
  double test00, test_x_00, test_y_00, f00;
  double ansatz_div;
  double ansatz_x_10, ansatz_y_01;

  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0,c1,c2,c3;

  MatrixA = LocMatrices[0];
  MatrixB = LocMatrices[1];

  Rhs0 = LocRhs[0];
  Rhs1 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];                          // u
  Orig1 = OrigValues[1];                          // u_x
  Orig2 = OrigValues[2];                          // u_y
  Orig3 = OrigValues[3];                          // p

  c0 = coeff[0];                                  // sigma
  c1 = coeff[1];                                  // f1
  c2 = coeff[2];                                  // f2
  c3 = coeff[3];                                  // g(x,y)

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];

    // HERE CHECK WHETHER (ROW-)SIGN SHOULD BE INVERSED
    int newsign_i = GetSignOfThisDOF(N_U,i);
    test_x_00 = newsign_i*Orig0[i];
    test_y_00 = newsign_i*Orig0[N_U+i];
    
    Rhs0[i] += Mult*(c1*test_x_00 + c2*test_y_00);
    //OutPut("rhs i " << i << " c1 " << c1 << " c2 " << c2 
    //       << "   test_x " << test_x_00 << "\ttest_y " << test_y_00
    //       << endl);
    for(j=0;j<N_U;j++)
    {
      // HERE CHECK WHETHER (COLUMN-)SIGN SHOULD BE INVERSED
      int newsign_j = GetSignOfThisDOF(N_U,j);
      ansatz_x_00 = newsign_j*Orig0[j];
      ansatz_y_00 = newsign_j*Orig0[N_U+j];

      // A: u_x v_x + u_y v_y

      val  = c0*(test_x_00*ansatz_x_00 + test_y_00*ansatz_y_00);
      MatrixRow[j] += Mult * val;
      
      //OutPut("N_U " << N_U << ", i " << i << ", j " << j << ", Entry " 
      //       << MatrixRow[j] <<endl);
    }                                             // endfor j
  }                                               // endfor i

  for(i=0;i<N_P;i++)
  {

    MatrixRow1 = MatrixB[i];

    test00 = Orig3[i];

    // assemble rhs: div u = g
    // rhs: -(g,q)
    // c3 = g(x,y)
    Rhs1[i] += -Mult*test00*c3;

    //OutPut("PressureDOF " << i << ", PressureVal " << test00 << ",\t g " << c3 
    //       << endl);
    for(j=0;j<N_U;j++)
    {
      int newsign_j = GetSignOfThisDOF(N_U,j);
      ansatz_x_10 = newsign_j*Orig1[j];
      ansatz_y_01 = newsign_j*Orig2[N_U+j];

      /*cout  << " N_U = " << N_U << " , newsign_j  " << newsign_j << endl;
      cout << "dU/dx base(" << j << ") =  " << Orig1[j] << endl
      << " dV/dy base(" << j << ") = " << Orig2[N_U+j] << endl;*/
      ansatz_div = ansatz_x_10 + ansatz_y_01;
      /*cout << " div: " << ansatz_div << endl;
      cout << " ** "  << endl; */

      val = -Mult*test00*ansatz_div;
      MatrixRow1[j] += val;

    }                                             // endfor j

  }                                               // endfor i
}


// ======================================================================
// (DarcyType 2)  CGLS Formulation [cfr. Correa,Loula (2007)]
// ======================================================================
/**
   sigma(u,v)-(p,divu)-(q,divu)
    -1/2* (1/sigma) * (sigma u + grad(p),sigma v + grad(q))
    +1/2*sigma (div(u),div(v))
    +1/2*(1/sigma)*(curl(u),curl(v))

    curl(u) = dx(u2)-dy(u1)

    A = sigma(u,v)
        -1/2* (1/sigma) * (sigma u,sigma v)
  +1/2*sigma (div(u),div(v))
+1/2*(1/sigma)*(curl(u),curl(v)

BT = -(p,divu)
-1/2* (1/sigma) (grad(p),sigma v)

B = -(q,divv)
-1/2* (1/sigma) (grad(q),sigma u)

C = -1/2* (1/sigma) * (grad(p),grad(q))

**/
void DarcyP1P1CurlStab(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T, **MatrixC;
  double *Rhs0, *Rhs1, *Rhs2;
  double val;
  double *MatrixRow1, *MatrixRow2, *MatrixRowC;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;

  double test00,test01,test10;
  double ansatz00, ansatz10, ansatz01;
  //double u00,u01,u10,p00,p10,p01;
  double q00,q01,q10,p00,p10,p01,u00,u01,u10,v00,v10,v01;

  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5;
  int i,j, N_U, N_P;
  double c0,c3;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixC   = LocMatrices[4];
  MatrixB1  = LocMatrices[5];
  MatrixB2  = LocMatrices[6];
  MatrixB1T  = LocMatrices[7];
  MatrixB2T  = LocMatrices[8];

  Rhs0 = LocRhs[0];
  Rhs1 = LocRhs[1];
  Rhs2 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];                          // u
  Orig1 = OrigValues[1];                          // u_x
  Orig2 = OrigValues[2];                          // u_y
  Orig3 = OrigValues[3];                          // p
  Orig4 = OrigValues[4];                          // p_x
  Orig5 = OrigValues[5];                          // p_y

  c0 = coeff[0];                                  // sigma
  c3 = coeff[3];                                  // g(x,y)

  // assembling for velocity test functions
  for(i=0;i<N_U;i++)
  {
    v00 = Orig0[i];
    v10 = Orig1[i];
    v01 = Orig2[i];

    // rhs: sigma/2*(g,div v)
    // quad_weigth * test_function * g
    // c3 = g(x,y)
    Rhs0[i] += Mult*c0/2*c3*v10;
    Rhs1[i] += Mult*c0/2*c3*v01;

    /*
    A = sigma(u,v)
        -1/2* (1/sigma) * (sigma u,sigma v)
    +1/2*sigma (div(u),div(v))
    +1/2*(1/sigma) (curl(u),curl(v))
    */
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    for(j=0;j<N_U;j++)
    {
      u00 = Orig0[j];
      u10 = Orig1[j];
      u01 = Orig2[j];
      val  = c0*v00*u00;

      // stabilization terms
      val = val - c0/2*v00*u00;
      val = val + c0/2*v10*u10;                   // div-div term
      val = val + 1./c0/2*v01*u01;                // curl-curl term

      Matrix11Row[j] += Mult * val;

      val  = c0*v00*u00;
      // stabilization terms
      val = val - c0/2*v00*u00;
      val = val + c0/2*v01*u01;                   // div-div term
      val = val + 1./c0/2*v10*u10;                // curl-curl term

      Matrix22Row[j] += Mult * val;

      // div-div term
      val = c0/2*v10*u01;
      val = val - 1./c0/2*v01*u10;                // curl-curl term
      Matrix12Row[j] += Mult * val;
      val = c0/2*v01*u10;
      val = val - (1./c0)/2*v10*u01;              // curl-curl term
      Matrix21Row[j] += Mult * val;

    }                                             // endfor j

    /*
      BT = -(p,divu)
     -1/2* (grad(p),v)
    */
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {

      p00 = Orig3[j];
      p10 = Orig4[j];
      p01 = Orig5[j];

      // p div(v)
      val  = -p00 * v10;                          // p dx(v_x)
      val  = val - 0.5*(p10 * v00);               // dx(p)*v_x
      MatrixRow1[j] += Mult*val;

      val  = -p00 * v01;                          // p dy(v_y)
      val  = val - 0.5*(p01 * v00);               // dy(p)*v_y
      MatrixRow2[j] += Mult*val;

    }                                             // endfor j

  }                                               // endfor i

  // assembling for pressure test functions
  for(i=0;i<N_P;i++)
  {
    q00 = Orig3[i];
    q10 = Orig4[i];
    q01 = Orig5[i];

    // assemble rhs: (g,q)
    // quad_weigth * test_function * g
    // c3 = g(x,y)
    Rhs2[i] += -Mult*q00*c3;

    // pressure-pressure block
    /*
      C = -1/2* (1/sigma) * (grad(p),grad(q))
    */
    MatrixRowC = MatrixC[i];
    for(j=0;j<N_P;j++)
    {
      p10 = Orig4[j];
      p01 = Orig5[j];
      val = - 0.5*(1./c0) * (p10*q10+p01*q01);
      MatrixRowC[j] += Mult*val;
    }

    /*
      B = -(q,divu)
    -1/2* (1/sigma) (grad(q),sigma u)
    */
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    for(j=0;j<N_U;j++)
    {

      u00 = Orig0[j];
      u10 = Orig1[j];
      u01 = Orig2[j];

      val = -q00*u10;                             // q dx(u)
      val = val - 0.5*(q10*u00);                  //dx(q)*u
      MatrixRow1[j] += Mult*val;

      val = -q00*u01;                             // q dy(u)
      val = val - 0.5*(q01*u00);                  //dy(q)*u
      MatrixRow2[j] += Mult*val;

    }                                             // endfor j

  }                                               // endfor i
}
