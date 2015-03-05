// ======================================================================
//  TCD2D.C 
// ======================================================================

#include <Database.h>
#include <ConvDiff2D_Routines.h>
#include <MooNMD_Io.h>
#include <math.h>
#include <stdlib.h>

double ComputeAlpha(double hK)
{
  double alpha;
  
  alpha = TDatabase::ParamDB->ARTIFICIAL_VISCOSITY_CONSTANT*
     pow(hK,TDatabase::ParamDB->ARTIFICIAL_VISCOSITY_POWER);
  return(alpha);

  // this is just for the response to the referee and the special example 
  // in [JKL05]
  double b, eps, Pe, t;

  b = sqrt(5.0);
  eps = 1/TDatabase::ParamDB->RE_NR;
  Pe = b*hK/(2*eps);
  t = 1/tanh(Pe) - 1/Pe;
  alpha = t*hK/(2*b);
  return(alpha);
}

void MatrixMRhsAssemble(double Mult, double *coeff, double *param,
                           double hK, 
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, *MatrixRow;
  double ansatz00;
  double test00;
  double *Orig0;
  int i,j, N_;
  double c4; 

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];

  c4 = coeff[4]; // f

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test00 = Orig0[i];

    Rhs[i] += Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig0[j];

      MatrixRow[j] += Mult*ansatz00*test00;
    } // endfor j
  } // endfor i
}


void MatrixMRhsAssemble_Axial3D(double Mult, double *coeff, double *param,
                           double hK, 
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, *MatrixRow;
  double ansatz00, x, r;
  double test00;
  double *Orig0;
  int i,j, N_;
  double c4; 

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];

  c4 = coeff[4]; // f
  x  = coeff[20]; // see DiscreteForm2D.C
  r  = fabs(x);
  
  if(r<1e-12)
   {
   OutPut("check MatrixMRhsAssemble_Axial3D x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as positive points"<<endl);
   }
   
//    cout<< "r " <<  r<< "  c4  " << c4 << endl;
//     exit(0);

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test00 = Orig0[i];

    Rhs[i] += r*Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig0[j];

      MatrixRow[j] += r*Mult*ansatz00*test00;
    } // endfor j
  } // endfor i
}


void MatrixMRhsAssemble_SUPG_Axial3D(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **Matrix, **MatrixS, *Rhs, val, *MatrixRow, *MatrixSRow;
  double ansatz00, x, r;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c00, c11, c22, c33; 
  double tau, bgradv, bb;
  double theta1 = TDatabase::TimeDB->THETA1;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  Matrix = LocMatrices[0];
  MatrixS= LocMatrices[1];
  
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // c
  c4 = coeff[4]; // f
  
  x  = coeff[20]; // see DiscreteForm2D.C
  r  = fabs(x);
  
  if(r<1e-12)
   {
   OutPut("check MatrixMRhsAssemble_SUPG_Axial3D x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as positive points"<<endl);
   }
   
   cout<< "r " <<  r<< "  c4  " << c4 << endl;
//     exit(0);  
  
  c00 = theta1 * time_step * c0;
  c11 = theta1 * time_step * c1;
  c22 = theta1 * time_step * c2;
  // reactive coefficient, inclusive term from the temporal derivative
  c33 = 1.0 + theta1 * time_step * c3;
  if (TDatabase::ParamDB->SDFEM_TYPE==8)
  {
      c33 = theta1 * time_step * c3;
  }
  if (fabs(c11) > fabs(c22))
      bb = fabs(c11);
  else
      bb = fabs(c22);
  // this is \tilde tau
  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, bb);
  // scale appropriately
  //OutPut(tau << " ");
  // do not apply for paper with J. Novo
  if ((TDatabase::ParamDB->SDFEM_TYPE!=9)&&(TDatabase::ParamDB->SDFEM_TYPE!=10)&&(TDatabase::ParamDB->SDFEM_TYPE!=11))
    tau *= theta1 * time_step;
  //OutPut(theta1 << " " << tau << " : ");
  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    MatrixSRow = MatrixS[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;

    Rhs[i] += r*Mult*(test00+tau*bgradv)*c4;

    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig2[j];

      MatrixRow[j] += r*Mult * ansatz00*test00;
      MatrixSRow[j] += r*Mult * ansatz00*bgradv;    
      
    } // endfor j
  } // endfor i
}

void MatricesAKRhsAssemble_SUPG_Axial3D(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixK, *Rhs, *MatrixRowA, *MatrixRowK;
  double **MatrixS, *MatrixRowS, ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, r, x;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  double val, val1, val2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5, Pe; 
  double c00, c11, c22, c33;
  double tau, bgradv, bb, res, sigma, norm_b;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  MatrixA = LocMatrices[0];
  MatrixK= LocMatrices[1];
  
  if (TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)  
   MatrixS= LocMatrices[2];
  
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  // coefficients of the problem
  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // c
  c4 = coeff[4]; // f

  x  = coeff[20]; // see DiscreteForm2D.C
  r  = fabs(x);
  
  if(r<1e-12)
   {
   OutPut("check MatrixMRhsAssemble_SUPG_Axial3D x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as positive points"<<endl);
   }
   
  // coefficients for the stabilization parameter
  val = theta1 * time_step;
  c00 = val * c0;
  c11 = val * c1;
  c22 = val * c2;
  // reactive coefficient, inclusive term from the temporal derivative
  c33 = 1.0 + val * c3;
  if (TDatabase::ParamDB->SDFEM_TYPE==8)
  {
      c33 = val * c3;
  }
  if (fabs(c11) > fabs(c22))
      bb = fabs(c11);
  else
      bb = fabs(c22);
  // this is tau
  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, bb);

  if (TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)
  {
      // rhs from previous time step 
      c5 = coeff[5];   
      // compute residual
      res = param[0] + theta1*time_step*(-c0*(param[3]+param[4]) + c1*param[1]
					 +c2*param[2] + c3*param[0])
	  -param[5] +  theta2*time_step*(-c0*(param[8]+param[9]) + c1*param[6]
					 +c2*param[7] + c3*param[5])
	  -theta3*time_step*c5 - theta4*time_step*c4;
      /*c00 =  time_step * theta1 * c0;
      c11 =  time_step * theta1 * c1;
      c22 =  time_step * theta1 * c2;
      c33 = 1.0 + time_step * theta1 *c3;*/
      c5 = time_step * theta4 * c4;
      // compute the parameter, c5 is just a dummy
      sigma = Compute_SOLD_sigma(hK, c00, c11, c22, c33, c5, bb, tau, param, res, 1,1);
      //OutPut( param[0] << " " << param[5] <<  " " << res << " " <<  sigma << endl);
      val2 = Mult * sigma;
      if (TDatabase::ParamDB->SOLD_TYPE==2)
      {
	  norm_b = c1*c1 + c2*c2;
	  if (norm_b >1e-10)
	      val2 /= norm_b;
	  else
	      val2 = 0.0;
      }
  }
  // scale appropriately, after it is used for the SOLD scheme
  // do not apply for paper with J. Novo
  // this is \tilde tau
  if ((TDatabase::ParamDB->SDFEM_TYPE!=9)&&(TDatabase::ParamDB->SDFEM_TYPE!=10)&&(TDatabase::ParamDB->SDFEM_TYPE!=11))
    tau *= val;

  // loop over the basis functions
  for(i=0;i<N_;i++)
  {
    MatrixRowA = MatrixA[i];
    MatrixRowK = MatrixK[i];
    if (TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)
	MatrixRowS = MatrixS[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;
    // scaling with the stabilization parameter
    // scaling with the time step is done in the main program
    bgradv *= tau;
    // THIS CHANGEs TEST00 !
    test00 += bgradv;
    Rhs[i] += r*Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];
	
      // Galerkin part of the bilinear form
      val1 = c1*ansatz10+c2*ansatz01;
      val1+= c3*ansatz00;

      val = c0*(test10*ansatz10+test01*ansatz01);
      // val1*test00 includes the SUPG part of the convective and reactive term
      val += val1*test00;
      // diffusion part of the SUPG stabilization
      val -= c0*(ansatz20 + ansatz02) * bgradv;
      
      MatrixRowA[j] += r*Mult * val;
      // time derivative part of the SUPG stabilization          
      MatrixRowK[j] += r*Mult * ansatz00*bgradv;
      
      // isotropic SOLD method
      if ((TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE) && (TDatabase::ParamDB->SOLD_TYPE==1))
      {
	  MatrixRowS[j] += val2 * (test10*ansatz10+test01*ansatz01);
      }
      if ((TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE) && (TDatabase::ParamDB->SOLD_TYPE==2))
      {

	  MatrixRowS[j] += val2 * (-c2*ansatz10+c1*ansatz01)*(-c2*test10+c1*test01);
      }
     
    } // endfor j
  } // endfor i
}


void MatrixMRhsAssemble_SUPG(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **Matrix, **MatrixS, *Rhs, val, *MatrixRow, *MatrixSRow;
  double ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c00, c11, c22, c33; 
  double tau, bgradv, bb;
  double theta1 = TDatabase::TimeDB->THETA1;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  Matrix = LocMatrices[0];
  MatrixS = LocMatrices[0];
    
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // c
  c4 = coeff[4]; // f
  
  c00 = theta1 * time_step * c0;
  c11 = theta1 * time_step * c1;
  c22 = theta1 * time_step * c2;
  // reactive coefficient, inclusive term from the temporal derivative
  c33 = 1.0 + theta1 * time_step * c3;
  if (TDatabase::ParamDB->SDFEM_TYPE==8)
  {
      c33 = theta1 * time_step * c3;
  }
  if (fabs(c11) > fabs(c22))
      bb = fabs(c11);
  else
      bb = fabs(c22);
  // this is \tilde tau
  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, bb);
  // scale appropriately
  //OutPut(tau << " ");
  // do not apply for paper with J. Novo
  if ((TDatabase::ParamDB->SDFEM_TYPE!=9)&&(TDatabase::ParamDB->SDFEM_TYPE!=10)&&(TDatabase::ParamDB->SDFEM_TYPE!=11))
    tau *= theta1 * time_step;
  //OutPut(theta1 << " " << tau << " : ");
  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    MatrixSRow = MatrixS[i];   
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;

    Rhs[i] += Mult*(test00+tau*bgradv)*c4;

    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig2[j];

      MatrixRow[j] += Mult * ansatz00*test00;
      MatrixSRow[j] += Mult * ansatz00*bgradv;     
    } // endfor j
  } // endfor i
}

void MatricesAKRhsAssemble_SUPG(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixK, *Rhs, *MatrixRowA, *MatrixRowK;
  double **MatrixS, *MatrixRowS, ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  double val, val1, val2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5, Pe; 
  double c00, c11, c22, c33;
  double tau, bgradv, bb, res, sigma, norm_b;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  MatrixA = LocMatrices[0];
  MatrixK= LocMatrices[1];
  
  if (TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)  
   MatrixS= LocMatrices[2];
  
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  // coefficients of the problem
  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // c
  c4 = coeff[4]; // f
  
  // coefficients for the stabilization parameter
  val = theta1 * time_step;
  c00 = val * c0;
  c11 = val * c1;
  c22 = val * c2;
  // reactive coefficient, inclusive term from the temporal derivative
  c33 = 1.0 + val * c3;
  if (TDatabase::ParamDB->SDFEM_TYPE==8)
  {
      c33 = val * c3;
  }
  if (fabs(c11) > fabs(c22))
      bb = fabs(c11);
  else
      bb = fabs(c22);
  // this is tau
  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, bb);

  if (TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)
  {
      // rhs from previous time step 
      c5 = coeff[5];   
      // compute residual
      res = param[0] + theta1*time_step*(-c0*(param[3]+param[4]) + c1*param[1]
					 +c2*param[2] + c3*param[0])
	  -param[5] +  theta2*time_step*(-c0*(param[8]+param[9]) + c1*param[6]
					 +c2*param[7] + c3*param[5])
	  -theta3*time_step*c5 - theta4*time_step*c4;
      /*c00 =  time_step * theta1 * c0;
      c11 =  time_step * theta1 * c1;
      c22 =  time_step * theta1 * c2;
      c33 = 1.0 + time_step * theta1 *c3;*/
      c5 = time_step * theta4 * c4;
      // compute the parameter, c5 is just a dummy
      sigma = Compute_SOLD_sigma(hK, c00, c11, c22, c33, c5, bb, tau, param, res, 1,1);
      //OutPut( param[0] << " " << param[5] <<  " " << res << " " <<  sigma << endl);
      val2 = Mult * sigma;
      if (TDatabase::ParamDB->SOLD_TYPE==2)
      {
	  norm_b = c1*c1 + c2*c2;
	  if (norm_b >1e-10)
	      val2 /= norm_b;
	  else
	      val2 = 0.0;
      }
  }
  // scale appropriately, after it is used for the SOLD scheme
  // do not apply for paper with J. Novo
  // this is \tilde tau
  if ((TDatabase::ParamDB->SDFEM_TYPE!=9)&&(TDatabase::ParamDB->SDFEM_TYPE!=10)&&(TDatabase::ParamDB->SDFEM_TYPE!=11))
    tau *= val;

  // loop over the basis functions
  for(i=0;i<N_;i++)
  {
    MatrixRowA = MatrixA[i];
    MatrixRowK = MatrixK[i];
    if (TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)
	MatrixRowS = MatrixS[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;
    // scaling with the stabilization parameter
    // scaling with the time step is done in the main program
    bgradv *= tau;
    // THIS CHANGEs TEST00 !
    test00 += bgradv;
    Rhs[i] += Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];
	
      // Galerkin part of the bilinear form
      val1 = c1*ansatz10+c2*ansatz01;
      val1+= c3*ansatz00;

      val = c0*(test10*ansatz10+test01*ansatz01);
      // val1*test00 includes the SUPG part of the convective and reactive term
      val += val1*test00;
      // diffusion part of the SUPG stabilization
      val -= c0*(ansatz20 + ansatz02) * bgradv;
      
      MatrixRowA[j] += Mult * val;
      // time derivative part of the SUPG stabilization          
      MatrixRowK[j] += Mult * ansatz00*bgradv;
      
      // isotropic SOLD method
      if ((TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE) && (TDatabase::ParamDB->SOLD_TYPE==1))
      {
	  MatrixRowS[j] += val2 * (test10*ansatz10+test01*ansatz01);
      }
      if ((TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE) && (TDatabase::ParamDB->SOLD_TYPE==2))
      {

	  MatrixRowS[j] += val2 * (-c2*ansatz10+c1*ansatz01)*(-c2*test10+c1*test01);
      }
     
    } // endfor j
  } // endfor i
}

void MatrixARhsAssemble(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixK, *Rhs, val, *MatrixRowA, *MatrixRowK;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, Pe, h; 

  MatrixA = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // c
  c4 = coeff[4]; // f

  if ((TDatabase::ParamDB->DISCTYPE==5)||(TDatabase::ParamDB->DISCTYPE==6)
      ||(TDatabase::ParamDB->DISCTYPE==7))
  {
    h = ComputeAlpha(hK);
    c0+= h;  
  }
  for(i=0;i<N_;i++)
  {
    MatrixRowA = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs[i] += Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      MatrixRowA[j] += Mult * val;
                
    } // endfor j
  } // endfor i
}



void MatrixARhsAssemble_Axial3D(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixK, *Rhs, val, *MatrixRowA, *MatrixRowK;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, Pe, h, x, r; 

  MatrixA = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // c
  c4 = coeff[4]; // f

  if ((TDatabase::ParamDB->DISCTYPE==5)||(TDatabase::ParamDB->DISCTYPE==6)
      ||(TDatabase::ParamDB->DISCTYPE==7))
  {
    h = ComputeAlpha(hK);
    c0+= h;  
  }
  
  x  = coeff[20]; // see DiscreteForm2D.C
  r  = fabs(x);
  
  if(r<1e-12)
   {
   OutPut("check MatrixMRhsAssemble_Axial3D x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as positive points"<<endl);
   }
   
//   cout<< "r " <<  r<< "  c4  " << c4 << endl;
//     exit(0);  */   
  
  for(i=0;i<N_;i++)
  {
    MatrixRowA = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs[i] += r*Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      MatrixRowA[j] += r*Mult * val;
                
    } // endfor j
  } // endfor i
}

// ==

// ======================================================================
// assemble matrices B,C,M from 
// ( A B )
// ( C M ) 
// ======================================================================

void MatricesAssemble_VMM(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixMcoarse, **MatrixB, **MatrixC, val, *MatrixRowMcoarse;
  double *MatrixRowB, *MatrixRowC;
  double ansatz10, ansatz01 ;
  double test10, test01, h;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_, N_coarse;

  MatrixMcoarse = LocMatrices[0];
  MatrixB = LocMatrices[1];
  MatrixC = LocMatrices[2];

  N_ = N_BaseFuncts[0];
  N_coarse = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // phi_x
  Orig1 = OrigValues[1]; // phi_y
  Orig2 = OrigValues[2]; // psi_x
  Orig3 = OrigValues[3]; // psi_y

  h = ComputeAlpha(hK);

  for(i=0;i<N_coarse;i++)
  {
    MatrixRowMcoarse = MatrixMcoarse[i];
    MatrixRowC = MatrixC[i];

    test10 = Orig2[i];
    test01 = Orig3[i];

    for(j=0;j<N_coarse;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      val = ansatz10*test10 + ansatz01*test01;
      MatrixRowMcoarse[j] += Mult * val;
    } // endfor j
    
    for (j=0;j<N_;j++)
    {
       ansatz10 = Orig0[j];
       ansatz01 = Orig1[j];
       val = ansatz10*test10 + ansatz01*test01;
       MatrixRowC[j] -= Mult * val;       
    }
  } // endfor i
  
  for (i=0;i<N_;i++)
  {
     MatrixRowB = MatrixB[i];
     
     test10 = Orig0[i];
     test01 = Orig1[i];
     for(j=0;j<N_coarse;j++)
     {
        ansatz10 = Orig2[j];
        ansatz01 = Orig3[j];
        val = ansatz10*test10+ansatz01*test01;
        MatrixRowB[j] -=  Mult*val*h;
     }
  }
}

// ======================================================================
// assemble matrices B1, B2, C1, C2, M from 
// ( A     B1 B2 )
// ( C1 c2   M ) 
// ======================================================================

void MatricesAssemble_VMM_KL02(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixMcoarse, **MatrixB1, **MatrixC1, val, *MatrixRowMcoarse;
  double **MatrixB2, **MatrixC2;
  double *MatrixRowB1, *MatrixRowC1, *MatrixRowB2, *MatrixRowC2;
  double ansatz10, ansatz01 ;
  double test10, test01;
  double *Orig0, *Orig1, h;
  int i,j, N_, N_coarse;
 
  MatrixMcoarse = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixC1 = LocMatrices[3];
  MatrixC2 = LocMatrices[4];

  N_ = N_BaseFuncts[0];
  N_coarse = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // phi_x
  Orig1 = OrigValues[1]; // phi_y

  h = ComputeAlpha(hK);

  for(i=0;i<N_coarse;i++)
  {
    MatrixRowMcoarse = MatrixMcoarse[i];
    MatrixRowC1 = MatrixC1[i];
    MatrixRowC2 = MatrixC2[i];

    for(j=0;j<N_coarse;j++)
    {
      MatrixRowMcoarse[j] += Mult;
    } // endfor j
    
    for (j=0;j<N_;j++)
    {
       ansatz10 = Orig0[j];
       ansatz01 = Orig1[j];
       MatrixRowC1[j] -= Mult * ansatz10;       
       MatrixRowC2[j] -= Mult * ansatz01;       
    }
  } // endfor i
  
  for (i=0;i<N_;i++)
  {
     MatrixRowB1 = MatrixB1[i];
     MatrixRowB2 = MatrixB2[i];
     
     test10 = Orig0[i];
     test01 = Orig1[i];
     for(j=0;j<N_coarse;j++)
     {
        MatrixRowB1[j] -=  Mult*test10*h;
        MatrixRowB2[j] -=  Mult*test01*h;       
     } 
  }
}

// ======================================================================
// MATRICES FOR REACTION PART OF BULK PRECIPITATION
// ======================================================================

// ======================================================================
// assemble mass matriX 
// ======================================================================
void MatrixMAssemble_Bulk(double Mult, double *coeff, double *param,
                          double hK, 
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz00;
  double test00;
  double *Orig0;
  int i,j, N_;
 
  Matrix = LocMatrices[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test00 = Orig0[i];
 
    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig0[j];

      if (!TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING)
	  MatrixRow[j] += Mult * ansatz00*test00;
      else
      {
	  // add to diagonal
          //if (TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING==1)
	      MatrixRow[i] += Mult * ansatz00*test00;
      }
    } // endfor j
  } // endfor i
}

// ======================================================================
// assemble matrix and rhs for upwind discretization
// ======================================================================
void MatricesA_Assemble_Bulk(double Mult, double *coeff, double *param,
			    double hK, 
			    double **OrigValues, int *N_BaseFuncts,
			    double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, val, *MatrixRowA;
  double ansatz00, ansatz10, ansatz01 ;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, Pe; 
  double k;

  MatrixA = LocMatrices[0];
 
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // u_1
  c2 = coeff[2]; // u_2
  c3 = coeff[3]; // other concentration

  for(i=0;i<N_;i++)
  {
      MatrixRowA = MatrixA[i];
      test10 = Orig0[i];
      test01 = Orig1[i];
      test00 = Orig2[i];
		
      for(j=0;j<N_;j++)
      {
	  ansatz10 = Orig0[j];
	  ansatz01 = Orig1[j];
	  ansatz00 = Orig2[j];
	      
	  val =  c0*(test10*ansatz10+test01*ansatz01);
	  if ((!TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING)
	      ||(i==j))
	  {
	      val += c3*ansatz00*test00;
	      // stiffness matrix
	      MatrixRowA[j] += Mult * val;
	  }
	  else
	  {
	      // stiffness matrix
	      MatrixRowA[j] += Mult * val;
	      //if (TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING==1)
	      {
		  val += c3*ansatz00*test00;
		  // add to diagonal
		  MatrixRowA[i] += Mult * val;
	      }	     
	  }
      } // endfor j
  } // endfor i
}

void Rhs_Assemble_Bulk(double Mult, double *coeff, double *param,
			    double hK, 
			    double **OrigValues, int *N_BaseFuncts,
			    double ***LocMatrices, double **LocRhs)
{
  double *Rhs;
  double test00; 
  double *Orig2;
  int i, N_;
  double c4; 

  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];
  Orig2 = OrigValues[0];
  c4 = coeff[4];

  for(i=0;i<N_;i++)
  {
      test00 = Orig2[i];		
      Rhs[i] += Mult*test00*c4;
  }
}

// ======================================================================
// assemble matrices and rhs for SUPGtest00
// matrix A contains the part of the equation which is multiplied with
//        theta_1 Delta t
// matrix K contains the part of the SUPG term which comes from the 
//        temporal discretization of the concentration
// ======================================================================
void MatricesA_Assemble_SUPG_Bulk(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs)
{
  double **MatrixA,  **MatrixK, val, *MatrixRowA, *MatrixRowK;
  double ansatz00, ansatz10, ansatz01 ;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_, sold_type=TDatabase::ParamDB->SOLD_TYPE;
  int sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
  double c0, c1, c2, c3, c4, Pe; 
  double tau, bgradv, bb, k, c00, c11, c22, c33, c44, param_sold[5];
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1, theta2;
  double sigma, conc_old, reaction_old, conc_old_x, conc_old_y;

  MatrixA = LocMatrices[0];
  MatrixK= LocMatrices[1];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // u_1
  c2 = coeff[2]; // u_2
  c3 = coeff[3]; // other concentration
  c4 = coeff[4];

  c00 = c0;
  c11 = c1;
  c22 = c2;
  // reactive coefficient, inclusive term from the temporal derivative
  c33 = 1.0/(time_step * theta1) + c3;
  if (fabs(c11) > fabs(c22))
      bb = fabs(c11);
  else
      bb = fabs(c22);
  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, bb);
  for(i=0;i<N_;i++)
  {
      MatrixRowA = MatrixA[i];
      MatrixRowK = MatrixK[i];
      test10 = Orig0[i];
      test01 = Orig1[i];
      test00 = Orig2[i];
      
      // supg test term
      bgradv = tau*(c1*test10+c2*test01);
      
      for(j=0;j<N_;j++)
      {
	  ansatz10 = Orig0[j];
	  ansatz01 = Orig1[j];
	  ansatz00 = Orig2[j];
	  
	  val =  c0*(test10*ansatz10+test01*ansatz01);
	  val += (c1*ansatz10+c2*ansatz01)* (test00+bgradv);
	  if ((!TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING)
	      ||(i==j))
	  {
	      val += c3*ansatz00*(test00+bgradv);
	      // stiffness matrix
	      MatrixRowA[j] += Mult * val;
	  }
	  else
	  {
	      // stiffness matrix
	      MatrixRowA[j] += Mult * val;
	      //if (TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING==1)
	      {
		  val = c3*ansatz00*(test00+bgradv);
		  // add to diagonal
		  MatrixRowA[i] += Mult * val;
	      }	     
	  }
	  // sdfem part which comes from time derivative 
	  MatrixRowK[j] += Mult * ansatz00*bgradv;
      } // endfor j
  } // endfor i
  
  // only SUPG term
  if  (!sold_type)
      return;

  MatrixA = LocMatrices[2];  // SOLD matrix
  theta2 = TDatabase::TimeDB->THETA2;

  // rhs of time-dependent equation for c_A, c_B
  // current velocity is used (as approximation)
  // strong diffusion term is neglected
  // the right hand side is treated in a simplified way
  conc_old = param[6];
  reaction_old = param[7];  
  conc_old_x = param[8];    
  conc_old_y = param[9];
  // rhs of time-dependent equation for c_C
  c44 = conc_old 
      - theta2 * time_step*(c1*conc_old_x + c2 *conc_old_y + reaction_old * conc_old)+ c4;

  //OutPut(conc_old << " " << reaction_old  << " "  << c44 << endl);


  if (sold_parameter_type > BE02_2)
  {
      OutPut("SOLD_PARAMETER_TYPE "<< sold_parameter_type << " not yet implemented!"<<endl);
      exit(4711);
  }
  
  time_step *= theta1;
  c00 = time_step * c0;
  c11 = time_step * c1;
  c22 = time_step * c2;
  c33 = 1.0 + time_step * c3;
  if (fabs(c11) > fabs(c22))
      bb = fabs(c11);
  else
      bb = fabs(c22);
  tau *=  time_step * theta1;

  param_sold[0] = param[10];  // concentration c
  param_sold[1] = param[11];  // c_x
  param_sold[2] = param[12];  // c_y
  param_sold[3] = 0;
  param_sold[4] = 0;
  sigma = Compute_SOLD_sigma(hK, c00, c11, c22, c33, c44, bb, tau, param_sold, 0, 0,1);
  for(i=0;i<N_;i++)
  {
      MatrixRowA = MatrixA[i];
      test10 = Orig0[i];
      test01 = Orig1[i];

      for(j=0;j<N_;j++)
      {
	  ansatz10 = Orig0[j];
	  ansatz01 = Orig1[j];
	  val = 0;
	  if  (sold_type==1)
	  {
	      // isotropic additional diffusion
	      if (bb >0)
		  val = sigma*(test10*ansatz10+test01*ansatz01);
	  }
	  else
	  {
	      // additional diffusion orthogonal to the streamlines
	      if (bb >0)
		  val = sigma * (-c22*ansatz10+c11*ansatz01)*(-c22*test10+c11*test01)/(c11*c11+c22*c22);
	  }
	  MatrixRowA[j] += Mult * val/time_step;
      } // endfor j
  } // endfor i
}
void Rhs_Assemble_SUPG_Bulk(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixK, val, *MatrixRowA, *MatrixRowK;
  double *Rhs;
  double ansatz00, ansatz10, ansatz01 ;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_, sold_type=TDatabase::ParamDB->SOLD_TYPE;
  double c0, c1, c2, c3, c4, Pe; 
  double tau, bgradv, bb, k, c00, c11, c22, c33;
  double theta1 = TDatabase::TimeDB->THETA1;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double sigma, conc_old;
 
  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // u_1
  c2 = coeff[2]; // u_2
  c3 = coeff[3]; // other concentration
  c4 = coeff[4];

  c00 = c0;
  c11 = c1;
  c22 = c2;
  // reactive coefficient, inclusive term from the temporal derivative
  c33 = 1.0/(time_step * theta1) + c3;
  if (fabs(c11) > fabs(c22))
      bb = fabs(c11);
  else
      bb = fabs(c22);
  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, bb);

  for(i=0;i<N_;i++) 
  {
      test10 = Orig0[i];
      test01 = Orig1[i];
      test00 = Orig2[i];
      
      // supg test term
      bgradv = tau*(c1*test10+c2*test01);
      Rhs[i] +=  Mult*(test00 + bgradv)*c4;      
  }
}

// ======================================================================
// assemble matrices and rhs for Galerkin (FCT-FEM)
// matrix A contains the part of the equation which is multiplied with
//        theta_1 Delta t
// ======================================================================
void MatricesA_Assemble_Galerkin_Bulk(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, val, *MatrixRowA;
  double ansatz00, ansatz10, ansatz01 ;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4;

  MatrixA = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // u_1
  c2 = coeff[2]; // u_2
  c3 = coeff[3]; // other concentration
  c4 = coeff[4];

  for(i=0;i<N_;i++)
  {
      MatrixRowA = MatrixA[i];
      test10 = Orig0[i];
      test01 = Orig1[i];
      test00 = Orig2[i];
      
      for(j=0;j<N_;j++)
      {
	  ansatz10 = Orig0[j];
	  ansatz01 = Orig1[j];
	  ansatz00 = Orig2[j];
	  
	  val =  c0*(test10*ansatz10+test01*ansatz01);
	  val += (c1*ansatz10+c2*ansatz01)*  test00;
	  if ((!TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING)
	      ||(i==j))
	  {
	      val += c3*ansatz00*test00;
	      // stiffness matrix
	      MatrixRowA[j] += Mult * val;
	  }
	  else
	  {
	      // stiffness matrix
	      MatrixRowA[j] += Mult * val;
	      //if (TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING==1)
	      {
		  val = c3*ansatz00*test00;
		  // add to diagonal
		  MatrixRowA[i] += Mult * val;
	      }	     
	  }
      } // endfor j
  } // endfor i
}


void MatrixMARhsAssemble_SUPG(double Mult, double *coeff, double *param,
                                 double hK, 
                                 double **OrigValues, int *N_BaseFuncts,
                                 double ***LocMatrices, double **LocRhs)
{
  int i,j, N_T;
  
  double c0, c1, c2, c3, c4;  
  double **MatrixA, **MatrixM, **MatrixK, *Rhs;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;  
  double *MatrixARow, *MatrixMRow, *MatrixKRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, val, val1;
  double c00, c11, c22, c33;
  double tau, bgradv, bb, res, sigma, norm_b;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  
  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixK = LocMatrices[2];
  
  Rhs = LocRhs[0];

  N_T = N_BaseFuncts[0]; 
  
  Orig0 = OrigValues[0]; // T_x
  Orig1 = OrigValues[1]; // T_y
  Orig2 = OrigValues[2]; // T
  Orig3 = OrigValues[3]; // 
  Orig4 = OrigValues[4]; //

  c0 = coeff[0]; // nu
  c1 = coeff[1] + param[0]; // u1-w1 (- to w is added in main prg)
  c2 = coeff[2] + param[1]; // u2-w2 (- to w is added in main prg)
  
  if(TDatabase::ParamDB->P6==1) // con-ALE
   {
    c3 = coeff[3] + param[3]; // concentration coeff - div w (- to w is added in main prg)
   }
  else// non-conservative ALE
   {
    c3 = coeff[3]; 
   }  
  
  c4 = coeff[4]; // rhs

//  cout<< " u1: "<< u1<< " u2: "<< u2<<endl;
//  cout<< "div w: "<< param[3]  <<endl;

  // coefficients for the stabilization parameter
  val = theta1 * time_step;
  c00 = val * c0;
  c11 = val * c1;
  c22 = val * c2;
  // reactive coefficient, inclusive term from the temporal derivative
  c33 = 1.0 + val * c3;
  // reactive coefficient, inclusive term from the temporal derivative
//   c33 = 1.0/(val) + c3;  
  
  if (fabs(c11) > fabs(c22))
      bb = fabs(c11);
  else
      bb = fabs(c22);
  // this is tau
  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, bb);

  for(i=0;i<N_T;i++)
   {
    MatrixARow = MatrixA[i];
    MatrixMRow  = MatrixM[i];
    MatrixKRow  = MatrixK[i];
    
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;
    // scaling with the stabilization parameter
    // scaling with the time step is done in the main program
    bgradv *= tau;
   
    // rhs
    Rhs[i] += Mult*(test00 + bgradv)*c4;

    for(j=0;j<N_T;j++)
     {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];    
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      // Galerkin part of the bilinear form
      val1 = c1*ansatz10+c2*ansatz01;
      val1+= c3*ansatz00;

      val = c0*(test10*ansatz10+test01*ansatz01);
      val += val1*(test00 + bgradv);
      
      // diffusion part of the SUPG stabilization
      val -= c0*(ansatz20 + ansatz02) * bgradv; 
        
      MatrixARow[j] += val*Mult;

      // mass mat
      val = Mult*ansatz00*test00;      
      MatrixMRow[j] += val;
           
      // time consistant term
      val = Mult*ansatz00*bgradv;      
      MatrixKRow[j] += val;        
      
      
     } // for(j=0;j<N_T;j++)     
   } //   for(i=0;i<N_T;i++)
   
}// MatrixMARhsAssemble_SUPG




void MatrixMARhsAssemble(double Mult, double *coeff, double *param,
                                 double hK, 
                                 double **OrigValues, int *N_BaseFuncts,
                                 double ***LocMatrices, double **LocRhs)
{
  int i,j, N_T;
  
  double c0, c1, c2, c3, c4;  
  double **MatrixA, **MatrixM, *Rhs;
  double *Orig0, *Orig1, *Orig2;  
  double *MatrixARow, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, val;
  
  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  
  Rhs = LocRhs[0];

  N_T = N_BaseFuncts[0]; 
  
  Orig0 = OrigValues[0]; // T_x
  Orig1 = OrigValues[1]; // T_y
  Orig2 = OrigValues[2]; // T

  c0 = coeff[0]; // nu
  c1 = coeff[1] + param[0]; // u1-w1 (- to w is added in main prg)
  c2 = coeff[2] + param[1]; // u2-w2 (- to w is added in main prg)  
  c4 = coeff[4]; // rhs  
  
  if(TDatabase::ParamDB->P6==1) // con-ALE
   {
    c3 = coeff[3] + param[3]; // concentration coeff - div w (- to w is added in main prg)
   }
  else// non-conservative ALE
   {
    c3 = coeff[3]; 
//     cout<< "div w: "<< param[0]  <<endl;
   }
//   cout<< " u1: "<< u1<< " u2: "<< u2<<endl;
//  cout<< "div w: "<< param[3]  <<endl;

  for(i=0;i<N_T;i++)
   {
    MatrixARow = MatrixA[i];
    MatrixMRow  = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    // rhs
    Rhs[i] += Mult*test00*c4;

    for(j=0;j<N_T;j++)
     {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];    
    
      //stiffness mat
      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;
   
      MatrixARow[j] += val*Mult;

      // mass mat
      val = Mult*ansatz00*test00;      
      MatrixMRow[j] += val;
      
     } // for(j=0;j<N_T;j++)     
   } //   for(i=0;i<N_T;i++)
   
}// MatrixMARhsAssemble



void MatrixMARhsAssemble_Axial3D(double Mult, double *coeff, double *param,
                                 double hK, 
                                 double **OrigValues, int *N_BaseFuncts,
                                 double ***LocMatrices, double **LocRhs)
{
  int i,j, N_T;
  
  double c0, c1, c2, c3, c4, x, r;  
  double **MatrixA, **MatrixM, *Rhs;
  double *Orig0, *Orig1, *Orig2;  
  double *MatrixARow, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, val;
  
  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  
  Rhs = LocRhs[0];

  N_T = N_BaseFuncts[0]; 
  
  Orig0 = OrigValues[0]; // T_x
  Orig1 = OrigValues[1]; // T_y
  Orig2 = OrigValues[2]; // T

  c0 = coeff[0]; // nu
//   c1 = coeff[1];  
//   c2 = coeff[2];  
  c3 = coeff[3]; // concentration coeff
  c4 = coeff[4]; // rhs
 
  c1 = param[0]; // u1-w1
  c2 = param[1]; // u2-w2
 
  x  = param[2]; // x
  r  = fabs(x);  
   
  //cout<< "x: "<< x<<" u1: "<< u1<< " u2: "<< u2<<endl;
    
  if(r<1e-12)
   {
   OutPut("check NSE2D Axisymmetric file x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as positive points"<<endl);
   }
   

  for(i=0;i<N_T;i++)
   {
    MatrixARow = MatrixA[i];
    MatrixMRow  = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    // rhs
    Rhs[i] += r*Mult*test00*c4;

    for(j=0;j<N_T;j++)
     {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];    
    
      //stiffness mat
      val  = r*c0*(test10*ansatz10+test01*ansatz01);
      val  += r* ((c1*ansatz10+c2*ansatz01)*test00);
      val += r*c3*ansatz00*test00;
  
      MatrixARow[j] += val*Mult;

      // mass mat
      val = r*Mult*ansatz00*test00;      
      MatrixMRow[j] += val;
      
     } // for(j=0;j<N_T;j++)     
   } //   for(i=0;i<N_T;i++)
   
}// MatrixMARhsAssemble_Axial3D


// this routine gives the fe value back
void TimeCDParamsVeloField(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  out[0] = in[2]; // value of first input fe function
  out[1] = in[3]; // value of second input fe function
  out[2] = in[0]; // x value for axial symmetric  
  // cout << out[0] << " " << out[1] <<  " " << out[2] << endl;
}

// this routine gives the fe value back
void TimeCDParamsVeloField_ALE(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  out[0] = in[2]; // value of first input fe function
  out[1] = in[3]; // value of second input fe function
  out[2] = in[0]; // x value for axial symmetric  
  out[3] = in[4] + in[5];  // \nabla\cdot w (divergence of the mesh velo, conservatiove ALE)
//   cout << out[0] << " " << out[1] <<  " " << out[2] << endl;
}

// this routine gives the fe value back
void TimeCDParamsSOLD(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  // cout << in[0] << " " << in[1] << endl;
  out[0] = in[2]; // current solution u
  out[1] = in[3]; // current solution u_x
  out[2] = in[4]; // current solution u_y
  out[3] = in[5]; // current solution u_xx
  out[4] = in[6]; // current solution u_yy
  out[5] = in[7]; // old solution u
  out[6] = in[8]; // old solution u_x
  out[7] = in[9]; // old solution u_y
  out[8] = in[10]; // old solution u_xx
  out[9] = in[11]; // old solution u_yy
  //cout << out[0] << " " << out[1] <<  " " << out[2] << endl; 
}

void TimeCDParamsBulk(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  if (in[2]>0)
     out[0] = in[2]; // value of first input fe function (other concentration)
  else
     out[0] = 0.0;
  out[1] = in[3]; // value of second input fe function (u1)
  out[2] = in[4]; // value of third input fe function (u2)
  // cout << out[0] << " " << out[1] <<  " " << out[2] << endl;
}

void TimeCDParamsBulk_SOLD(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  if (in[2]>0)
     out[0] = in[2]; // value of first input fe function (other concentration)
  else
     out[0] = 0.0;
  out[1] = in[3]; // value of second input fe function (u1)
  out[2] = in[4]; // value of third input fe function (u2)
  // this have to be the same indices as in TimeCDParamsBulk_SOLD_Cc
  out[6] = in[5]; // value of old concentration
  out[7] = in[6]; // value of old reaction (other concentration)
  out[8] = in[7]; // x-deriv of old concentration
  out[9] = in[8]; // y_deriv old concentration
  out[10] = in[9];  // concentration
  out[11] = in[10]; // x-deriv of concentration
  out[12] = in[11]; // y-deriv of concentration
  // cout << out[0] << " " << out[1] <<  " " << out[2] << endl;
}

void TimeCDParamsBulk_Cc(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  out[0] = in[4]; // u1
  out[1] = in[5]; // u2
  out[2] = in[2]; // C_a
  out[3] = in[3]; // C_b
  out[4] = in[6]; // C_c_old
  out[5] = in[7]; // integral_value
}

void TimeCDParamsBulk_SOLD_Cc(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  out[0] = in[4]; // u1
  out[1] = in[5]; // u2
  out[2] = in[2]; // C_a
  out[3] = in[3]; // C_b
  out[4] = in[6]; // C_c_old
  out[5] = in[7]; // integral_value
  out[6] = in[8]; // value of old concentration
  out[7] = 0; // explicit treatment
  out[8] = in[11]; // x-deriv of old concentration
  out[9] = in[12]; // y_deriv of old concentration
  out[10] = in[4]; // x-deriv of concentration
  out[11] = in[13]; // x-deriv of concentration
  out[12] = in[14]; // y_deriv of concentration
}

void TimeCDParamsBulk_mom(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  out[0] = in[2]; // u1
  out[1] = in[3]; // u2
  out[2] = in[4]; // c_C
  out[3] = in[5]; // mom_{k-1}
}

void TimeCDParamsBulk_SOLD_mom(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  out[0] = in[2]; // u1
  out[1] = in[3]; // u2
  out[2] = in[4]; // C_c
  out[3] = in[5]; // C_c_old
}
