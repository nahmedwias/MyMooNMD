// Second order Brezzi-Douglas-Marini vector element, nonconforming, 2D

// coefficient matrix for the degrees of freedom
static double N_T_BDM2_2D_CM[144] = {
 0.,   0. ,   0. ,   0.,   0.,   0.,  -1.,  -3. ,  -5. ,   0.,   0.,    0.,
-1.,   3. ,  -5. ,   0.,   0.,   0.,   0.,   0. ,   0. ,   0.,   0.,    0.,
 0.,  -1.5,  -7.5,  -3.,   3.,   0.,   6.,  10.5,   7.5,  18.,   6.,   90.,
 0.,  -6. ,  30. ,   0.,   0.,   0.,   0.,   0. ,   0. ,   0.,   0.,    0.,
 0.,   0. ,   0. ,   0.,   0.,   0.,   0.,   6. ,  30. ,   0.,   0.,    0.,
 6., -10.5,   7.5,  -3.,  -3.,   0.,   0.,   1.5,  -7.5,   6.,  18.,  -90.,
 1.,   4.5,  12.5,   4.,  -6.,   5.,  -5.,  -7.5,  -2.5, -18.,  -6.,  -90.,
 0.,   0. , -30. ,   0.,   0.,   0.,   0.,   0. ,   0. ,   0.,   0.,    0.,
-2.,  -3. ,   5. ,   4.,   0., -10.,  -2., -15. , -25. , -12., -12., -180.,
-2.,  15. , -25. ,   4.,   0., -10.,  -2.,   3. ,   5. , -12., -12.,  180.,
 0.,   0. ,   0. ,   0.,   0.,   0.,   0.,   0. , -30. ,   0.,   0.,    0.,
-5.,   7.5,  -2.5,   4.,   6.,   5.,   1.,  -4.5,  12.5,  -6., -18.,   90.
};

static void N_T_BDM2_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={1,0,  xi,0,  eta,0,  xi*xi,0,  xi*eta,0,  eta*eta,0};
  double mon_y[12]={0,1,  0,xi,  0,eta,  0,xi*xi,  0,xi*eta,  0,eta*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in xi direction
static void N_T_BDM2_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0,0,  1,0,  0,0,  2*xi,0,  eta,0,  0,0};
  double mon_y[12]={0,0,  0,1,  0,0,  0,2*xi,  0,eta,  0,0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in eta direction
static void N_T_BDM2_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0,0,  0,0,  1,0,  0,0,  xi,0,  2*eta,0};
  double mon_y[12]={0,0,  0,0,  0,1,  0,0,  0,xi,  0,2*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-xi direction
static void N_T_BDM2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0,0,  0,0,  0,0,  2,0,  0,0,  0,0};
  double mon_y[12]={0,0,  0,0,  0,0,  0,2,  0,0,  0,0};
 
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in eta-eta direction
static void N_T_BDM2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0,0,  0,0,  0,0,  0,0,  0,0,  2,0};
  double mon_y[12]={0,0,  0,0,  0,0,  0,0,  0,0,  0,2};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-eta direction
static void N_T_BDM2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0,0,  0,0,  0,0,  0,0,  1,0,  0,0};
  double mon_y[12]={0,0,  0,0,  0,0,  0,0,  0,1,  0,0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// ***********************************************************************

// the first dof on each edge is the mean flux, the second on each edge is
// an integral where the integrand is multiplied by x, therefore it has to be
// changed with TBaseFunct2D::ChangeBF
static int N_T_BDM2_2D_ChangeJ0[1] = { 1 };
static int N_T_BDM2_2D_ChangeJ1[1] = { 4 };
static int N_T_BDM2_2D_ChangeJ2[1] = { 7 };

static int *N_T_BDM2_2D_Change[3] = {N_T_BDM2_2D_ChangeJ0, N_T_BDM2_2D_ChangeJ1,
                                     N_T_BDM2_2D_ChangeJ2 };

TBaseFunct2D *BF_N_T_BDM2_2D_Obj = new TBaseFunct2D
        (12, BF_N_T_BDM2_2D, BFUnitTriangle,
         N_T_BDM2_2D_Funct, N_T_BDM2_2D_DeriveXi,
         N_T_BDM2_2D_DeriveEta, N_T_BDM2_2D_DeriveXiXi,
         N_T_BDM2_2D_DeriveXiEta, N_T_BDM2_2D_DeriveEtaEta, 2, 2,
         1, N_T_BDM2_2D_Change, 2);
