// Third order Brezzi-Douglas-Marini vector element, nonconforming, 2D

// coefficient matrix for the degrees of freedom (this is in fact 32 times that 
// matrix)
static double N_Q_BDM3_2D_CM[484] = {
  0.,  0.,  0.,  -7., -4.,  0.,-20.,  0.,  0.,  0.,  0.,   7.,  4.,  0., 20.,  0., 12.,  0.,  0.,  0.,  0.,  0.,
  4.,  0., 20.,   0.,  0.,  0.,  0., -7., -4.,  0.,-20.,   0.,  0.,  0.,  0.,  7.,  0., 12.,  0.,  0.,  0.,  0.,
  0.,  0.,  0.,   0.,-12.,  0.,-20.,  0.,  0.,  0.,  0.,   0.,-12.,  0.,-20.,  0.,  0.,  0., 60.,  0.,  0.,  0.,
  0., 12.,  0.,  84.,  0.,  0.,  0.,  0.,  0., 12.,  0.,  84.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 36.,  0.,  0.,
  0.,  0.,  0.,   0.,  0.,-12.,  0.,-84.,  0.,  0.,  0.,   0.,  0.,-12.,  0.,-84.,  0.,  0.,  0.,  0., 36.,  0.,
-12.,  0.,-20.,   0.,  0.,  0.,  0.,  0.,-12.,  0.,-20.,   0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 60.,
  0.,  0.,  0.,  42., 12.,  0.,  0.,  0.,  0.,  0.,  0., -42.,-12.,  0.,  0.,  0.,-12.,  0.,  0.,  0.,  0.,  0.,
  0.,  0.,-60.,   0.,  0.,  0.,  0.,  0.,  0.,  0., 60.,   0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
  0.,  0.,  0.,   0.,  0., 24.,  0.,-84.,  0.,  0.,  0.,   0.,  0.,-24.,  0., 84.,  0.,  0.,  0.,  0.,  0.,  0.,
  0., 24.,  0., -84.,  0.,  0.,  0.,  0.,  0.,-24.,  0.,  84.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
  0.,  0.,  0.,   0.,  0.,  0., 60.,  0.,  0.,  0.,  0.,   0.,  0.,  0.,-60.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
-12.,  0.,  0.,   0.,  0.,  0.,  0., 42., 12.,  0.,  0.,   0.,  0.,  0.,  0.,-42.,  0.,-12.,  0.,  0.,  0.,  0.,
  0.,  0.,  0.,   0., 20.,  0.,  0.,  0.,  0.,  0.,  0.,   0., 20.,  0.,  0.,  0.,  0.,  0.,-60.,  0.,  0.,  0.,
  0.,  0.,  0.,-140.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,-140.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
  0.,  0.,  0.,   0.,  0., 36.,  0.,  0.,  0.,  0.,  0.,   0.,  0., 36.,  0.,  0.,  0.,  0.,  0.,  0.,-36.,  0.,
  0.,  0., 60.,   0.,  0.,  0.,  0.,  0.,  0.,  0., 60.,   0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
  0.,  0.,  0.,   0.,  0.,  0., 60.,  0.,  0.,  0.,  0.,   0.,  0.,  0., 60.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
  0.,-36.,  0.,   0.,  0.,  0.,  0.,  0.,  0.,-36.,  0.,   0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,-36.,  0.,  0.,
  0.,  0.,  0.,   0.,  0.,  0.,  0.,140.,  0.,  0.,  0.,   0.,  0.,  0.,  0.,140.,  0.,  0.,  0.,  0.,  0.,  0.,
 20.,  0.,  0.,   0.,  0.,  0.,  0.,  0., 20.,  0.,  0.,   0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,-60.,
  0.,  0.,  0., -35.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  35.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
  0.,  0.,  0.,   0.,  0.,  0.,  0., 35.,  0.,  0.,  0.,   0.,  0.,  0.,  0.,-35.,  0.,  0.,  0.,  0.,  0.,  0.
};


static void N_Q_BDM3_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 22; // number of basis functions
  // monomials x-component and y-component
  double mon_x[22]={1,0,  xi,0,  eta,0,  xi*xi,0,  xi*eta,0,  eta*eta,0,  xi*xi*xi,0,  xi*xi*eta,0,  xi*eta*eta,0,  eta*eta*eta,0,  xi*xi*xi*xi,4*xi*eta*eta*eta};
  double mon_y[22]={0,1,  0,xi,  0,eta,  0,xi*xi,  0,xi*eta,  0,eta*eta,  0,xi*xi*xi,  0,xi*xi*eta,  0,xi*eta*eta,  0,eta*eta*eta,  -4*xi*xi*xi*eta,-eta*eta*eta*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM3_2D_CM[i+j*nBF]*mon_x[j] / 32.;
      values[i+nBF] += N_Q_BDM3_2D_CM[i+j*nBF]*mon_y[j] / 32.;
    }
  }
}

// values of the derivatives in xi direction
static void N_Q_BDM3_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 22; // number of basis functions
  // monomials x-component and y-component
  double mon_x[22]={0,0,  1,0,  0,0,  2*xi,0,  eta,0,  0,0,  3*xi*xi,0,  2*xi*eta,0,  eta*eta,0,  0,0,  4*xi*xi*xi,4*eta*eta*eta};
  double mon_y[22]={0,0,  0,1,  0,0,  0,2*xi,  0,eta,  0,0,  0,3*xi*xi,  0,2*xi*eta,  0,eta*eta,  0,0,  -12*xi*xi*eta,0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM3_2D_CM[i+j*nBF]*mon_x[j] / 32.;
      values[i+nBF] += N_Q_BDM3_2D_CM[i+j*nBF]*mon_y[j] / 32.;
    }
  }
}

// values of the derivatives in eta direction
static void N_Q_BDM3_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 22; // number of basis functions
  // monomials x-component and y-component
  double mon_x[22]={0,0,  0,0,  1,0,  0,0,  xi,0,  2*eta,0,  0,0,  xi*xi,0,  xi*2*eta,0,  3*eta*eta,0,  0,12*xi*eta*eta};
  double mon_y[22]={0,0,  0,0,  0,1,  0,0,  0,xi,  0,2*eta,  0,0,  0,xi*xi,  0,xi*2*eta,  0,3*eta*eta,  -4*xi*xi*xi,-4*eta*eta*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM3_2D_CM[i+j*nBF]*mon_x[j] / 32.;
      values[i+nBF] += N_Q_BDM3_2D_CM[i+j*nBF]*mon_y[j] / 32.;
    }
  }
}

// values of derivatives in xi-xi direction
static void N_Q_BDM3_2D_DeriveXiXi(double xi, double eta, double *values)
{
  int nBF = 22; // number of basis functions
  // monomials x-component and y-component
  double mon_x[22]={0,0,  0,0,  0,0,  2,0,  0,0,  0,0,  6*xi,0,  2*eta,0,  0,0,  0,0,  12*xi*xi,0};
  double mon_y[22]={0,0,  0,0,  0,0,  0,2,  0,0,  0,0,  0,6*xi,  0,2*eta,  0,0,  0,0,  -24*xi*eta,0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM3_2D_CM[i+j*nBF]*mon_x[j] / 32.;
      values[i+nBF] += N_Q_BDM3_2D_CM[i+j*nBF]*mon_y[j] / 32.;
    }
  }
}

// values of derivatives in eta-eta direction
static void N_Q_BDM3_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  int nBF = 22; // number of basis functions
  // monomials x-component and y-component
  double mon_x[22]={0,0,  0,0,  0,0,  0,0,  0,0,  2,0,  0,0,  0,0,  xi*2,0,  6*eta,0,  0,24*xi*eta};
  double mon_y[22]={0,0,  0,0,  0,0,  0,0,  0,0,  0,2,  0,0,  0,0,  0,xi*2,  0,6*eta,  0,-12*eta*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM3_2D_CM[i+j*nBF]*mon_x[j] / 32.;
      values[i+nBF] += N_Q_BDM3_2D_CM[i+j*nBF]*mon_y[j] / 32.;
    }
  }
}

// values of derivatives in xi-eta direction
static void N_Q_BDM3_2D_DeriveXiEta(double xi, double eta, double *values)
{
  int nBF = 22; // number of basis functions
  // monomials x-component and y-component
  double mon_x[22]={0,0,  0,0,  0,0,  0,0,  1,0,  0,0,  0,0,  2*xi,0,  2*eta,0,  0,0,  0,12*eta*eta};
  double mon_y[22]={0,0,  0,0,  0,0,  0,0,  0,1,  0,0,  0,0,  0,2*xi,  0,2*eta,  0,0,  -12*xi*xi,0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM3_2D_CM[i+j*nBF]*mon_x[j] / 32.;
      values[i+nBF] += N_Q_BDM3_2D_CM[i+j*nBF]*mon_y[j] / 32.;
    }
  }
}

// the first dof on each edge is the mean flux, the second on each edge is
// an integral where the integrand is multiplied by x, therefore it has to be
// changed with TBaseFunct2D::ChangeBF, similarly the fourth is an integral 
// where the integrand is multiplied by the third Legendre polynomial which 
// also has the property P_3(-x) = -P_3(x)
static int N_Q_BDM3_2D_ChangeJ0[2] = { 1, 3 };
static int N_Q_BDM3_2D_ChangeJ1[2] = { 5, 7 };
static int N_Q_BDM3_2D_ChangeJ2[2] = { 9, 11 };
static int N_Q_BDM3_2D_ChangeJ3[2] = { 13, 15 };

static int *N_Q_BDM3_2D_Change[4] = {N_Q_BDM3_2D_ChangeJ0,N_Q_BDM3_2D_ChangeJ1,
                                     N_Q_BDM3_2D_ChangeJ2,N_Q_BDM3_2D_ChangeJ3};


// ***********************************************************************

TBaseFunct2D *BF_N_Q_BDM3_2D_Obj = new TBaseFunct2D
        (22, BF_N_Q_BDM3_2D, BFUnitSquare,
         N_Q_BDM3_2D_Funct, N_Q_BDM3_2D_DeriveXi,
         N_Q_BDM3_2D_DeriveEta, N_Q_BDM3_2D_DeriveXiXi,
         N_Q_BDM3_2D_DeriveXiEta, N_Q_BDM3_2D_DeriveEtaEta, 4, 4,
         2, N_Q_BDM3_2D_Change, 2);
