// Second order Brezzi-Douglas-Marini vector element, nonconforming, 2D

// coefficient matrix for the degrees of freedom (this is in fact 8 times that 
// matrix)
static double N_Q_BDM2_2D_CM[196] = {
 0.,  0.,  0., -1.,  0., -5.,  0.,  0.,  0.,  1.,  0.,  5.,  3.,  0.,
 1.,  0.,  5.,  0.,  0.,  0., -1.,  0., -5.,  0.,  0.,  0.,  0.,  3.,
 0.,  0.,  5.,  2.,  0., -5.,  0.,  0.,  5.,  2.,  0., -5.,  0.,  0.,
 0., -6.,  0.,  0.,  0.,  0.,  0., -6.,  0.,  0.,  0.,  0.,  0.,  0.,
 0.,  0.,  0.,  0.,  6.,  0.,  0.,  0.,  0.,  0.,  6.,  0.,  0.,  0.,
 2.,  0., -5.,  0.,  0.,  5.,  2.,  0., -5.,  0.,  0.,  5.,  0.,  0.,
 0.,  0.,  0.,  3.,  0.,  0.,  0.,  0.,  0., -3.,  0.,  0., -3.,  0.,
 0.,  0.,-15.,  0.,  0.,  0.,  0.,  0., 15.,  0.,  0.,  0.,  0.,  0.,
 0.,  0.,  0.,  0.,  6.,  0.,  0.,  0.,  0.,  0., -6.,  0.,  0.,  0.,
 0.,  6.,  0.,  0.,  0.,  0.,  0., -6.,  0.,  0.,  0.,  0.,  0.,  0.,
 0.,  0.,  0.,  0.,  0., 15.,  0.,  0.,  0.,  0.,  0.,-15.,  0.,  0.,
-3.,  0.,  0.,  0.,  0.,  0.,  3.,  0.,  0.,  0.,  0.,  0.,  0., -3.,
 0.,  0., -5.,  0.,  0.,  0.,  0.,  0., -5.,  0.,  0.,  0.,  0.,  0.,
 0.,  0.,  0.,  0.,  0.,  5.,  0.,  0.,  0.,  0.,  0.,  5.,  0.,  0.
};

static void N_Q_BDM2_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={1,0,  xi,0,  eta,0,  xi*xi,0,  xi*eta,0,  eta*eta,0,  xi*xi*xi,3*xi*eta*eta};
  double mon_y[14]={0,1,  0,xi,  0,eta,  0,xi*xi,  0,xi*eta,  0,eta*eta,  -3*xi*xi*eta,-eta*eta*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j] / 8.;
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j] / 8.;
    }
  }
}

// values of the derivatives in xi direction
static void N_Q_BDM2_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0,  1,0,  0,0,  2*xi,0, eta,0,  0,0,   3*xi*xi,3*eta*eta};
  double mon_y[14]={0,0,  0,1,  0,0,  0,2*xi, 0,eta,  0,0,  -6*xi*eta,0      };
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j] / 8.;
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j] / 8.;
    }
  }
}

// values of the derivatives in eta direction
static void N_Q_BDM2_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0,  0,0,  1,0,  0,0,  xi,0,  2*eta,0,  0,6*xi*eta};
  double mon_y[14]={0,0,  0,0,  0,1,  0,0,  0,xi,  0,2*eta, -3*xi*xi,-3*eta*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j] / 8.;
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j] / 8.;
    }
  }
}

// values of derivatives in xi-xi direction
static void N_Q_BDM2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0, 0,0, 0,0, 2,0, 0,0, 0,0, 6*xi,0};
  double mon_y[14]={0,0, 0,0, 0,0, 0,2, 0,0, 0,0, -6*eta,0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j] / 8.;
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j] / 8.;
    }
  }
}

// values of derivatives in eta-eta direction
static void N_Q_BDM2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0, 0,0, 0,0, 0,0, 0,0, 2,0, 0,6*xi};
  double mon_y[14]={0,0, 0,0, 0,0, 0,0, 0,0, 0,2, 0,-6*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j] / 8.;
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j] / 8.;
    }
  }
}

// values of derivatives in xi-eta direction
static void N_Q_BDM2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0, 0,0, 0,0, 0,0, 1,0, 0,0,  0,6*eta};
  double mon_y[14]={0,0, 0,0, 0,0, 0,0, 0,1, 0,0, -6*xi,0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j] / 8.;
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j] / 8.;
    }
  }
}

// the first dof on each edge is the mean flux, the second on each edge is
// an integral where the integrand is multiplied by x, therefore it has to be
// changed with TBaseFunct2D::ChangeBF
static int N_Q_BDM2_2D_ChangeJ0[1] = { 1 };
static int N_Q_BDM2_2D_ChangeJ1[1] = { 4 };
static int N_Q_BDM2_2D_ChangeJ2[1] = { 7 };
static int N_Q_BDM2_2D_ChangeJ3[1] = { 10 };

static int *N_Q_BDM2_2D_Change[4] = {N_Q_BDM2_2D_ChangeJ0,N_Q_BDM2_2D_ChangeJ1,
                                     N_Q_BDM2_2D_ChangeJ2,N_Q_BDM2_2D_ChangeJ3};

// ***********************************************************************

TBaseFunct2D *BF_N_Q_BDM2_2D_Obj = new TBaseFunct2D
        (14, BF_N_Q_BDM2_2D, BFUnitSquare,
         N_Q_BDM2_2D_Funct, N_Q_BDM2_2D_DeriveXi,
         N_Q_BDM2_2D_DeriveEta, N_Q_BDM2_2D_DeriveXiXi,
         N_Q_BDM2_2D_DeriveXiEta, N_Q_BDM2_2D_DeriveEtaEta, 3, 3,
         1, N_Q_BDM2_2D_Change, 2);
