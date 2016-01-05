// Second order Raviart-Thomas vector element, nonconforming, 2D

// coefficient matrix for the degrees of freedom (this is in fact 2 times that 
// matrix)
static double N_T_RT2_2D_CM[225] = { 
  0.,   0.,   0.,   0.,   0.,   0.,  -2.,  -6., -10.,   0.,   0.,    0.,    0.,    0.,    0.,
 -2.,   6., -10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,    0.,    0.,    0.,    0.,
  2.,  -9.,  25.,  12.,  -6.,   0.,  22.,  51.,  35., 360., 120., -540., -240., -480., -180.,
  0., -12.,  60.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,    0.,    0.,    0.,    0.,
  0.,   0.,   0.,   0.,   0.,   0.,   0.,  12.,  60.,   0.,   0.,    0.,    0.,    0.,    0.,
 22., -51.,  35.,  12.,   6.,   0.,   2.,   9.,  25., 120., 360., -180., -480., -240., -540.,
  0.,  15.,-105., -40.,  45.,  -5., -50., -90., -40.,-900.,-300., 1620.,  780., 1020.,  360.,
  0.,   0., -60.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,    0.,    0.,    0.,    0.,
-20.,  60., -40., -40., -30.,  10.,   0., -90.,-150.,-600.,-600.,  720.,  840., 1320., 1080.,
  0.,  90.,-150., -40.,  30.,  10., -20., -60., -40.,-600.,-600., 1080., 1320.,  840.,  720.,
  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., -60.,   0.,   0.,    0.,    0.,    0.,    0.,
-50.,  90., -40., -40., -45.,  -5.,   0., -15.,-105.,-300.,-900.,  360., 1020.,  780., 1620.,
  0.,   0.,  90.,  30., -45.,  15.,  30.,  45.,  15., 540., 180.,-1080., -540., -540., -180.,
  0., -90.,  90.,  60.,   0., -60.,   0.,  90.,  90., 720., 720.,-1080.,-1440.,-1440.,-1080.,
 30., -45.,  15.,  30.,  45.,  15.,   0.,   0.,  90., 180., 540., -180., -540., -540.,-1080.
};

static void N_T_RT2_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component and y-component
  double mon_x[15]={1,0,  xi,0,  eta,0,  xi*xi,0,  xi*eta,0,  eta*eta,0,  xi*xi*xi,  xi*xi*eta,  xi*eta*eta};
  double mon_y[15]={0,1,  0,xi,  0,eta,  0,xi*xi,  0,xi*eta,  0,eta*eta,  xi*xi*eta, xi*eta*eta, eta*eta*eta };
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT2_2D_CM[i+j*nBF]*mon_x[j] / 2.;
      values[i+nBF] += N_T_RT2_2D_CM[i+j*nBF]*mon_y[j] / 2.;
    }
  }
}

// values of the derivatives in xi direction
static void N_T_RT2_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 15; // number of basis functions
  //Dscal(nBF*nBF,1.0,mat);
  // monomials x-component and y-component
  double mon_x[15]={0,0,  1,0,  0,0,  2*xi,0,  eta,0,  0,0,  3*xi*xi, 2*xi*eta, eta*eta};
  double mon_y[15]={0,0,  0,1,  0,0,  0,2*xi,  0,eta,  0,0,  2*xi*eta, eta*eta, 0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT2_2D_CM[i+j*nBF]*mon_x[j] / 2.;
      values[i+nBF] += N_T_RT2_2D_CM[i+j*nBF]*mon_y[j] / 2.;
    }
  }
}

// values of the derivatives in eta direction
static void N_T_RT2_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component and y-component
  double mon_x[15]={0,0, 0,0, 1,0,  0,0,  xi,0,  2*eta,0,  0,     xi*xi,    2*xi*eta};
  double mon_y[15]={0,0, 0,0, 0,1,  0,0,  0,xi,  0,2*eta,  xi*xi, 2*xi*eta, 3*eta*eta };
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT2_2D_CM[i+j*nBF]*mon_x[j] / 2.;
      values[i+nBF] += N_T_RT2_2D_CM[i+j*nBF]*mon_y[j] / 2.;
    }
  }
}

// values of derivatives in xi-xi direction
static void N_T_RT2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component and y-component
  double mon_x[15]={0,0, 0,0, 0,0, 2,0, 0,0, 0,0, 6*xi,  2*eta, 0};
  double mon_y[15]={0,0, 0,0, 0,0, 0,2, 0,0, 0,0, 2*eta, 0,     0};
 
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT2_2D_CM[i+j*nBF]*mon_x[j] / 2.;
      values[i+nBF] += N_T_RT2_2D_CM[i+j*nBF]*mon_y[j] / 2.;
    }
  }
}

// values of derivatives in eta-eta direction
static void N_T_RT2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component and y-component
  double mon_x[15]={0,0, 0,0, 0,0, 0,0, 0,0, 2,0, 0, 0,    2*xi};
  double mon_y[15]={0,0, 0,0, 0,0, 0,0, 0,0, 0,2, 0, 2*xi, 6*eta };
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT2_2D_CM[i+j*nBF]*mon_x[j] / 2.;
      values[i+nBF] += N_T_RT2_2D_CM[i+j*nBF]*mon_y[j] / 2.;
    }
  }
}

// values of derivatives in xi-eta direction
static void N_T_RT2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component and y-component
  double mon_x[15]={0,0, 0,0, 0,0, 0,0, 1,0, 0,0, 0,    2*xi,  2*eta};
  double mon_y[15]={0,0, 0,0, 0,0, 0,0, 0,1, 0,0, 2*xi, 2*eta, 0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT2_2D_CM[i+j*nBF]*mon_x[j] / 2.;
      values[i+nBF] += N_T_RT2_2D_CM[i+j*nBF]*mon_y[j] / 2.;
    }
  }
}

// the first dof on each edge is the mean flux, the second on each edge is
// an integral where the integrand is multiplied by x, therefore it has to be
// changed with TBaseFunct2D::ChangeBF
static int N_T_RT2_2D_ChangeJ0[1] = { 1 };
static int N_T_RT2_2D_ChangeJ1[1] = { 4 };
static int N_T_RT2_2D_ChangeJ2[1] = { 7 };

static int *N_T_RT2_2D_Change[3] = { N_T_RT2_2D_ChangeJ0, N_T_RT2_2D_ChangeJ1,
                                     N_T_RT2_2D_ChangeJ2};


// ***********************************************************************

TBaseFunct2D *BF_N_T_RT2_2D_Obj = new TBaseFunct2D
        (15, BF_N_T_RT2_2D, BFUnitTriangle, 
         N_T_RT2_2D_Funct, N_T_RT2_2D_DeriveXi,
         N_T_RT2_2D_DeriveEta, N_T_RT2_2D_DeriveXiXi,
         N_T_RT2_2D_DeriveXiEta, N_T_RT2_2D_DeriveEtaEta, 3, 2,
         1, N_T_RT2_2D_Change, 2);
