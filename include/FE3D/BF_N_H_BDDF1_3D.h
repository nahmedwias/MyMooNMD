// ***********************************************************************
// Brezzi-Douglas-Duran-Fortin element of first order, 3D
// ***********************************************************************

static double N_H_BDDF1_3D_CM[324] = {
//using tschebyscheff points (see NF_N_H_BDDF1_3D.h)
    -0.0441941738,0.0441941738,0,-0.0441941738,-0,0.0441941738,0,0.0625,0.0625,-0.0441941738,0,0.0441941738,0,-0.0625,-0.0625,-0.0441941738,0,0.0441941738,
    -0.0441941738,0,0.0441941738,0,-0.0625,-0.0625,-0.0441941738,0,0.0441941738,0,0.0625,0.0625,-0.0441941738,0.0441941738,0,-0.0441941738,0.0441941738,-0,
    0,-0.0625,-0.0625,-0.0441941738,0.0441941738,0,-0.0441941738,0.0441941738,0,-0.0441941738,0.0441941738,0,-0.0441941738,-0,0.0441941738,0,0.0625,0.0625,
    0,0,0,0,-0,0,0,0.0625,0.0625,0,0,0,0,0.0625,0.0625,0,0,-0,
    0,0,0,0.0883883476,-0,-0.0883883476,0,0,0,-0.0883883476,0,0.0883883476,0,0,0,0,0,-0,
    0.0883883476,-0.0883883476,0,0,-0,0,0,0,0,0,0,0,0,0,0,-0.0883883476,0,0.0883883476,
    0,0,0,0,-0,0,-0.0883883476,0,0.0883883476,0,0,0,0.0883883476,-0.0883883476,0,0,0,-0,
    0,0,0,0,0.0625,0.0625,0,0,0,0,0.0625,0.0625,0,-0,0,0,0,0,
    0.0883883476,0,-0.0883883476,0,-0,0,0,0,0,0,0,0,0,0,0,-0.0883883476,0.0883883476,-0,
    0,0,0,0,-0,0,-0.0883883476,0.0883883476,0,0,0,0,0.0883883476,0,-0.0883883476,0,0,-0,
    0,0,0,0.0883883476,-0.0883883476,0,0,0,0,-0.0883883476,0.0883883476,0,0,0,0,0,0,-0,
    0,0.0625,0.0625,0,-0,0,0,0,0,0,0,0,0,0,0,0,0.0625,0.0625,
    0,0,0,0.0441941738,-0,-0.0441941738,0,0,0,0.0441941738,0,-0.0441941738,0,0,0,0,0,-0,
    0,0,0,0,-0,0,0.0441941738,-0.0441941738,0,0,0,0,0.0441941738,0,-0.0441941738,0,0,-0,
    0.0441941738,0,-0.0441941738,0,-0,0,0,0,0,0,0,0,0,0,0,0.0441941738,-0.0441941738,-0,
    0,0,0,0,-0,0,-0.0441941738,0,0.0441941738,0,0,0,-0.0441941738,0.0441941738,0,0,0,-0,
    -0.0441941738,0.0441941738,0,0,-0,0,0,0,0,0,0,0,0,0,0,-0.0441941738,0,0.0441941738,
    0,0,0,-0.0441941738,0.0441941738,0,0,0,0,-0.0441941738,0.0441941738,0,0,0,0,0,0,0
};

static void N_H_BDDF1_3D_Funct(double xi, double eta, double zeta,
                               double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={1,0,0,xi,0,0,eta,0,0,zeta,0,0,
                  xi*xi,-2*xi*zeta,0,2*xi*eta,-xi*xi,0};
  double mon_y[]={0,1,0,0,xi,0,0,eta,0,0,zeta,0,
                  -2*xi*eta,0,eta*eta,-eta*eta,0,2*eta*zeta};
  double mon_z[]={0,0,1,0,0,xi,0,0,eta,0,0,zeta,
                  0,zeta*zeta,-2*eta*zeta,0,2*xi*zeta,-zeta*zeta};
  
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveXi(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,1,0,0,0,0,0,0,0,0, 2*xi,-2*zeta,0,2*eta,-2*xi,0};
  double mon_y[]={0,0,0,0,1,0,0,0,0,0,0,0, -2*eta,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,1,0,0,0,0,0,0, 0,0,0,0,2*zeta,0};
  
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveEta(double xi, double eta, double zeta,
                                   double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,1,0,0,0,0,0, 0,0,0,2*xi,0,0};
  double mon_y[]={0,0,0,0,0,0,0,1,0,0,0,0, -2*xi,0,2*eta,-2*eta,0,2*zeta};
  double mon_z[]={0,0,0,0,0,0,0,0,1,0,0,0, 0,0,-2*zeta,0,0,0};
  
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveZeta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,1,0,0, 0,-2*xi,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,2*eta};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,1, 0,2*zeta,-2*eta,0,2*xi,-2*zeta};
  
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveXiXi(double, double, double, double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0, 2,0,0,0,-2,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveXiEta(double, double, double, double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,2,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0, -2,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveXiZeta(double, double, double, double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,-2,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,2,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveEtaEta(double, double, double, double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,2,-2,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveEtaZeta(double, double, double, double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,2};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,-2,0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveZetaZeta(double, double, double, double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,2,0,0,0,-2};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

TBaseFunct3D *BF_N_H_BDDF1_3D_Obj =
new TBaseFunct3D(18, BF_N_H_BDDF1_3D, BFUnitHexahedron,
                 N_H_BDDF1_3D_Funct, N_H_BDDF1_3D_DeriveXi,
                 N_H_BDDF1_3D_DeriveEta, N_H_BDDF1_3D_DeriveZeta,
                 N_H_BDDF1_3D_DeriveXiXi, N_H_BDDF1_3D_DeriveXiEta,
                 N_H_BDDF1_3D_DeriveXiZeta, N_H_BDDF1_3D_DeriveEtaEta,
                 N_H_BDDF1_3D_DeriveEtaZeta, N_H_BDDF1_3D_DeriveZetaZeta,
                 2, 1, 0, nullptr, 3);
