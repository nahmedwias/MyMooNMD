// ***********************************************************************
// P1 element, discontinuous, 3D
// ***********************************************************************

static void D_H_P1_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  values[0] = 1;
  values[1] = xi;
  values[2] = eta;
  values[3] = zeta;
}

static void D_H_P1_3D_DeriveXi(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 1;
  values[2] = 0;
  values[3] = 0;
}

static void D_H_P1_3D_DeriveEta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 1;
  values[3] = 0;
}

static void D_H_P1_3D_DeriveZeta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 1;
}

static void D_H_P1_3D_DeriveXiXi(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_H_P1_3D_DeriveXiEta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_H_P1_3D_DeriveXiZeta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_H_P1_3D_DeriveEtaEta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_H_P1_3D_DeriveEtaZeta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_H_P1_3D_DeriveZetaZeta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

TBaseFunct3D *BF_D_H_P1_3D_Obj = 
new TBaseFunct3D(4, BF_D_H_P1_3D, BFUnitHexahedron, 
                 D_H_P1_3D_Funct, D_H_P1_3D_DeriveXi,
                 D_H_P1_3D_DeriveEta, D_H_P1_3D_DeriveZeta,
                 D_H_P1_3D_DeriveXiXi, D_H_P1_3D_DeriveXiEta,
                 D_H_P1_3D_DeriveXiZeta, D_H_P1_3D_DeriveEtaEta,
                 D_H_P1_3D_DeriveEtaZeta, D_H_P1_3D_DeriveZetaZeta,
                 1, 1,
                 0, nullptr);
