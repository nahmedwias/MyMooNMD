// ***********************************************************************
// P00 element, conforming, 2D
// ***********************************************************************

// base function values
static void C_T_P00_2D_Funct(double, double, double *values)
{
  values[0]=0;
}

// values of the derivatives in xi direction
static void C_T_P00_2D_DeriveXi(double, double, double *values)
{
  values[0]=0;
}

// values of the derivatives in eta direction
static void C_T_P00_2D_DeriveEta(double, double, double *values)
{
  values[0]=0;
}
// values of the derivatives in xi-xi  direction
static void C_T_P00_2D_DeriveXiXi(double, double, double *values)
{
  values[0]=0;
}
// values of the derivatives in xi-eta direction
static void C_T_P00_2D_DeriveXiEta(double, double, double *values)
{
  values[0]=0;
}
// values of the derivatives in eta-eta direction
static void C_T_P00_2D_DeriveEtaEta(double, double, double *values)
{
  values[0]=0;
}
// ***********************************************************************

TBaseFunct2D *BF_C_T_P00_2D_Obj = new TBaseFunct2D
        (1, BF_C_T_P00_2D, BFUnitTriangle, 
         C_T_P00_2D_Funct, C_T_P00_2D_DeriveXi, C_T_P00_2D_DeriveEta,
         C_T_P00_2D_DeriveXiXi, C_T_P00_2D_DeriveXiEta, 
         C_T_P00_2D_DeriveEtaEta, 0, 0,
         0, nullptr);
