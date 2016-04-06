
static double NF_C_L_P3_1D_Xi[] = { -1., -0.5, 0.5, 1. };
static double NF_C_L_P3_1D_Eta[] = { 0, 0, 0 };
static double NF_C_L_P3_1D_T[] = { -1, 1 };

void NF_C_L_P3_1D_EvalAll(double *PointValues, double *Functionals)
{

  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[0];
}

void NF_C_L_P3_1D_EvalEdge( double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
}

TNodalFunctional1D *NF_C_L_P3_1D_Obj = new TNodalFunctional1D
        (NF_C_L_P3_1D, 4, 2, 4, 2, NF_C_L_P3_1D_Xi, NF_C_L_P3_1D_Eta,
         NF_C_L_P3_1D_T, NF_C_L_P3_1D_EvalAll, NF_C_L_P3_1D_EvalEdge);
