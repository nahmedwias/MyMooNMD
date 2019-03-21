/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

static double NF_C_Q_UL2S_2D_Xi[] = { -1.0, 0.0, 1.0, 1.0, 1.0, 0.0, -1.0, -1.0, 0.0 };
static double NF_C_Q_UL2S_2D_Eta[] = { -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0 };
static double NF_C_Q_UL2S_2D_T[] = { -1, 0, 1 };

void NF_C_Q_UL2S_2D_EvalAll(const TCollection *, const TBaseCell *,
                            const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4];
  Functionals[5] = PointValues[5];
  Functionals[6] = PointValues[6];
  Functionals[7] = PointValues[7];
  Functionals[8] = PointValues[8];
}

void NF_C_Q_UL2S_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                             const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
}

TNodalFunctional2D *NF_C_Q_UL2S_2D_Obj = new TNodalFunctional2D
        (NF_C_Q_UL2S_2D, 9, 3, 9, 3, NF_C_Q_UL2S_2D_Xi, NF_C_Q_UL2S_2D_Eta,
         NF_C_Q_UL2S_2D_T, NF_C_Q_UL2S_2D_EvalAll, NF_C_Q_UL2S_2D_EvalEdge);
