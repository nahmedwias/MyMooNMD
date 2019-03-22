/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

static double NF_C_Q_UL3S_2D_Xi[] = {
-1.0, -0.3333333333333333, 0.3333333333333333, 1.0, 1.0, 1.0, 1.0,
0.3333333333333333, -0.3333333333333333, -1.0, -1.0, -1.0,
-0.3333333333333333, 0.3333333333333333, -0.3333333333333333 };

static double NF_C_Q_UL3S_2D_Eta[] = {
-1.0, -1.0, -1.0, -1.0, -0.3333333333333333, 0.3333333333333333, 1.0,
1.0, 1.0, 1.0, 0.3333333333333333, -0.3333333333333333,
-0.3333333333333333, -0.3333333333333333, 0.3333333333333333 };

static double NF_C_Q_UL3S_2D_T[] = { -1, -0.33333333333333333333,
                                 0.33333333333333333333, 1 };

void NF_C_Q_UL3S_2D_EvalAll(const TCollection *, const TBaseCell *,
                            const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 15*SizeOfDouble);
}

void NF_C_Q_UL3S_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                             const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 4*SizeOfDouble);
}

TNodalFunctional2D *NF_C_Q_UL3S_2D_Obj = new TNodalFunctional2D
        (NF_C_Q_UL3S_2D, 15, 4, 15, 4, NF_C_Q_UL3S_2D_Xi, NF_C_Q_UL3S_2D_Eta,
         NF_C_Q_UL3S_2D_T, NF_C_Q_UL3S_2D_EvalAll, NF_C_Q_UL3S_2D_EvalEdge);
