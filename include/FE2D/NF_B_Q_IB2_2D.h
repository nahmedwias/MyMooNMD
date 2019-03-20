// ***********************************************************************
// internal bubble of degree 2 (in the sense of Q2)
// ***********************************************************************
/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

static double NF_B_Q_IB2_2D_Xi[] = { 0 };
static double NF_B_Q_IB2_2D_Eta[] = { 0 };
static double *NF_B_Q_IB2_2D_T = nullptr;

void NF_B_Q_IB2_2D_EvalAll(TCollection *, const TBaseCell *,
                           const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
}

void NF_B_Q_IB2_2D_EvalEdge(TCollection *, const TBaseCell *, int,
                            const double *, double *)
{
}

TNodalFunctional2D *NF_B_Q_IB2_2D_Obj = new TNodalFunctional2D
        (NF_B_Q_IB2_2D, 1, 0, 1, 0, NF_B_Q_IB2_2D_Xi, NF_B_Q_IB2_2D_Eta,
         NF_B_Q_IB2_2D_T, NF_B_Q_IB2_2D_EvalAll, NF_B_Q_IB2_2D_EvalEdge);
