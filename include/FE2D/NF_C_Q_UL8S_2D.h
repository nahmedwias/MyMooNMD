/*
TNodalFunctional2D(NodalFunctional2D id,
         int n_allfunctionals, int n_edgefunctionals, 
         int n_pointsall, int n_pointSdge, 
         double *xi, double *eta, double *t, 
         DoubleFunctVect *evalall, 
         DoubleFunctVect *evaledge); 
*/

static double NF_C_Q_UL8S_2D_Xi[] = {
-1.0 ,-3.0/4.0 ,-1.0/2.0 ,-1.0/4.0 ,0.0 ,1.0/4.0 ,1.0/2.0 ,3.0/4.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,3.0/4.0 ,1.0/2.0 ,1.0/4.0 ,0.0 ,-1.0/4.0 ,-1.0/2.0 ,-3.0/4.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-3.0/4.0 ,-3.0/4.0 ,-3.0/4.0 ,-3.0/4.0 ,-3.0/4.0 ,-3.0/4.0 ,-1.0/2.0 ,-1.0/2.0 ,-1.0/2.0 ,-1.0/2.0 ,-1.0/2.0 ,-1.0/2.0 ,-1.0/4.0 ,-1.0/4.0 ,-1.0/4.0 ,-1.0/4.0 ,-1.0/4.0 ,-1.0/4.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,1.0/4.0 ,1.0/4.0 ,1.0/4.0 ,1.0/4.0 ,1.0/4.0 ,1.0/4.0 ,1.0/2.0 ,1.0/2.0 ,1.0/2.0 ,1.0/2.0 ,1.0/2.0 ,1.0/2.0 ,3.0/4.0 ,-3.0/4.0
};

static double NF_C_Q_UL8S_2D_Eta[] = {
-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-3.0/4.0 ,-1.0/2.0 ,-1.0/4.0 ,0.0 ,1.0/4.0 ,1.0/2.0 ,3.0/4.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,3.0/4.0 ,1.0/2.0 ,1.0/4.0 ,0.0 ,-1.0/4.0 ,-1.0/2.0 ,-3.0/4.0 ,-3.0/4.0 ,-1.0/2.0 ,-1.0/4.0 ,0.0 ,1.0/4.0 ,1.0/2.0 ,-3.0/4.0 ,-1.0/2.0 ,-1.0/4.0 ,0.0 ,1.0/4.0 ,1.0/2.0 ,-3.0/4.0 ,-1.0/2.0 ,-1.0/4.0 ,0.0 ,1.0/4.0 ,1.0/2.0 ,-3.0/4.0 ,-1.0/2.0 ,-1.0/4.0 ,0.0 ,1.0/4.0 ,1.0/2.0 ,-3.0/4.0 ,-1.0/2.0 ,-1.0/4.0 ,0.0 ,1.0/4.0 ,1.0/2.0 ,-3.0/4.0 ,-1.0/2.0 ,-1.0/4.0 ,0.0 ,1.0/4.0 ,1.0/2.0 ,-3.0/4.0 ,3.0/4.0
};

static double NF_C_Q_UL8S_2D_T[] = { 
-1.0 ,-3.0/4.0 ,-1.0/2.0 ,-1.0/4.0 ,0.0 ,1.0/4.0 ,1.0/2.0 ,3.0/4.0 ,1.0 };

void NF_C_Q_UL8S_2D_EvalAll(TCollection *, TBaseCell *,
                            const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 70*SizeOfDouble);
}

void NF_C_Q_UL8S_2D_EvalEdge(TCollection *, TBaseCell *, int,
                             const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 9*SizeOfDouble);
}

TNodalFunctional2D *NF_C_Q_UL8S_2D_Obj = new TNodalFunctional2D
        (NF_C_Q_UL8S_2D, 70, 9, 70, 9, NF_C_Q_UL8S_2D_Xi, NF_C_Q_UL8S_2D_Eta,
         NF_C_Q_UL8S_2D_T, NF_C_Q_UL8S_2D_EvalAll, NF_C_Q_UL8S_2D_EvalEdge);
