// ***********************************************************************
// internal bubble of degree 2 (in the sense of Q2)
// ***********************************************************************
/*
    TNodalFunctional3D(NodalFunctional3D id,
                       int n_allfunctionals, int *n_facefunctionals,
                       int n_pointsall, int *n_pointsface,
                       double *xi, double *eta, double *zeta,
                       double **xiarray, double **etaarray,
                       double **zetaarray,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evalface);
*/

/* for all functionals */
static double NF_B_H_IB2_3D_Xi[]   = {  0 };
static double NF_B_H_IB2_3D_Eta[]  = {  0 };
static double NF_B_H_IB2_3D_Zeta[] = {  0 };

/* face 0                               0 */
static double *NF_B_H_IB2_3D_F0_Xi = nullptr;
static double *NF_B_H_IB2_3D_F0_Eta = nullptr;
static double *NF_B_H_IB2_3D_F0_Zeta = nullptr;

/* face 1                               1 */
static double *NF_B_H_IB2_3D_F1_Xi = nullptr;
static double *NF_B_H_IB2_3D_F1_Eta = nullptr;
static double *NF_B_H_IB2_3D_F1_Zeta = nullptr;

/* face 2                               2 */
static double *NF_B_H_IB2_3D_F2_Xi = nullptr;
static double *NF_B_H_IB2_3D_F2_Eta = nullptr;
static double *NF_B_H_IB2_3D_F2_Zeta = nullptr;

/* face 3                               3 */
static double *NF_B_H_IB2_3D_F3_Xi = nullptr;
static double *NF_B_H_IB2_3D_F3_Eta = nullptr;
static double *NF_B_H_IB2_3D_F3_Zeta = nullptr;

/* face 4                               4 */
static double *NF_B_H_IB2_3D_F4_Xi = nullptr;
static double *NF_B_H_IB2_3D_F4_Eta = nullptr;
static double *NF_B_H_IB2_3D_F4_Zeta = nullptr;

/* face 5                               5 */
static double *NF_B_H_IB2_3D_F5_Xi = nullptr;
static double *NF_B_H_IB2_3D_F5_Eta = nullptr;
static double *NF_B_H_IB2_3D_F5_Zeta = nullptr;

static double *NF_B_H_IB2_3D_XiArray[6] = { 
                        NF_B_H_IB2_3D_F0_Xi,
                        NF_B_H_IB2_3D_F1_Xi,
                        NF_B_H_IB2_3D_F2_Xi,
                        NF_B_H_IB2_3D_F3_Xi,
                        NF_B_H_IB2_3D_F4_Xi,
                        NF_B_H_IB2_3D_F5_Xi };

static double *NF_B_H_IB2_3D_EtaArray[6] = { 
                        NF_B_H_IB2_3D_F0_Eta,
                        NF_B_H_IB2_3D_F1_Eta,
                        NF_B_H_IB2_3D_F2_Eta,
                        NF_B_H_IB2_3D_F3_Eta,
                        NF_B_H_IB2_3D_F4_Eta,
                        NF_B_H_IB2_3D_F5_Eta };

static double *NF_B_H_IB2_3D_ZetaArray[6] = { 
                        NF_B_H_IB2_3D_F0_Zeta,
                        NF_B_H_IB2_3D_F1_Zeta,
                        NF_B_H_IB2_3D_F2_Zeta,
                        NF_B_H_IB2_3D_F3_Zeta,
                        NF_B_H_IB2_3D_F4_Zeta,
                        NF_B_H_IB2_3D_F5_Zeta };

static double *NF_B_H_IB2_3D_T = nullptr;
static double *NF_B_H_IB2_3D_S = nullptr;

void NF_B_H_IB2_3D_EvalAll(const TCollection *, const TBaseCell *,
                           const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
}

void NF_B_H_IB2_3D_EvalFace(const TCollection *, const TBaseCell *, int,
                            const double *, double *)
{
}

static int NF_B_H_IB2_3D_N_AllFunctionals = 1;
static int NF_B_H_IB2_3D_N_PointsAll = 1;
static int NF_B_H_IB2_3D_N_FaceFunctionals[] = { 0, 0, 0, 0, 0, 0 };
static int NF_B_H_IB2_3D_N_PointsFace[] = { 0, 0, 0, 0, 0, 0 };

TNodalFunctional3D *NF_B_H_IB2_3D_Obj = new TNodalFunctional3D
        (NF_B_H_IB2_3D, NF_B_H_IB2_3D_N_AllFunctionals,
         NF_B_H_IB2_3D_N_FaceFunctionals, NF_B_H_IB2_3D_N_PointsAll,
         NF_B_H_IB2_3D_N_PointsFace,
         NF_B_H_IB2_3D_Xi, NF_B_H_IB2_3D_Eta, NF_B_H_IB2_3D_Zeta,
         NF_B_H_IB2_3D_XiArray, NF_B_H_IB2_3D_EtaArray,
         NF_B_H_IB2_3D_ZetaArray,
         NF_B_H_IB2_3D_T, NF_B_H_IB2_3D_S,
         NF_B_H_IB2_3D_EvalAll, NF_B_H_IB2_3D_EvalFace);
