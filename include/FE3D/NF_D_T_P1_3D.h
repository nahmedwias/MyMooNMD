// ***********************************************************************
// P1 element, discontinuous, 3D
// 
// Author:     Sashikumaar Ganesan
//
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
static double NF_D_T_P1_3D_Xi[]   = {0, 1, 0, 0, 0.25 };
static double NF_D_T_P1_3D_Eta[]  = {0, 0, 1, 0, 0.25 };
static double NF_D_T_P1_3D_Zeta[] = {0, 0, 0, 1, 0.25 };

/* face 0                               0 */
static double *NF_D_T_P1_3D_F0_Xi = nullptr;
static double *NF_D_T_P1_3D_F0_Eta = nullptr;
static double *NF_D_T_P1_3D_F0_Zeta = nullptr;

/* face 1                               1 */
static double *NF_D_T_P1_3D_F1_Xi = nullptr;
static double *NF_D_T_P1_3D_F1_Eta = nullptr;
static double *NF_D_T_P1_3D_F1_Zeta = nullptr;

/* face 2                               2 */
static double *NF_D_T_P1_3D_F2_Xi = nullptr;
static double *NF_D_T_P1_3D_F2_Eta = nullptr;
static double *NF_D_T_P1_3D_F2_Zeta = nullptr;

/* face 3                               3 */
static double *NF_D_T_P1_3D_F3_Xi = nullptr;
static double *NF_D_T_P1_3D_F3_Eta = nullptr;
static double *NF_D_T_P1_3D_F3_Zeta = nullptr;

static double *NF_D_T_P1_3D_XiArray[4] = {
                        NF_D_T_P1_3D_F0_Xi,
                        NF_D_T_P1_3D_F1_Xi,
                        NF_D_T_P1_3D_F2_Xi,
                        NF_D_T_P1_3D_F3_Xi };

static double *NF_D_T_P1_3D_EtaArray[4] = {
                        NF_D_T_P1_3D_F0_Eta,
                        NF_D_T_P1_3D_F1_Eta,
                        NF_D_T_P1_3D_F2_Eta,
                        NF_D_T_P1_3D_F3_Eta };

static double *NF_D_T_P1_3D_ZetaArray[4] = {
                        NF_D_T_P1_3D_F0_Zeta,
                        NF_D_T_P1_3D_F1_Zeta,
                        NF_D_T_P1_3D_F2_Zeta,
                        NF_D_T_P1_3D_F3_Zeta };

static double *NF_D_T_P1_3D_T = nullptr;
static double *NF_D_T_P1_3D_S = nullptr;

void NF_D_T_P1_3D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues, double *Functionals)
{
  Functionals[0] = ( 16*PointValues[4] + PointValues[0] + PointValues[1]
                     + PointValues[2] + PointValues[3] )/120;
  Functionals[1] =( -30*PointValues[0] +90*PointValues[1]
                    -30*PointValues[2] -30*PointValues[3] )/120;
  Functionals[2] =( -30*PointValues[0] -30*PointValues[1]
                    +90*PointValues[2] -30*PointValues[3] )/120;
  Functionals[3] =( -30*PointValues[0] -30*PointValues[1]
                    -30*PointValues[2] +90*PointValues[3] )/120;
}

void NF_D_T_P1_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int Joint, 
			   double *PointValues, double *Functionals)
{
}

static int NF_D_T_P1_3D_N_AllFunctionals = 4;
static int NF_D_T_P1_3D_N_PointsAll = 5;
static int NF_D_T_P1_3D_N_FaceFunctionals[] = { 0, 0, 0, 0 };
static int NF_D_T_P1_3D_N_PointsFace[] = { 0, 0, 0, 0 };

TNodalFunctional3D *NF_D_T_P1_3D_Obj = new TNodalFunctional3D
        (NF_D_T_P1_3D, NF_D_T_P1_3D_N_AllFunctionals,
         NF_D_T_P1_3D_N_FaceFunctionals, NF_D_T_P1_3D_N_PointsAll,
         NF_D_T_P1_3D_N_PointsFace,
         NF_D_T_P1_3D_Xi, NF_D_T_P1_3D_Eta, NF_D_T_P1_3D_Zeta,
         NF_D_T_P1_3D_XiArray, NF_D_T_P1_3D_EtaArray,
         NF_D_T_P1_3D_ZetaArray,
         NF_D_T_P1_3D_T, NF_D_T_P1_3D_S,
         NF_D_T_P1_3D_EvalAll, NF_D_T_P1_3D_EvalFace);
