// nonconforming P1 and Q1 like behaviour on considered edge
static char N2_2_N2_2_1Reg_Name[] = "N2_2_N2_2_1Reg";
static char N2_2_N2_2_1Reg_Desc[] = "nonconforming P2 or Q2 element, one regular grid";
static int N2_2_N2_2_1Reg_N0 = 2;
//static int N2_2_N2_2_1Reg_N1 = 2;
static int N2_2_N2_2_1Reg_N2 = 2;
static int N2_2_N2_2_1Reg_NMid = 0;
static int *N2_2_N2_2_1Reg_Mid = nullptr;
static int N2_2_N2_2_1Reg_NPairs = 0;
static int *N2_2_N2_2_1Reg_Pairs = nullptr;
static int N2_2_N2_2_1Reg_NHanging = 2;
static int N2_2_N2_2_1Reg_Hanging[] = { 0, 1 };
static HNDesc N2_2_N2_2_1Reg_HangingTypes[] = { HN_N_P1_2D_0, HN_N_P2_2D_0 };
static int N2_2_N2_2_1Reg_Coupling_0[] = { 2, 4, };
static int N2_2_N2_2_1Reg_Coupling_1[] = { 2, 3, 4, 5 };
static int *N2_2_N2_2_1Reg_Coupling[] = { N2_2_N2_2_1Reg_Coupling_0,
                                          N2_2_N2_2_1Reg_Coupling_1 };
static int N2_2_N2_2_1Reg_NFarHanging = 0;
static int *N2_2_N2_2_1Reg_FarHanging = nullptr;
static HNDesc *N2_2_N2_2_1Reg_FarHangingTypes = nullptr;
static int ****N2_2_N2_2_1Reg_FarCoupling = nullptr;
static int N2_2_N2_2_1Reg_NNoopposite = 4;
static int N2_2_N2_2_1Reg_Nopposite[] = { 2, 3, 4, 5 };
static int N2_2_N2_2_1Reg_NNodes = 6;

TFE2DMapper1Reg *N2_2_N2_2_1Reg = new TFE2DMapper1Reg(
                N2_2_N2_2_1Reg_Name, N2_2_N2_2_1Reg_Desc,
                N2_2_N2_2_1Reg_N0, N2_2_N2_2_1Reg_N2, N2_2_N2_2_1Reg_N2,
                N2_2_N2_2_1Reg_NPairs, (int *)N2_2_N2_2_1Reg_Pairs,
                N2_2_N2_2_1Reg_NMid, (int *)N2_2_N2_2_1Reg_Mid,
                N2_2_N2_2_1Reg_NHanging, N2_2_N2_2_1Reg_Hanging,
                N2_2_N2_2_1Reg_HangingTypes, N2_2_N2_2_1Reg_Coupling,
                N2_2_N2_2_1Reg_NFarHanging, N2_2_N2_2_1Reg_FarHanging,
                N2_2_N2_2_1Reg_FarHangingTypes, N2_2_N2_2_1Reg_FarCoupling,
                N2_2_N2_2_1Reg_NNoopposite, N2_2_N2_2_1Reg_Nopposite,
                N2_2_N2_2_1Reg_NNodes);
