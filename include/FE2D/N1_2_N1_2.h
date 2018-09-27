// nonconfoming P1 and Q1 like behaviour on considered edge
static char N1_2_N1_2_Name[] = "N1_2_N1_2";
static char N1_2_N1_2_Desc[] = "nonconforming P1 or Q1 element";
static int N1_2_N1_2_N0 = 1;
static int N1_2_N1_2_N1 = 1;
static int N1_2_N1_2_NPairs = 1;
static int N1_2_N1_2_Pairs[][2] = { {0,1} };
static int N1_2_N1_2_NHanging = 0;
static int *N1_2_N1_2_Hanging = nullptr;
static HNDesc *N1_2_N1_2_HangingTypes = nullptr;
static int **N1_2_N1_2_Coupling = nullptr;
static int N1_2_N1_2_NFarHanging = 0;
static int *N1_2_N1_2_FarHanging = nullptr;
static HNDesc *N1_2_N1_2_FarHangingTypes = nullptr;
static int ****N1_2_N1_2_FarCoupling = nullptr;
static int N1_2_N1_2_NNoopposite = 0;
static int *N1_2_N1_2_Nopposite = nullptr;
static int N1_2_N1_2_NNodes = 2;

TFE2DMapper *N1_2_N1_2 = new TFE2DMapper(N1_2_N1_2_Name, N1_2_N1_2_Desc,
                             N1_2_N1_2_N0, N1_2_N1_2_N1,
                             N1_2_N1_2_NPairs, (int *)N1_2_N1_2_Pairs,
                             N1_2_N1_2_NHanging, N1_2_N1_2_Hanging,
                             N1_2_N1_2_HangingTypes, N1_2_N1_2_Coupling,
                             N1_2_N1_2_NFarHanging, N1_2_N1_2_FarHanging,
                             N1_2_N1_2_FarHangingTypes, N1_2_N1_2_FarCoupling,
                             N1_2_N1_2_NNoopposite, N1_2_N1_2_Nopposite,
                             N1_2_N1_2_NNodes);
