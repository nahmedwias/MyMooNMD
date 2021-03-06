/*
    TFE3DMapper1Reg(char *name, char *description, int nfine, int ncoarse,
              int n_pairs, int **pairs,
              int n_hanging, int *hanging,
              HNDesc *hangingtypes, int **coupling,
              int n_nodes, int **twistpermutation);
*/
static char P3_P3_1Reg_Name[] = "P3_P3_1Reg";
static char P3_P3_1Reg_Desc[] = "conforming P3 element, 1-regular";
static int P3_P3_1Reg_NFine = 10;
static int P3_P3_1Reg_NCoarse = 10;
static int P3_P3_1Reg_N_Pairs = 22;
static int P3_P3_1Reg_Pairs0[][2] = { {0,40}, {2,44}, {3,19}, {6,34},
                                      {7,41}, {8,37}, {9,23},
                                      {10,49}, {12,48}, {13,29}, {16,32},
                                      {17,47}, {18,31}, {19,30},
                                      {20,43}, {22,42}, {23,39}, {26,38},
                                      {27,46}, {28,36}, {29,33},
                                      {35,45} };
static int *P3_P3_1Reg_Pairs[1] = { (int *)P3_P3_1Reg_Pairs0 };

static int P3_P3_1Reg_NNodes = 50;

static int P3_P3_1Reg_NHanging = 18;
static int P3_P3_1Reg_Hanging[18] = { 1, 4, 11, 14, 21, 24, 3, 13, 23, 
                                     5, 15, 25, 6, 8, 16, 18, 26, 28 };

static HNDesc P3_P3_1Reg_HangingTypes[18] = {
        HN_C_P3_3D_E, HN_C_P3_3D_E,
        HN_C_P3_3D_E, HN_C_P3_3D_E,
        HN_C_P3_3D_E, HN_C_P3_3D_E,
        HN_C_P3_3D_M, HN_C_P3_3D_M,
        HN_C_P3_3D_M,
        HN_C_P3_3D_F, HN_C_P3_3D_F,
        HN_C_P3_3D_F,
        HN_C_P3_3D_G, HN_C_P3_3D_G,
        HN_C_P3_3D_G, HN_C_P3_3D_G,
        HN_C_P3_3D_G, HN_C_P3_3D_G };

static int P3_P3_1Reg_HN0[] = { 0, 2, 17, 10 };
static int P3_P3_1Reg_HN1[] = { 0, 7, 22, 20 };
static int P3_P3_1Reg_HN2[] = { 10, 12, 27, 20 };
static int P3_P3_1Reg_HN3[] = { 10, 17, 2, 0 };
static int P3_P3_1Reg_HN4[] = { 20, 22, 7, 0 };
static int P3_P3_1Reg_HN5[] = { 20, 27, 12, 10 };
static int P3_P3_1Reg_HN6[] = { 0, 2, 17, 10 };
static int P3_P3_1Reg_HN7[] = { 10, 12, 27, 20 };
static int P3_P3_1Reg_HN8[] = { 0, 7, 22, 20 };
static int P3_P3_1Reg_HN9[] = { 2, 7, 17, 35, 22, 10, 12, 27, 20 };
static int P3_P3_1Reg_HN10[] = { 12, 17, 27, 35, 2, 20, 22, 7, 0 };
static int P3_P3_1Reg_HN11[] = { 22, 27, 7, 35, 12, 0, 2, 17, 10 };
static int P3_P3_1Reg_HN12[] = { 0, 2, 17, 10, 7, 35, 12, 22, 27, 20 };
static int P3_P3_1Reg_HN13[] = { 0, 7, 22, 20, 2, 35, 27, 17, 12, 10 };
static int P3_P3_1Reg_HN14[] = { 10, 12, 27, 20, 17, 35, 22, 2, 7, 0 };
static int P3_P3_1Reg_HN15[] = { 10, 17, 2, 0, 12, 35, 7, 27, 22, 20 };
static int P3_P3_1Reg_HN16[] = { 20, 22, 7, 0, 27, 35, 2, 12, 17, 10 };
static int P3_P3_1Reg_HN17[] = { 20, 27, 12, 10, 22, 35, 17, 7, 2, 0 };

static int *P3_P3_1Reg_Coupling[18] = { P3_P3_1Reg_HN0, P3_P3_1Reg_HN1,
                                       P3_P3_1Reg_HN2, P3_P3_1Reg_HN3,
                                       P3_P3_1Reg_HN4, P3_P3_1Reg_HN5,
                                       P3_P3_1Reg_HN6, P3_P3_1Reg_HN7,
                                       P3_P3_1Reg_HN8, P3_P3_1Reg_HN9,
                                       P3_P3_1Reg_HN10, P3_P3_1Reg_HN11,
                                       P3_P3_1Reg_HN12, P3_P3_1Reg_HN13,
                                       P3_P3_1Reg_HN14, P3_P3_1Reg_HN15,
                                       P3_P3_1Reg_HN16, P3_P3_1Reg_HN17 };

static int P3_P3_1Reg_TwistPerm0[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
static int P3_P3_1Reg_TwistPerm1[] = { 3, 6, 8, 9, 2, 5, 7, 1, 4, 0 };
static int P3_P3_1Reg_TwistPerm2[] = { 9, 7, 4, 0, 8, 5, 1, 6, 2, 3 };

static int *P3_P3_1Reg_TwistPerm[3] = { P3_P3_1Reg_TwistPerm0,
                                        P3_P3_1Reg_TwistPerm1,
                                        P3_P3_1Reg_TwistPerm2 };

static int P3_P3_1Reg_NNoOpposite = 0;
static int **P3_P3_1Reg_NoOpposite = nullptr;
        
TFE3DMapper1Reg *P3_P3_1Reg = new TFE3DMapper1Reg(
        P3_P3_1Reg_Name, P3_P3_1Reg_Desc,
        P3_P3_1Reg_NFine, P3_P3_1Reg_NCoarse,
        P3_P3_1Reg_N_Pairs, P3_P3_1Reg_Pairs,
        P3_P3_1Reg_NNoOpposite, P3_P3_1Reg_NoOpposite,
        P3_P3_1Reg_NHanging, P3_P3_1Reg_Hanging,
        P3_P3_1Reg_HangingTypes, P3_P3_1Reg_Coupling,
        P3_P3_1Reg_NNodes, P3_P3_1Reg_TwistPerm);
