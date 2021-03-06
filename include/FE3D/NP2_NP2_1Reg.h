/*
    TFE3DMapper1Reg(char *name, char *description, int nfine, int ncoarse,
              int n_pairs, int **pairs,
              int n_hanging, int *hanging,
              HNDesc *hangingtypes, int **coupling,
              int n_nodes, int **twistpermutation);
*/

static char NP2_NP2_1Reg_Name[] = "NP2_NP2_1Reg";
static char NP2_NP2_1Reg_Desc[] = "conforming NP2 element, 1-regular";
static int NP2_NP2_1Reg_NFine = 3;
static int NP2_NP2_1Reg_NCoarse = 3;
static int NP2_NP2_1Reg_N_Pairs = 0;
static int *NP2_NP2_1Reg_Pairs0 = nullptr;
static int *NP2_NP2_1Reg_Pairs[1] = { (int *)NP2_NP2_1Reg_Pairs0 };

static int NP2_NP2_1Reg_NNodes = 15;

static int NP2_NP2_1Reg_NHanging = 3;
static int NP2_NP2_1Reg_Hanging[3] = { 12, 13, 14 };
static HNDesc NP2_NP2_1Reg_HangingTypes[3] = { HN_N_P2_3D_0, 
                                               HN_N_P2_3D_1,
                                               HN_N_P2_3D_2 };

static int NP2_NP2_1Reg_HN0[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
static int NP2_NP2_1Reg_HN1[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
static int NP2_NP2_1Reg_HN2[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

static int *NP2_NP2_1Reg_Coupling[3] = { NP2_NP2_1Reg_HN0, 
                                         NP2_NP2_1Reg_HN1,
                                         NP2_NP2_1Reg_HN2 };

static int NP2_NP2_1Reg_TwistPerm0[] = { 0, 1, 2 };
static int NP2_NP2_1Reg_TwistPerm1[] = { 1, 2, 0 };
static int NP2_NP2_1Reg_TwistPerm2[] = { 2, 0, 1 };

static int *NP2_NP2_1Reg_TwistPerm[3] = { NP2_NP2_1Reg_TwistPerm0,
                                          NP2_NP2_1Reg_TwistPerm1,
                                          NP2_NP2_1Reg_TwistPerm2 };

static int NP2_NP2_1Reg_NNoOpposite = 12;
static int NP2_NP2_1Reg_NNoOpposite1[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
static int *NP2_NP2_1Reg_NoOpposite[] = { NP2_NP2_1Reg_NNoOpposite1 };
        
TFE3DMapper1Reg *NP2_NP2_1Reg = new TFE3DMapper1Reg(
        NP2_NP2_1Reg_Name, NP2_NP2_1Reg_Desc,
        NP2_NP2_1Reg_NFine, NP2_NP2_1Reg_NCoarse,
        NP2_NP2_1Reg_N_Pairs, NP2_NP2_1Reg_Pairs,
        NP2_NP2_1Reg_NNoOpposite, NP2_NP2_1Reg_NoOpposite,
        NP2_NP2_1Reg_NHanging, NP2_NP2_1Reg_Hanging,
        NP2_NP2_1Reg_HangingTypes, NP2_NP2_1Reg_Coupling,
        NP2_NP2_1Reg_NNodes, NP2_NP2_1Reg_TwistPerm);
