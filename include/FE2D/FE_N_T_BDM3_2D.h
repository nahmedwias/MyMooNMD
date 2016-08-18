// ***********************************************************************
// BDM3 BDM vector element, nonconforming , 2D
// History:  26.09.2013 implementation (Markus Wolff)
// ***********************************************************************

// number of degrees of freedom
static int N_T_BDM3_2D_NDOF = 20;

// number of dofs on the closure of the joints
static int N_T_BDM3_2D_JointDOF = 4;

// which local dofs are on the joints
static int N_T_BDM3_2D_J0[4] = { 0, 1, 2, 3};
static int N_T_BDM3_2D_J1[4] = { 4, 5, 6, 7};
static int N_T_BDM3_2D_J2[4] = { 8, 9,10,11};

 
static int *N_T_BDM3_2D_J[3] = { N_T_BDM3_2D_J0,
                                N_T_BDM3_2D_J1,
                                N_T_BDM3_2D_J2
                              };

// number of inner dofs
static int N_T_BDM3_2D_NInner = 8;

// array containing the numbers for the inner dofs
static int N_T_BDM3_2D_Inner[8] = {12,13,14,15,16,17,18,19,};

// number of outer dofs (dofs on edges)
static int N_T_BDM3_2D_NOuter = 12;

// array containing the numbers for the outer dofs
static int N_T_BDM3_2D_Outer[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

static char N_T_BDM3_2D_String[] = "N_T_BDM3_2D";

TFEDesc2D *FE_N_T_BDM3_2D_Obj=new TFEDesc2D(N_T_BDM3_2D_String, N_T_BDM3_2D_NDOF,
                                        N_T_BDM3_2D_JointDOF, N_T_BDM3_2D_J,
                                        N_T_BDM3_2D_NInner, N_T_BDM3_2D_Inner,
                                        N_T_BDM3_2D_NOuter, N_T_BDM3_2D_Outer);
