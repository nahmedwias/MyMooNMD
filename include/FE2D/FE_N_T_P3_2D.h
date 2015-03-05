// ***********************************************************************
// P3 element, nonconforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int N_T_P3_2D_NDOF = 12;

// number of dofs on the closure of the joints
static int N_T_P3_2D_JointDOF = 3;

// which local dofs are on the joints
static int N_T_P3_2D_J0[3] = { 0, 3, 6 };
static int N_T_P3_2D_J1[3] = { 1, 4, 7 };
static int N_T_P3_2D_J2[3] = { 2, 5, 8 };
 
static int *N_T_P3_2D_J[3] = { N_T_P3_2D_J0, N_T_P3_2D_J1, N_T_P3_2D_J2 };

// number of inner dofs
static int N_T_P3_2D_NInner = 3;

// array containing the numbers for the inner dofs
static int N_T_P3_2D_Inner[3] = { 9, 10, 11 };

// number of outer dofs
static int N_T_P3_2D_NOuter = 9;

// array containing the numbers for the outer dofs
static int N_T_P3_2D_Outer[9] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };

static char N_T_P3_2D_String[] = "N_T_P3_2D";

TFEDesc2D *FE_N_T_P3_2D_Obj=new TFEDesc2D(N_T_P3_2D_String, N_T_P3_2D_NDOF,
                                        N_T_P3_2D_JointDOF, N_T_P3_2D_J,
                                        N_T_P3_2D_NInner, N_T_P3_2D_Inner,
                                        N_T_P3_2D_NOuter, N_T_P3_2D_Outer);
