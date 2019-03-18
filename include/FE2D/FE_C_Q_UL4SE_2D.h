// ***********************************************************************
// UL4S element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_UL4SE_2D_NDOF = 20;

// number of dofs on the closure of the joints
static int C_Q_UL4SE_2D_JointDOF = 5;

// which local dofs are on the joints
static int C_Q_UL4SE_2D_J0[5] = {  0,  1,  2,  3,  4 };
static int C_Q_UL4SE_2D_J1[5] = {  4,  5,  6,  7,  8 };
static int C_Q_UL4SE_2D_J2[5] = {  8,  9, 10, 11, 12 };
static int C_Q_UL4SE_2D_J3[5] = { 12, 13, 14, 15,  0 };

static int *C_Q_UL4SE_2D_J[4] = { C_Q_UL4SE_2D_J0, C_Q_UL4SE_2D_J1,
                                  C_Q_UL4SE_2D_J2, C_Q_UL4SE_2D_J3 };

// number of inner dofs
static int C_Q_UL4SE_2D_NInner = 4;

// array containing the numbers for the inner dofs
static int C_Q_UL4SE_2D_Inner[4] = { 16, 17, 18, 19 };

// number of outer dofs
static int C_Q_UL4SE_2D_NOuter = 16;

// array containing the numbers for the outer dofs
static int C_Q_UL4SE_2D_Outer[16] = { 0,  1,  2,  3,  4,  5, 6, 7, 8, 9,
                                     10, 11, 12, 13, 14, 15 };

static char C_Q_UL4SE_2D_String[] = "C_Q_UL4SE_2D";

TFEDesc2D *FE_C_Q_UL4SE_2D_Obj=new TFEDesc2D(C_Q_UL4SE_2D_String, C_Q_UL4SE_2D_NDOF,
                                C_Q_UL4SE_2D_JointDOF, C_Q_UL4SE_2D_J,
                                C_Q_UL4SE_2D_NInner, C_Q_UL4SE_2D_Inner,
                                C_Q_UL4SE_2D_NOuter, C_Q_UL4SE_2D_Outer);
