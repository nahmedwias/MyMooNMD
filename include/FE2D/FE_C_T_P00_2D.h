// ***********************************************************************
// P00 element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_T_P00_2D_NDOF = 1;

// number of dofs on the closure of the joints
static int C_T_P00_2D_JointDOF = 0;

// which local dofs are on the joints
static int *C_T_P00_2D_J[3] = { NULL, NULL, NULL };

// number of inner dofs
static int C_T_P00_2D_NInner = 1;

// array containing the numbers for the inner dofs 
static int C_T_P00_2D_Inner[1]={ 0 };

static char C_T_P00_2D_String[] = "C_T_P00_2D";

TFEDesc2D *FE_C_T_P00_2D_Obj=new TFEDesc2D(C_T_P00_2D_String, C_T_P00_2D_NDOF, C_T_P00_2D_JointDOF,
                              C_T_P00_2D_J, C_T_P00_2D_NInner, C_T_P00_2D_Inner);
