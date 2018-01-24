// ***********************************************************************
// Q1 element, discontinuous, 3D
// ***********************************************************************

// number of degrees of freedom
static int D_H_Q1_3D_NDOF = 8;

// number of dofs on the closure of the joints
static int D_H_Q1_3D_JointDOF = 0;

// which local dofs are on the joints
static int *D_H_Q1_3D_J0 = nullptr;
static int *D_H_Q1_3D_J1 = nullptr;
static int *D_H_Q1_3D_J2 = nullptr;
static int *D_H_Q1_3D_J3 = nullptr;
static int *D_H_Q1_3D_J4 = nullptr;
static int *D_H_Q1_3D_J5 = nullptr;

static int *D_H_Q1_3D_J[6] = { D_H_Q1_3D_J0, D_H_Q1_3D_J1,
                             D_H_Q1_3D_J2, D_H_Q1_3D_J3,
                             D_H_Q1_3D_J4, D_H_Q1_3D_J5};

// number of inner dofs
static int D_H_Q1_3D_NInner = 8;

// array containing the numbers for the inner dofs
static int D_H_Q1_3D_Inner[] = { 0, 1, 2, 3, 4, 5, 6, 7 };

static char D_H_Q1_3D_String[] = "D_H_Q1_3D";

TFEDesc3D *FE_D_H_Q1_3D_Obj=new TFEDesc3D(D_H_Q1_3D_String, D_H_Q1_3D_NDOF, 
                                D_H_Q1_3D_JointDOF,
                                D_H_Q1_3D_J, D_H_Q1_3D_NInner,
                                D_H_Q1_3D_Inner);
