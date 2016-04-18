// ***********************************************************************
// Brezzi-Douglas-Duran-Fortin element of third order on hexahedra, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_H_BDDF3_3D_NDOF = 72;

// number of dofs on the closure of each joints
static int N_H_BDDF3_3D_JointDOF = 10;

// which local dofs are on the joints
static int N_H_BDDF3_3D_J0[10] = {  0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
static int N_H_BDDF3_3D_J1[10] = { 10,11,12,13,14,15,16,17,18,19 };
static int N_H_BDDF3_3D_J2[10] = { 20,21,22,23,24,25,26,27,28,29 };
static int N_H_BDDF3_3D_J3[10] = { 30,31,32,33,34,35,36,37,38,39 };
static int N_H_BDDF3_3D_J4[10] = { 40,41,42,43,44,45,46,47,48,49 };
static int N_H_BDDF3_3D_J5[10] = { 50,51,52,53,54,55,56,57,58,59 };

static int *N_H_BDDF3_3D_J[6] = { N_H_BDDF3_3D_J0, N_H_BDDF3_3D_J1,
                                N_H_BDDF3_3D_J2, N_H_BDDF3_3D_J3,
                                N_H_BDDF3_3D_J4, N_H_BDDF3_3D_J5 };

// number of inner dofs
static int N_H_BDDF3_3D_NInner = 12;

// array containing the numbers for the inner dofs
static int N_H_BDDF3_3D_Inner[12] = { 60,61,62,63,64,65,66,67,68,69,70,71 };

// number of outer dofs
static int N_H_BDDF3_3D_NOuter = 60;

// array containing the numbers for the outer dofs
static int N_H_BDDF3_3D_Outer[60] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59 };

static char N_H_BDDF3_3D_String[] = "N_H_BDDF3_3D";

TFEDesc3D *FE_N_H_BDDF3_3D_Obj=new TFEDesc3D(N_H_BDDF3_3D_String, N_H_BDDF3_3D_NDOF,
                                           N_H_BDDF3_3D_JointDOF, N_H_BDDF3_3D_J,
                                           N_H_BDDF3_3D_NInner, N_H_BDDF3_3D_Inner,
                                           N_H_BDDF3_3D_NOuter, N_H_BDDF3_3D_Outer);
