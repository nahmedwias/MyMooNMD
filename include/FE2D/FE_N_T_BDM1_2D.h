// ***********************************************************************
// P1 Raviart-Thomas vector element, nonconforming , 2D
// History:  12.09.2013 implementation (Markus Wolff)
// ***********************************************************************

// number of degrees of freedom
static int N_T_BDM1_2D_NDOF = 6;

// number of dofs on the closure of the joints
static int N_T_BDM1_2D_JointDOF = 2;

// which local dofs are on the joints
static int N_T_BDM1_2D_J0[2] = { 0, 1 };
static int N_T_BDM1_2D_J1[2] = { 2, 3 };
static int N_T_BDM1_2D_J2[2] = { 4, 5 };

 
static int *N_T_BDM1_2D_J[3] = { N_T_BDM1_2D_J0,
                                N_T_BDM1_2D_J1,
                                N_T_BDM1_2D_J2
                              };

// number of inner dofs
static int N_T_BDM1_2D_NInner = 0;

// array containing the numbers for the inner dofs
static int *N_T_BDM1_2D_Inner = nullptr;

// number of outer dofs (dofs on edges)
static int N_T_BDM1_2D_NOuter = 6;

// array containing the numbers for the outer dofs
static int N_T_BDM1_2D_Outer[6] = { 0, 1, 2, 3, 4, 5};

static char N_T_BDM1_2D_String[] = "N_T_BDM1_2D";

TFEDesc2D *FE_N_T_BDM1_2D_Obj=new TFEDesc2D(N_T_BDM1_2D_String, N_T_BDM1_2D_NDOF,
                                        N_T_BDM1_2D_JointDOF, N_T_BDM1_2D_J,
                                        N_T_BDM1_2D_NInner, N_T_BDM1_2D_Inner,
                                        N_T_BDM1_2D_NOuter, N_T_BDM1_2D_Outer);
