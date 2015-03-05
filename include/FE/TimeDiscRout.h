// ======================================================================
// TimeDiscRout.C       07/12/18
//
// Volker John and Joachim Rang
//
// routines which are used for the different temporal discretizations
//
// ======================================================================

/******************************************************************************/
/*                                                                            */
/* GENERAL                                                                    */
/*                                                                            */
/******************************************************************************/
void SetTimeDiscParameters();

int GetN_SubSteps();

/******************************************************************************/
/*                                                                            */
/* Output of a square matrix                                                  */
/*                                                                            */
/******************************************************************************/
void mat_ausgabe(int n_sqmatrices, TSquareMatrix2D **sqmatrices);

/******************************************************************************/
/*                                                                            */
/* Output of rectangular matrices                                             */
/*                                                                            */
/******************************************************************************/
void rmat_ausgabe(int n_sqmatrices, TMatrix2D **matrices);

/******************************************************************************/
/*                                                                            */
/* Output of square matrix on level k                                         */
/*                                                                            */
/******************************************************************************/
void mat_ausgabe1(int k, TSquareMatrix2D **sqmatrices);

/*                                                                            */
/* ROSENBROCK                                                                 */
/*                                                                            */
/******************************************************************************/

void AllocateAuxiliaryVectorsRB(int &rb_order, int &RB_s,
				double* &RB_m, double* &RB_ms,
				double* &RB_alpha, double* &RB_gamma,
				double* &RB_sigma, double* &RB_RHS_YN,
				double* &old_sol_rbU, double* &rb_mein,
				double* &rb_diff, double* &sol_tilde,				
				double* &B1, double* &B2, double* &dq,				
				double* RB_A, double* RB_C, double *RB_S,
				int N_Unknowns);

void AssembleRHS_RB_DIRK(TFESpace2D **fesp, TFEFunction2D **fefct, TFESpace2D **ferhs,
TFESpace2D **USpaces, TFEFunction2D **U1Array, TFEFunction2D **U2Array,
TDiscreteForm2D *DiscreteForm, TDiscreteForm2D *DiscreteFormRHS,
double **RHSs, double *rhs, TAuxParam2D *aux,
BoundCondFunct2D **BoundaryConditions, BoundValueFunct2D **my_BoundValues,
			 int mg_type, int mg_level, int N_U, int N_Unknowns);


/******************************************************************************/
/*                                                                            */
/* DIRK                                                                       */
/*                                                                            */
/******************************************************************************/

void AllocateAuxiliaryVectorsDIRK(int &rb_order, int &RB_s,
				  int &stiff_acc1, int &stiff_acc2,
				  double* &RB_m, double* &RB_ms,
				  double* &RB_alpha, 
				  double* &RB_RHS_YN,
				  double* &old_sol_rbU, double* &rb_mein,
				  double* &old_sol_rbK, double* &old_rhs,				  
				  double* &rb_diff, double* &sol_tilde,				
				  double* &B1, double* &B2,				
				  double* RB_A,
				  int N_Unknowns);
