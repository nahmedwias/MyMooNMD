/****************************************************************************************
 *                                                                                      *
 *                         Urea_3d4d.h                                                  *
 *                        -------------                                                 *
 *                                                                                      *
 ***************************************************************************************/

#ifndef __UREA_3D4D__
#define __UREA_3D4D__
double growth_rate(double c, double temp);
double b_nuc(double c, double temp);
double InletPSD(double a);
double psd_boundary_urea(double x, double y, double z, double a,
                    double c, double temp);
int PSD_bound_cond_from_velo_inflow_urea(double x, double y, double z);
void Urea_FWE_FDM_Upwind_4D(TCollection *coll,
			    TFEFunction3D *velocity1, TFEFunction3D *velocity2, TFEFunction3D *velocity3,
                            TFEFunction3D *concent_C,
                            TFEFunction3D *Temp,
                            double *f_old,  double *f_new,
                            double *rhs_psd,
                            int N_x, int N_y, int N_z, int N_a,
                            double *x_coord, double *y_coord, double *z_coord, double *a_coord,
                            double x_min, double x_max, double y_min, double y_max,
                            double z_min, double z_max, double a_min, double a_max,
                            double *velo1, double *velo2, double *velo3,
                            double *concent_C_array, double *Temp_array,
			    int *correspond_3dgrid);


void Urea_BWE_FDM_Upwind_4D(TCollection *coll,
                            TFEFunction3D *velocity1, TFEFunction3D *velocity2, TFEFunction3D *velocity3,
                            TFEFunction3D *concent_C,
                            TFEFunction3D *Temp,
                            double *sol,
                            double *rhs_psd,
                            int *correspond_3dgrid,
                            int N_x, int N_y, int N_z, int N_a,
                            double *x_coord, double *y_coord, double *z_coord, double *a_coord,
                            TSquareMatrix3D *mat);


void Urea_Compute_Neum_To_Diri_FEM_FCT(int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord,
double *z_coord, double *a_coord,
int &N_neum_to_diri,
int* &neum_to_diri,
double* &neum_to_diri_x,
double* &neum_to_diri_y,
double* &neum_to_diri_z,
				       double* &neum_to_diri_a);

void Urea_Build_4D_FEM_FCT_Matrix_Q1(TCollection *coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2,
TFEFunction3D *velocity3,
TFEFunction3D *concent_C,
TFEFunction3D *Temp,
double *sol, double *oldsol,double *rhs_psd,
double *lump_mass_PSD, double *matrix_D_Entries_PSD,
int *correspond_3dgrid,
int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord, double *a_coord,
TSquareMatrix3D *mat,
TSquareMatrix3D *matM_cons,
TSquareMatrix3D *matM,
int *index_test_ansatz,
int N_neum_to_diri,
int *neum_to_diri,
double *neum_to_diri_x,
double *neum_to_diri_y,
double *neum_to_diri_z,
				     double *neum_to_diri_a);
double calculate_q_3(int N, double a_layers_coord, double sol_psd);

void Evaluate_f_at_outflow1(int N_x, int N_y, int N_z, int N_a,
                           double *x_coord, double *y_coord,double *z_coord,
                           double *a_coord, double *f);


void Calculate_Volume_Distribution(TCollection *Coll,int N_x, int N_y, int N_z, int N_a, 
			   double *x_coord, double *y_coord, double *z_coord,
				   double *a_layers_coord, double *f);

void Calculate_PSD_on_node(int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord,
			   double *a_layers_coord, double *sol_psd, double x, double y, double z);

void PrepareAgglomerationBreakage(TCollection *Coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2,
TFEFunction3D *velocity3,
TFEFunction3D *temperature,
int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord,
				  double *a_layers_coord, double *f, double *rhs_new);



void call_apply_integral_operator(int nx, int ny, int nz, int na, double* input, double* v, 
                                  double* grad_v, double* temp, double* output, double* grid,
                                  double* params);
  

#endif


