
#include <FESpace2D.h>
#include <TriaAffin.h>
#include <QuadAffin.h>
class TSquareMatrix2D;

int GetSignOfThisDOF(int N_DOF, int DOF);


void projection_matrices(int current_cell, const TFESpace2D* ansatzSpace, 
                         const TFESpace2D* testSpace, double ***locMatrix);

void ProjectionMatricesNSE2D(int current_cell, const TFESpace2D* ansatzSpace, 
                         const TFESpace2D* testSpace, double ***locMatrix);


void ProjectionMatricesTNSE2D(int current_cell, const TFESpace2D* ansatzSpace, 
                         const TFESpace2D* testSpace, double ***locMatrix);

void MatVectMult(double ***inputMat, std::pair<int,int>size, double *inputRhs, 
                 double **outputRhs);
void matrices_reconstruction(double ***inputMat, int *nrowInput, int *ncolInput,
                                double ***outputMat, int *nrowOut, int *ncolOut,
                                double **inputrhs,  int *ndimInput,
                                double **outputrhs, int *ndimOutput);

/**************************************************************************** */
void nonlinear_term_reconstruct(double ***inputMat, int *nrowInput, int *ncolInput,
                                double ***outputMat, int *nrowOut, int *ncolOut,
                                double **inputrhs,  int *ndimIn,
                                double **outputrhs, int *ndimOut);
