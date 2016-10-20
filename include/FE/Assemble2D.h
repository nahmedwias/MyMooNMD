// =======================================================================
// @(#)Assemble2D.h        1.5 04/13/00
// 
// Purpose:     assemble matrix and right-hand side
//
// Author:      Gunar Matthies (10.08.98)
//
// History:     start of implementation 10.08.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __ASSEMBLE2D__
#define __ASSEMBLE2D__

#include <Constants.h>
#include <FEDatabase2D.h>
#include <LocalAssembling2D.h>

#ifdef __3D__
  #include <Aux2D3D.h>
#endif

/** a function from a finite element space */
void Assemble2D(int n_fespaces, const TFESpace2D **fespaces,
                int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                int n_matrices, TMatrix2D **matrices,
                int n_rhs, double **rhs, const TFESpace2D **ferhs,
                BoundCondFunct2D **BoundaryConditions,
                BoundValueFunct2D **BoundaryValues,
                LocalAssembling2D& la
#ifdef __3D__
                , TAux2D3D *Aux2D3D
#endif
, int AssemblePhaseID = -1 
               );





/** a function from a finite element space */
void Assemble2D( int n_fespaces, const TFESpace2D** fespaces, int n_sqmatrices, 
TSquareMatrix2D** sqmatrices, int n_matrices, TMatrix2D** matrices, int n_rhs, 
double** rhs, const TFESpace2D** ferhs, TDiscreteForm2D* DiscreteForm, 
BoundCondFunct2D** BoundaryConditions, BoundValueFunct2D** BoundaryValues, 
TAuxParam2D* Parameters, int AssemblePhaseID = -1 
               );

/** assembling of matrices multiplied by a factor */
void Assemble2D_FCT(int n_fespaces, TFESpace2D **fespaces, int n_sqmatrices, 
                    TSquareMatrix2D **sqmatrices, int n_matrices, 
                    TMatrix2D **matrices, int n_rhs, double **rhs, 
                    TFESpace2D **ferhs, TDiscreteForm2D *DiscreteForm,
                    BoundCondFunct2D **BoundaryConditions,
                    BoundValueFunct2D **BoundaryValues, TAuxParam2D *Parameters,
                    double factor
#ifdef __3D__
                    , TAux2D3D *Aux2D3D
#endif
                   );



/** assembling of slip type bc */
void Assemble2DSlipBC(int n_fespaces, TFESpace2D **fespaces,
                      int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                      int n_matrices, TMatrix2D **matrices,
                      int n_rhs, double **rhs, TFESpace2D **ferhs,
                      TDiscreteForm2D *DiscreteForm,
                      BoundCondFunct2D **BoundaryConditions,
                      BoundValueFunct2D **BoundaryValues,
                      TAuxParam2D *parameters,
                      TFEFunction2D *u1, TFEFunction2D *u2);

/** assembling for methods which need values on neighbour cells */
void Assemble2D_neigh(int n_fespaces, TFESpace2D **fespaces,
                int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                int n_matrices, TMatrix2D **matrices,
                int n_rhs, double **rhs, TFESpace2D **ferhs,
                TDiscreteForm2D *DiscreteForm,
                BoundCondFunct2D **BoundaryConditions,
                BoundValueFunct2D **BoundaryValues,
                TAuxParam2D *Parameters
#ifdef __3D__
                , TAux2D3D *Aux2D3D
#endif
                );

/** assembling for discontinuous Galerkin discretization */
void Assemble2D_DG(CoeffFct2D *Coeff, int n_fespaces, TFESpace2D **fespaces,
		   int n_sqmatrices, TSquareMatrix2D **sqmatrices,
		   int n_matrices, TMatrix2D **matrices,
		   int n_rhs, double **rhs, TFESpace2D **ferhs,
		   BoundCondFunct2D **BoundaryConditions,
		   BoundValueFunct2D **BoundaryValues,
		   TAuxParam2D *Parameters);

/** assembling for continuous interior penalty discretization */
void Assemble2D_CIP(CoeffFct2D *Coeff,int n_fespaces, TFESpace2D **fespaces,
		    int n_sqmatrices, TSquareMatrix2D **sqmatrices,
		    int n_matrices, TMatrix2D **matrices,
		    int n_rhs, double **rhs, TFESpace2D **ferhs,
		    BoundCondFunct2D **BoundaryConditions,
		    BoundValueFunct2D **BoundaryValues,
		    TAuxParam2D *Parameters);

/** assembling for vector finite elements (Raviart-Thomas (RT) and 
 * Brezzi-Douglas-Marini (BDM)) */
void Assemble2D_VectFE(int n_fespaces, const TFESpace2D **fespaces,
           int n_sqmatrices, TSquareMatrix2D **sqmatrices,
           int n_matrices, TMatrix2D **matrices,
           int n_rhs, double **rhs, const TFESpace2D **ferhs,
           LocalAssembling2D& la,
           BoundCondFunct2D **BoundaryConditions,
           BoundValueFunct2D * const * const BoundaryValues
           );

#ifdef __MORTAR__
void Assemble(TMatrix2D *matrix);

  /** add link term (SDFEM) */
  #ifdef __ADD_LINK_SDFEM__
  void AddLinkSDFEM(TFESpace2D *Space2D, TSquareMatrix2D *Matrix,
  #endif
#endif // __MORTAR__

void Assemble2D(int n_fespaces, TFESpace2D **fespaces,
                int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                int n_matrices, TMatrix2D **matrices,
                int n_rhs, double **rhs, TFESpace2D **ferhs,
                TDiscreteForm2D *DiscreteForm,
                BoundCondFunct2D **BoundaryConditions,
                BoundValueFunct2D **BoundaryValues,
                TAuxParam2D *Parameters,
		TAuxParam2D *ParametersBound, 
		TypeBoundSwitchFunct2D *TypeBoundSwitcher,
		int *CounterBoundaryParam
#ifdef __3D__
                , TAux2D3D *Aux2D3D
#endif
                );
#ifdef __2D__
void ComputeInterfaceConditionDirichlet(TSquareMatrix2D  **sqmatrices,
					TMatrix2D **matrices,
					TFEVectFunct2D *u_other,
					double *rhs,
					int *interface_dof_Array,
					int *interface_dof_other_Array,
					int interface_dof, double density, 
					double density_other);


void ComputeInterfaceConditionStress(TSquareMatrix2D  **sqmatrices,
				     TMatrix2D **matrices,
				     TFEVectFunct2D *u_other,
				     double *rhs,
				     int *interface_dof_Array,
				     int *interface_dof_other_Array,
				     int interface_dof, double density, 
				     double density_other, 
				     double *rhs_for_stress_other);
     
#endif // __2D__

void Assemble2D_VectFE(int n_fespaces, const TFESpace2D** fespaces,
           int n_matrices_assemble,
           std::vector<int> row_space, std::vector<int>col_space, 
           int n_rhs_assemble, std::vector<int> row_space_rhs,
           int n_sqmatrices_stored, TSquareMatrix2D** sqmatrices_stored,
           int n_matrices_stored, TMatrix2D** matrices_stored, int n_rhs_stored,
           double** rhs_stored, const TFESpace2D** ferhs_stored,
           LocalAssembling2D& la_assmble, 
           ManipulateMatrices manipulateMatrices, MatrixVector matrixvector,
           ProjectionMatrix projection_matrix, 
           BoundCondFunct2D** BoundaryConditions, 
           BoundValueFunct2D * const * const BoundaryValues);

void Assemble2D_MixedFEM(int n_fespaces, const TFESpace2D** fespaces,
           int n_sqmatrices_assemble, int n_recmatrices_assemble,
           std::vector<int> row_space, std::vector<int>col_space, 
           int n_rhs_assemble, std::vector<int> row_space_rhs,
           int n_sqmatrices_stored, TSquareMatrix2D** sqmatrices_stored,
           int n_matrices_stored, TMatrix2D** matrices_stored, 
           int n_rhs_stored, double** rhs_stored, const TFESpace2D** ferhs_stored,
           LocalAssembling2D& la_assmble, ManipulateMatrices manipulateMatrices,
           MatrixVector matrixvector, ProjectionMatrix projection_matrix,
           BoundCondFunct2D** BoundaryConditions, 
           BoundValueFunct2D * const * const BoundaryValues);

#endif // __ASSEMBLE2D__

