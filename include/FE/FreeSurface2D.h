#ifdef __2D__

// ====================================================================
// solve grid equation -To be copied to freesurface2d.h
// ====================================================================
void SolveGridEquation(double **Entries, double *sol, double *rhs,
                       int *KCol, int *RowPtr, int N_DOF);

void Solver_3dia(int N_Splines, double *a, double *b, double *c,
                 double *rhs, double *sol);

// ====================================================================
// modify matrices and rhs due to integrals on free surface
// ====================================================================
void FreeSurfInt(TSquareMatrix2D *A11, TSquareMatrix2D *A12,
                 TSquareMatrix2D *A21, TSquareMatrix2D *A22,
                 double *rhs1, double *rhs2,
                 BoundCondFunct2D *BoundaryCondition,
                 double dt, double factor);

// ====================================================================
// determine grid velocity in whole domain
// ====================================================================
void GetGridVelocity(TMultiGrid2D *GridMG, TFEVectFunct2D *GridPos,
                     TFEVectFunct2D *AuxGridPos,
                     double *Nx, double *Ny,
                     TFEVectFunct2D *Velocity, double dt,
                     TFEVectFunct2D *GridVelocity);

// ====================================================================
// determine new grid position
// ====================================================================
void MoveGrid(TMultiGrid2D *GridMG, TFEVectFunct2D *GridPos,
              double *IsoX, double *IsoY,
              double *Nx, double *Ny,
              TFEVectFunct2D *Velocity, double dt,
              TFEVectFunct2D *NewGridPos, 
              double *NewIsoX, double *NewIsoY,
              TFEVectFunct2D *GridVelocity);

void GetGridVelocity(double **Entries, double *Sol, double *Rhs,
                     int *KCol, int *RowPtr,
                     TFEVectFunct2D *GridPos,
                     TFEVectFunct2D *AuxGridPos,
                     TFEVectFunct2D *Velocity, double dt,
                     TFEVectFunct2D *GridVelocity, int *Velo_CellNo);      

void GetGridVelo_outer(double **Entries, double *Sol, double *Rhs,
                       int *KCol, int *RowPtr,
                       TFEVectFunct2D *GridPos,
                       TFEVectFunct2D *AuxGridPos,
                       TFEVectFunct2D *Velocity, double dt,
                       TFEVectFunct2D *GridVelocity, int *Velo_CellNo);
                       
void MoveGrid_2Phase(double **Entries, double *Sol, double *Rhs,
                     int *KCol, int *RowPtr,
                     TFEVectFunct2D *GridPos,
                     TFEVectFunct2D *NewGridPos,
                     TFEVectFunct2D *Velocity, double dt,  
                     int *Velo_CellNo, int isoupdate);

/** Get the inner angles of the cells in whole domain */
void Getcellangle(TFESpace2D *Space, double *MinMaxAngle);
                                                 
double Volume(TFESpace2D *FESpace);

void ReParametrize_pts(int &N_Edges, TBaseCell **cell, int *EdgeNo, double h_min, double **FreePts);

void IntUn(TFEVectFunct2D *u, double *Nx, double *Ny);

#endif
