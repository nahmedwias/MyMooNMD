// =======================================================================
// %W% %G%
//
// Class:       TMultiGrid3D
// Purpose:     store all data for a multi grid method in 3d
//
// Author:      Gunar Matthies 26.06.2000
//
// History:     26.06.2000 start of implementation
//
// =======================================================================

#include <MultiGrid3D.h>
#include <FEDatabase3D.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <LinAlg.h>

#include <stdlib.h>
#include <string.h>
#ifdef _MPI  
#include <ParFECommunicator3D.h>
#include <FEFunction3D.h>
#endif

/** constructor */
TMultiGrid3D::TMultiGrid3D(int n_problems, int n_parameters, 
                       double *parameters)
{
  N_Levels = 0;
  
  N_Problems = n_problems;

  N_Parameters = n_parameters;

  Parameters = parameters;
}

/** add new level as finest */
void TMultiGrid3D::AddLevel(TMGLevel3D *MGLevel)
{ 
  MultiGridLevels[N_Levels] = MGLevel;

  FESpaces[N_Levels] = MGLevel->GetFESpace();

  N_Levels++;
}
/** add new level as finest */
void TMultiGrid3D::ReplaceLevel(int i,TMGLevel3D *MGLevel)
{ 
  TMGLevel3D *ret;

  ret = MultiGridLevels[i];
  MultiGridLevels[i] = MGLevel;

//  MultiGridLevels[N_Levels] = MGLevel;

  FESpaces[i] = MGLevel->GetFESpace();

  if (i>=N_Levels)
    N_Levels = i+1;
}

/** restrict solution from finest grid to all coarser grids */
void TMultiGrid3D::RestrictToAllGrids()
{
  int lev, j;
  TMGLevel3D *CurrentLevel, *CoarserLevel;
  double *X1, *R1;
  double *X2, *R2, *Aux;

  for(lev=N_Levels-1;lev>0;lev--)
  {
    CurrentLevel = MultiGridLevels[lev];
    X1 = CurrentLevel->GetSolution();

    CoarserLevel = MultiGridLevels[lev-1];
    X2 = CoarserLevel->GetSolution();
    Aux = CoarserLevel->GetAuxVector(0);

    RestrictFunction(FESpaces[lev-1], FESpaces[lev], X2, 
                     X1, Aux);
  } // endfor lev
} // RestrictToAllGrids


#ifdef _MPI  
/***************************************************************************************
/** Copy Values (own+Hallo) to Own FeFunct sol  
************************************************************************************** */
double tC1=0.0,tC2=0.0;
void CopyValuesToOwnFeFunct(TFEFunction3D *C, TFEFunction3D *OwnC, double *value)
{
int i,j,N_owncells,N_cells,N_LocDof,*OwnBeginIndex,*BeginIndex,*OwnGN,*GN,*DOF,*OwnDOF;
double *Solution,*OwnSolution;	
TCollection *coll,*owncoll;
TBaseCell *cell,*owncell;
TFESpace3D *ScalarSpace,*OwnScalarSpace;
double t1,t2;
#ifdef _MPI
t1 = MPI_Wtime();
#endif
ScalarSpace = C->GetFESpace3D();
OwnScalarSpace = OwnC->GetFESpace3D();
coll = ScalarSpace->GetCollection();
owncoll = OwnScalarSpace->GetCollection();

BeginIndex = ScalarSpace->GetBeginIndex();		
GN = ScalarSpace->GetGlobalNumbers();         

OwnBeginIndex = OwnScalarSpace->GetBeginIndex();		
OwnGN = OwnScalarSpace->GetGlobalNumbers();         

N_owncells = owncoll->GetN_OwnCells();
N_cells = coll->GetN_Cells();

#ifdef _HYBRID
#pragma omp parallel default(shared) private(i,cell)
#pragma omp for schedule(guided) nowait 
#endif
  for(i=0;i<N_cells;i++)
   {
   cell = coll->GetCell(i);
    cell->SetLocalCellNo(i);
   }

  Solution = C->GetValues();
  OwnSolution = OwnC->GetValues();
  N_LocDof = BeginIndex[1]-BeginIndex[0];
  //printf("NLOC %d\n",OwnScalarSpace->GetN_DegreesOfFreedom());

#ifdef _HYBRID
  #pragma omp parallel default(shared) private(i,owncell,DOF,OwnDOF,j)
  #pragma omp for schedule(guided) nowait 
#endif
  for(i=0;i<N_owncells;i++)
  {
    owncell = owncoll->GetCell(i);
    DOF = GN + BeginIndex[owncell->GetLocalCellNo()];
    OwnDOF = OwnGN + OwnBeginIndex[i];

    for(j=0;j<N_LocDof;j++)
     OwnSolution[OwnDOF[j]] = value[DOF[j]];		
    }
#ifdef _MPI
    t2= MPI_Wtime();
    tC1 += (t2-t1);
#endif
} // CopyValuesToOwnFeFunct 


/** Copy Own FeFunct sol To Values (own+Hallo) */
void CopyOwnFeFunctToValues(TFEFunction3D *C,TFEFunction3D *OwnC,double *value, double *aux, int flag)
{
int i,j,N_owncells,N_cells,N_LocDof,*OwnBeginIndex,*BeginIndex,*OwnGN,*GN,*DOF,*OwnDOF, N_FineDOFs;
double *Solution,*OwnSolution;	
TCollection *coll,*owncoll;
TBaseCell *cell,*owncell;
TFESpace3D *ScalarSpace,*OwnScalarSpace;
double t1,t2;
#ifdef _MPI
t1 = MPI_Wtime();
#endif
ScalarSpace = C->GetFESpace3D();
OwnScalarSpace = OwnC->GetFESpace3D();
coll = ScalarSpace->GetCollection();
owncoll = OwnScalarSpace->GetCollection();

BeginIndex = ScalarSpace->GetBeginIndex();		
GN = ScalarSpace->GetGlobalNumbers();         

OwnBeginIndex = OwnScalarSpace->GetBeginIndex();		
OwnGN = OwnScalarSpace->GetGlobalNumbers();         

N_owncells = owncoll->GetN_OwnCells();
N_cells = coll->GetN_Cells();
N_FineDOFs = ScalarSpace->GetN_DegreesOfFreedom();
 
#ifdef _HYBRID
#pragma omp parallel default(shared) private(i,cell)
//{
#pragma omp for schedule(guided) nowait 
#endif
for(i=0;i<N_cells;i++)
{
  cell = coll->GetCell(i);
  cell->SetLocalCellNo(i);
}

Solution = C->GetValues();
OwnSolution = OwnC->GetValues();
N_LocDof = BeginIndex[1]-BeginIndex[0];
//printf("NLOC %d\n",OwnScalarSpace->GetN_DegreesOfFreedom());

if(flag!=0)
  memset(aux, 0, SizeOfDouble*N_FineDOFs);

if(flag == 0)
{
 #ifdef _HYBRID
   #pragma omp parallel default(shared) private(i,owncell,DOF,OwnDOF,j)
   #pragma omp for schedule(guided) nowait 
 #endif
 for(i=0;i<N_owncells;i++)
 {
	owncell = owncoll->GetCell(i);
	DOF = GN + BeginIndex[owncell->GetLocalCellNo()];
	OwnDOF = OwnGN + OwnBeginIndex[i];
	
	for(j=0;j<N_LocDof;j++)
	 {
	    value[DOF[j]] = OwnSolution[OwnDOF[j]];
//          printf("OwnDof %d DOF %d GlobalCellNo %d\n", OwnDOF[j], DOF[j], owncell->GetGlobalCellNo());
	 }	
 }
}
else
{
 #ifdef _HYBRID
//    #pragma omp parallel default(shared) private(i,owncell,DOF,OwnDOF,j)
//    #pragma omp for schedule(guided) nowait 
 #endif
 for(i=0;i<N_owncells;i++)
 {
	owncell = owncoll->GetCell(i);
	DOF = GN + BeginIndex[owncell->GetLocalCellNo()];
	OwnDOF = OwnGN + OwnBeginIndex[i];
	
	for(j=0;j<N_LocDof;j++)
	 {
//  	     #pragma omp critical 
	     if( fabs(aux[DOF[j]])<1.e-10)
	     {
	      {
// 	      #pragma omp atomic 
               value[DOF[j]] += TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SCALAR*OwnSolution[OwnDOF[j]];  
// 	       #pragma omp atomic
	       aux[DOF[j]] += 1.;
	      }
	     }
	 }	
 }
}
#ifdef _MPI
t2 = MPI_Wtime();
tC2 += (t2-t1);
#endif
}// 

double tC = 0.0; 
/** defect restriction from level+1 to level */
void CellsPerDOF(TFESpace3D *FineSpace, double *aux)
{
  int i,j;
  TBaseCell *cell;
  TCollection *CoarseColl, *FineColl;
  FE3D CoarseId, FineId;
  TFE3D *CoarseElement, *FineElement;
  BaseFunct3D FineBF;
  int N_FineCells;
  int N_FineDOFs;
  int  *FineBeginIndex;
  int *FineGlobalNumbers;
  int *FineDOF;
  int N_Fine;
  int *DOF;
  double t1,t2;
#ifdef _MPI
  t1 = MPI_Wtime();
#endif  
  FineColl = FineSpace->GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  FineBeginIndex = FineSpace->GetBeginIndex();
  FineGlobalNumbers = FineSpace->GetGlobalNumbers();
  N_FineDOFs = FineSpace->GetN_DegreesOfFreedom();

   memset(aux, 0, SizeOfDouble*N_FineDOFs);
   
 #ifdef _HYBRID
#pragma omp parallel default(shared) private(i,j,cell,DOF,FineId,FineElement,FineBF,N_Fine)
#pragma omp for schedule(guided) nowait 
#endif
  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
  
    DOF = FineGlobalNumbers+FineBeginIndex[i];
    FineId = FineSpace->GetFE3D(i, cell);
    FineElement = TFEDatabase3D::GetFE3D(FineId);
    FineBF = FineElement->GetBaseFunct3D_ID();
    N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();
    for(j=0;j<N_Fine;j++)
    #pragma omp atomic 
      aux[DOF[j]] += 1;
  }
#ifdef _MPI
  t2 = MPI_Wtime();
  tC += (t2-t1);
#endif
}

#endif  



/** one cycle on level i */
void TMultiGrid3D::Cycle(int i, double &res)
{
  double s;
  
  TMGLevel3D *CurrentLevel, *CoarserLevel;
  double *CurrentSol, *CoarserSol, *CoarserRhs;
  double *CurrentRhs, *CurrentDefect, *CurrentAux;
  double *CurrentAux2, *OldSol, *OldDefect;
  double oldres,reduction, alpha;
  int j, N_DOF, maxit, it, slc,gam;
  double initres, normsol, firstres;

  CurrentLevel = MultiGridLevels[i];
  CurrentDefect = CurrentLevel->GetAuxVector(0);
  CurrentAux = CurrentLevel->GetAuxVector(1);
  slc =0;
  if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR)||
      (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR))
    slc = 1;
  if (slc)
  {
    OldSol = CurrentLevel->GetAuxVector(2);
    OldDefect = CurrentLevel->GetAuxVector(3);
  }
  CurrentAux2 = CurrentDefect;
  CurrentSol = CurrentLevel->GetSolution();
  CurrentRhs = CurrentLevel->GetRhs();
  N_DOF = CurrentLevel->GetN_DOF();
  
#ifdef _MPI  
  TParFECommunicator3D *ParComm, *CoarseParComm; 
  int rank;
 
   ParComm = CurrentLevel->GetParComm();  

   MPI_Comm_rank(ParComm->GetComm(), &rank); 
#endif  
 // OutPut("Norm of B rhs in cycle " <<  sqrt(Ddot(N_DOF,CurrentRhs,CurrentRhs)) <<"i is:"<<i <<endl); 
  
  if(i==0)
  {
    // coarse grid
    // cout << "coarse grid" << endl;
    res = 1;
    maxit =  TDatabase::ParamDB->SC_COARSE_MAXIT_SCALAR;
    it = 0;
    CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, res);
    
   if(TDatabase::ParamDB->SC_VERBOSE>=2 
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif  
     )
     {
      OutPut("residual before on coarse "<<res << endl);
     }
     
    reduction = TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SCALAR*res;
    while ((res>reduction)&&(it<maxit))
    {
      switch(TDatabase::ParamDB->SC_COARSE_SMOOTHER_SCALAR)
      {
        case 1: // Jacobi
          CurrentLevel->Jacobi(CurrentSol, CurrentRhs, CurrentAux,
            N_Parameters, Parameters);
#ifdef _MPI  
         ParComm->CommUpdate(CurrentSol,CurrentRhs);
#endif
          break;
        case 2: // SOR
          CurrentLevel->SOR(CurrentSol, CurrentRhs, CurrentAux,
            N_Parameters, Parameters);
#ifdef _MPI  
         ParComm->CommUpdate(CurrentSol,CurrentRhs);
#endif
          break;
        case 3: // SSOR
          CurrentLevel->SSOR(CurrentSol, CurrentRhs, CurrentAux,
            N_Parameters, Parameters);
#ifdef _MPI  
         ParComm->CommUpdate(CurrentSol,CurrentRhs);
#endif
          break;
        case 4: // ILU
          CurrentLevel->ILU(CurrentSol, CurrentRhs, CurrentDefect,
            N_Parameters, Parameters);
          break;
        case 17: // solution with Gaussian elimination
           CurrentLevel->SolveExact(CurrentSol, CurrentRhs);
           break;
#ifdef _MPI
  #ifdef _HYBRID
	case 5: //SOR_Reorder
	for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->SOR_Color(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
	}
        break;
  #else
       case 5: //SOR_Reorder
	for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->SOR_Re(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
	}
        break;
  #endif
#endif
        default:
           OutPut("Coarse smoother not implemented !! Use coarse smoother 3" << endl);
           CurrentLevel->SSOR(CurrentSol, CurrentRhs, CurrentAux,
                              N_Parameters, Parameters);
      } // endswitch SC_COARSE_SMOOTHER_SCALAR

#ifdef _MPI  
         ParComm->CommUpdate(CurrentSol,CurrentRhs);
#endif
   
      CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, res);
      it++;
      if(TDatabase::ParamDB->SC_VERBOSE>=2
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif
        )
         OutPut("res on coarse: " << res << endl);
    }
  }
  else
  {
    slc =0;

    if (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR)
      slc = 1;
    else if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR)
             &&(i==N_Levels-1))
      slc = 1;
    
    CoarserLevel = MultiGridLevels[i-1];
    CoarserSol = CoarserLevel->GetSolution();
    CoarserRhs = CoarserLevel->GetRhs();
  
#ifdef _MPI  
  CoarseParComm = CoarserLevel->GetParComm();      
#endif  
 
    // smoothing
    CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, oldres);  
    firstres = initres = oldres;  
    normsol = sqrt(Ddot(N_DOF, CurrentSol, CurrentSol));
    if(TDatabase::ParamDB->SC_VERBOSE>=2
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif   
      )
      {
       OutPut("level " << i << " ");
       OutPut("res before presmoothing: " << oldres << endl);
      }
   
    if (slc)
    {
      memcpy(OldSol, CurrentSol, N_DOF*SizeOfDouble);
      memcpy(OldDefect, CurrentDefect, N_DOF*SizeOfDouble);
    }

    switch(TDatabase::ParamDB->SC_SMOOTHER_SCALAR)
    {
      case 1: // Jacobi
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->Jacobi(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
#ifdef _MPI  
         ParComm->CommUpdate(CurrentSol,CurrentRhs);
#endif
	}
        break;
      case 2: // SOR
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->SOR(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
#ifdef _MPI  
         ParComm->CommUpdate(CurrentSol,CurrentRhs);
#endif
	}
        break;
      case 3: // SSOR
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->SSOR(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
#ifdef _MPI  
         ParComm->CommUpdate(CurrentSol,CurrentRhs);
#endif
	}
        break;
      case 4: // ILU
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
        {
          CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, 
                oldres);
          CurrentLevel->ILU(CurrentSol, CurrentRhs, CurrentDefect,
                N_Parameters, Parameters);
        }
        break;
#ifdef _MPI
  #ifdef _HYBRID
	case 5: //SOR_Reorder
	for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->SOR_Color(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
	}
        break;
  #else	
	case 5: //SOR_Reorder
	for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->SOR_Re(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
	}
        break;
  #endif
#endif
      default:
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
          CurrentLevel->SOR(CurrentSol, CurrentRhs, CurrentAux,
            N_Parameters, Parameters);
    } // endswitch SC_SMOOTHER

    // calculate defect
    CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, oldres);
if (TDatabase::ParamDB->SC_VERBOSE>=2
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif 
       )
      {
         OutPut("normsol: " << normsol << " oldres: " << oldres << endl);
        OutPut("level " << i << " ");
        OutPut("res after presmoothing: " << oldres << endl);
         OutPut("Smoothing (" << i << "): " << oldres/normsol << endl);
      }
    // restrict defect
#ifdef _MPI  
        ParComm->CommUpdate(CurrentDefect,CurrentRhs);
        CellsPerDOF(CurrentLevel->GetFESpace(),CurrentAux);
       // CopyValuesToOwnFeFunct(CurrentLevel->GetFEFunction(),  CurrentLevel->GetOwnFEFunction(), CurrentAux);
        memcpy(CurrentAux,CurrentLevel->GetOwnFEFunction()->GetValues(),(CurrentLevel->GetOwnFESpace()->GetN_DegreesOfFreedom())*sizeof(double));

        CopyValuesToOwnFeFunct(CurrentLevel->GetFEFunction(),  CurrentLevel->GetOwnFEFunction(), CurrentDefect);

        // first restrict defect to coarse sol, then copy coarse own sol to all (hallo) coarse rhs  
         DefectRestriction(CoarserLevel->GetOwnFESpace(), CurrentLevel->GetOwnFESpace(),
                           CoarserLevel->GetOwnSolution(), CurrentLevel->GetOwnSolution(),
                           CurrentAux);

         memset(CoarserRhs, 0, SizeOfDouble*CoarserLevel->GetN_DOF());
         CopyOwnFeFunctToValues(CoarserLevel->GetFEFunction(),  CoarserLevel->GetOwnFEFunction(), CoarserRhs, CurrentAux, 0); 
         CoarseParComm->CommUpdateReduce(CoarserSol,CoarserRhs);       
#else
    DefectRestriction(FESpaces[i-1], FESpaces[i],CoarserRhs, CurrentDefect, CurrentAux);
#endif

    CoarserLevel->CorrectDefect(CoarserRhs);
    CoarserLevel->Reset(CoarserSol);
    // coarse grid correction
    // coarse grid correction, apply mg recursively*/
    for(j=0;j<mg_recursions[i];j++)
       Cycle(i-1, res);
    if (TDatabase::ParamDB->SC_MG_CYCLE_SCALAR<1) mg_recursions[i] = 1;              // F--cycle 

    // prolongate correction
#ifdef _MPI  
      // copy sol from own+Hallo fespace to own fespace
      CopyValuesToOwnFeFunct(CoarserLevel->GetFEFunction(),  CoarserLevel->GetOwnFEFunction(),  CoarserSol);
      
      // first restrict to coarse sol, then copy coarse sol to all (hallo) coarse rhs 
      Prolongate(CoarserLevel->GetOwnFESpace(), CurrentLevel->GetOwnFESpace(),
                       CoarserLevel->GetOwnSolution(), CurrentLevel->GetOwnSolution(),
                       CurrentLevel->GetAuxVector(1));

      CopyOwnFeFunctToValues(CurrentLevel->GetFEFunction(),  CurrentLevel->GetOwnFEFunction(), CurrentSol, CurrentLevel->GetAuxVector(1), 1);    
      
      ParComm->CommUpdate(CurrentSol,CurrentRhs);
#else 
    Prolongate(FESpaces[i-1], FESpaces[i], 
                   CoarserSol, CurrentAux2, CurrentAux);

    CurrentLevel->CorrectNodes(CurrentAux2);

    CurrentLevel->Update(CurrentSol, CurrentAux2);
#endif  
    
    CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, oldres);
    initres = oldres;
    normsol = sqrt(Ddot(N_DOF, CurrentSol, CurrentSol));
 /*   if (TDatabase::ParamDB->SC_VERBOSE>=2
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif  
       )
      {
        OutPut("level " << i << " ");
        OutPut("res before postsmoothing: " << oldres << endl);
      }*/
    // smoothing
    switch(TDatabase::ParamDB->SC_SMOOTHER_SCALAR)
    {
      case 1: // Jacobi
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->Jacobi(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
#ifdef _MPI  
          // communicate the values (sol & rhs) to the slave DOFs from master DOF
          ParComm->CommUpdate(CurrentSol,CurrentRhs);
#endif 
	}
        break;
      case 2: // SOR
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->SOR(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
#ifdef _MPI  
          // communicate the values (sol & rhs) to the slave DOFs from master DOF
          ParComm->CommUpdate(CurrentSol,CurrentRhs);
#endif 
	}
        break;
      case 3: // SSOR
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->SSOR(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
#ifdef _MPI  
          // communicate the values (sol & rhs) to the slave DOFs from master DOF
          ParComm->CommUpdate(CurrentSol,CurrentRhs);
#endif 
	}
        break;
      case 4: // ILU
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
        {
          CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, 
                oldres);
          CurrentLevel->ILU(CurrentSol, CurrentRhs, CurrentDefect,
                N_Parameters, Parameters);
        }
        break;
#ifdef _MPI
  #ifdef _HYBRID
	case 5: //SOR_Reorder
	for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->SOR_Color(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
	}
        break;
  #else	
	case 5: //SOR_Reorder
	for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->SOR_Re(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
	}
        break;
  #endif
#endif
      default:
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->SOR(CurrentSol, CurrentRhs, CurrentAux,
            N_Parameters, Parameters);
#ifdef _MPI  
          // communicate the values (sol & rhs) to the slave DOFs from master DOF
          ParComm->CommUpdate(CurrentSol,CurrentRhs);
#endif 
	}
    } // endswitch SC_SMOOTHER_SCALAR

    // calculate defect

    if (slc)
    {
      alpha = CurrentLevel->StepLengthControl(CurrentSol, OldSol,
                          OldDefect,                  
                          N_Parameters,Parameters);       
      
      for (j=0;j<N_DOF;j++)
        CurrentSol[j] = OldSol[j] + alpha *( CurrentSol[j]-OldSol[j]);
    }

    CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, res);
  /*  if (TDatabase::ParamDB->SC_VERBOSE>=2
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif  
    )
      {
        OutPut("level " << i << " ");
        OutPut("res after postsmoothing: " << res);
        OutPut(" rate: " << res/firstres << endl);
        // OutPut("Smoothing2 (" << i << "): " << initres/normsol << endl);
      }*/
  }
}

void TMultiGrid3D::SetDirichletNodes(int i)
{
  int HangingNodeBound, N_Dirichlet;
  TMGLevel3D *CurrentLevel;
  double *X, *R;

  if(i>=N_Levels) return;

  CurrentLevel = MultiGridLevels[i];

  X = CurrentLevel->GetSolution();
  R = CurrentLevel->GetRhs();

  HangingNodeBound = CurrentLevel->GetHangingNodeBound();
  N_Dirichlet = CurrentLevel->GetN_Dirichlet();

  memcpy(X+HangingNodeBound, R+HangingNodeBound, SizeOfDouble*N_Dirichlet);
}

/** set recursion for MultiGrid3D */ 
void TMultiGrid3D::SetRecursion(int levels)
{
  int gam = TDatabase::ParamDB->SC_MG_CYCLE_SCALAR,k;

  // coarsest grid 
  mg_recursions[1] = 1;
  if (gam>0)
    for (k=2;k<=levels;k++)        
      mg_recursions[k] = gam;
  else                /* F -- cycle */
    for (k=2;k<=levels;k++)        
      mg_recursions[k] = 2;
}

