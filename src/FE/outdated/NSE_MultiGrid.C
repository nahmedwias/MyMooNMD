// =======================================================================
// @(#)NSE_MultiGrid.C        1.24 06/27/00
//
// Class:       TNSE_MultiGrid
// Purpose:     store all data for a multi grid method fo
//              Stokes or Navier-Stokes problems
//
// Author:      Volker John 25.07.2000
//
// History:      25.07.2000 start of implementation
//
// =======================================================================

#include <NSE_MultiGrid.h>
#include <Database.h>
#include <LinAlg.h>
#include <MooNMD_Io.h>
#include <Constants.h>

#ifdef __2D__
  #include <FEDatabase2D.h>
#endif  
#ifdef __3D__
  #include <FEDatabase3D.h>
#endif  

#include <stdlib.h>
#include <string.h>

/** constructor */
TNSE_MultiGrid::TNSE_MultiGrid(int n_problems, int n_parameters, 
                                 double *parameters)
{
  N_Levels = 0;
  
  N_Problems = n_problems;

  N_Parameters = n_parameters;
  
  Parameters = parameters;
}

/** add new level as finest */
void TNSE_MultiGrid::AddLevel(TNSE_MGLevel *MGLevel)
{ 
  MultiGridLevels[N_Levels] = MGLevel;

  USpaces[N_Levels] = MGLevel->GetUSpace();
  PSpaces[N_Levels] = MGLevel->GetPSpace();

  N_Levels++;
}

/** replace level i by given MGLevel and return old level */
TNSE_MGLevel *TNSE_MultiGrid::ReplaceLevel(int i, 
                                           TNSE_MGLevel *MGLevel)
{
  TNSE_MGLevel *ret;

  ret = MultiGridLevels[i];

  MultiGridLevels[i] = MGLevel;

  USpaces[i] = MGLevel->GetUSpace();
  PSpaces[i] = MGLevel->GetPSpace();

  if (i>=N_Levels)
    N_Levels = i+1;
  return ret;
}

/** get parameter */
double TNSE_MultiGrid::GetParam(int i)
{
    if (i<N_Parameters)
	return Parameters[i];
    else
	return -4711;
}

/** get parameter */
void TNSE_MultiGrid::SetParam(int i, double a)
{
    if (i<N_Parameters)
	Parameters[i] = a;
    else
    {
	OutPut("TNSE_MultiGrid::SetParam: not enough parameters " << endl);
	exit(4711);
    }
}


/** restrict u1, u2 from finest grid to all coarser grids */
void TNSE_MultiGrid::RestrictToAllGrids()
{
  TNSE_MGLevel *CurrentLevel, *CoarserLevel; 
  int lev;
  double *CurrentU1, *CoarserU1, *CoarserAux;

  for(lev=N_Levels-1;lev>0;lev--)
  {
    CurrentLevel = MultiGridLevels[lev];
    CurrentLevel->GetSolutionVector(CurrentU1);

    CoarserLevel = MultiGridLevels[lev-1];
    CoarserLevel->GetSolutionVector(CoarserU1);
    CoarserAux = CoarserLevel->GetAuxVector(1);

    RestrictFunction(USpaces[lev-1], USpaces[lev], GEO_DIM,
                     CoarserU1, CurrentU1, CoarserAux);
  } // endfor lev
} // RestrictToAllGrids

/** one cycle on level i */
void TNSE_MultiGrid::Cycle(int i, double &res)
{
  int verbosity = 1;

  int j, N_UDOF, N_PDOF, N_DOF, maxit, smoother; // NLevels=N_Levels;
  int slc, defect_calc, umfpack_flag, ii;
  double oldres; // res2, s;
  
  TNSE_MGLevel *CurrentLevel, *CoarserLevel;
  double *CurrentU;  //*CurrentP;
  double *OldU; // *OldP;
  double *CoarserU, *CoarserP;
  double *CoarserRhsU, *CoarserRhsP;
  double *CurrentRhsU; // *CurrentRhsP;
  double *CurrentOldDefU; // *CurrentOldDefP;
  double *CurrentDefectU, *CurrentDefectP;
  double *CurrentAux, *CurrentCounter;
  double alpha;
  double divfactor = TDatabase::ParamDB->SC_DIV_FACTOR;

  CurrentLevel = MultiGridLevels[i];
  N_UDOF = CurrentLevel->GetN_UDOF();
  N_PDOF = CurrentLevel->GetN_PDOF();
  N_DOF = GEO_DIM*N_UDOF + N_PDOF;
  CurrentLevel->GetSolutionVector(CurrentU);
//  CurrentP = CurrentU + GEO_DIM*N_UDOF;
  CurrentLevel->GetRhsVector(CurrentRhsU);
//  CurrentRhsP = CurrentRhsU + GEO_DIM*N_UDOF;
  CurrentDefectU = CurrentLevel->GetAuxVector(0);
  CurrentDefectP = CurrentDefectU+GEO_DIM*N_UDOF;
  CurrentAux = CurrentLevel->GetAuxVector(1);

  slc =0;
  if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE)||
      (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE))
    slc = 1;
  if (slc)
  {
    CurrentCounter = CurrentLevel->GetAuxVector(2);
    CurrentOldDefU = CurrentCounter;
//    CurrentOldDefP = CurrentOldDefU+GEO_DIM*N_UDOF;
    CurrentCounter = CurrentLevel->GetAuxVector(3);
    OldU = CurrentCounter;
//    OldP = OldU+GEO_DIM*N_UDOF;
  }

  CurrentLevel->CorrectNodes(CurrentU);
 
  // coarsest level ********************************************************* 
  if(i==0)                     
  {
    smoother = TDatabase::ParamDB->SC_COARSE_SMOOTHER_SADDLE;
    switch(smoother)
    {
      case 17 : 
        CurrentLevel->SolveExact(CurrentU, CurrentRhsU);
        CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, res);
        if (verbosity==2)
        {
          if (res>1e-12)
          {
            OutPut("MESSAGE: residual not zero !! ("<<res);
            OutPut(") Check boundary conditions !! "<< endl);
          }
        }
        break;
      case 18 : // Direct Solver with UMFPACK
        //OutPut("This should be done by UMFPACK in the future" << endl);
        umfpack_flag = (int) Parameters[9];
        CurrentLevel->SolveExactUMFPACK(CurrentU, CurrentRhsU, umfpack_flag);
        Parameters[9] = umfpack_flag;
        CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, res);
        if (verbosity>=2)
        {
          if (res>1e-12)
          {
            OutPut("MESSAGE: residual not zero !! ("<<res);
            OutPut(") Check boundary conditions !! "<< endl);
          }
        }
        break;
      case 1 :
      case 2 :
        slc = 0;
        if (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE)
          slc = 1;
        else 
          if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE)
                 &&(i==N_Levels-1))
            slc = 1;

        maxit = TDatabase::ParamDB->SC_COARSE_MAXIT_SADDLE; // max # of iterations
        j = 0;

        // compute defect 
        CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, oldres);
        if (verbosity>=2)
        {
          OutPut("smoother " <<smoother <<" , level 0 res before smoothing: " << oldres << endl);
        }
        // res2 = res = oldres;  // not used later in this function
        divfactor *= res;
        if (slc)
        {
          ii = N_DOF*SizeOfDouble;
          memcpy(OldU, CurrentU, ii);
          memcpy(CurrentOldDefU, CurrentDefectU, ii);
        }
        // iterate, do at least one iteration
        while (((res>TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SADDLE*oldres)
             || (j==0)))
        { 
          // apply smother
          CurrentLevel->CellVanka(CurrentU, CurrentRhsU, CurrentAux,
                                  N_Parameters, Parameters, smoother, N_Levels);
          j++;
          // maxit reached             
          if(j>=maxit) break;
          // compute defect
          CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, res);
          if (res > divfactor)
          {
            OutPut("mesh cell Vanka : coarse grid solver diverged " <<  res  << endl);
            exit(4711);
          }
          // cout << "residual " << j << " " << res << endl;
        }
        if (verbosity >=2)
        {
          // compute defect
          CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, res);
          OutPut("smoother " << smoother <<", level 0 res after smoothing: "<< res << " reduction " <<  res/oldres );
          OutPut(" number of Vanka iters: " << j << endl);
        }
        if (slc)
        { 
          alpha = CurrentLevel->
          StepLengthControl(CurrentU, OldU, CurrentOldDefU,
                            N_Parameters,Parameters);
          for (j=0;j<N_DOF;j++)
            CurrentU[j] = OldU[j] + alpha *( CurrentU[j]-OldU[j]);
        }
        break;
	
      case 3:
      case 4:
        slc =0;
        if (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE)
          slc = 1;
        else
          if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE)
            &&(i==N_Levels-1))
          slc = 1;
            
        maxit = TDatabase::ParamDB->SC_COARSE_MAXIT_SADDLE; // max # of iterations
        j = 0;
        // compute defect 
        CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, oldres);
        // res2 = res = oldres;   // not used later in this function
        if (verbosity>=2)
        {
          OutPut("smoother " << smoother <<"level 0 res before smoothing: " << oldres << endl);
        }
        if (slc)
        {
          memcpy(OldU, CurrentU, N_DOF*SizeOfDouble);
          memcpy(CurrentOldDefU, CurrentDefectU, N_DOF*SizeOfDouble);
        }
        // iterate
        while (((res>TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SADDLE*oldres)
             || (j==0)))
        {
          // apply smother
          CurrentLevel->NodalVanka(CurrentU, CurrentRhsU, CurrentAux, 
                                   N_Parameters, Parameters,smoother,N_Levels);
          j++;
          // maxit reached             
          if(j>=maxit) break;
          // compute defect
          CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, res);
          if (verbosity >=2)
          {
            OutPut("level 0: Vanka ite: " << j << " " << res << endl);
          }

          if (res > divfactor)
          {
            OutPut("nodal Vanka : coarse grid solver diverged " <<  res  << endl);
            exit(4711);
          }
          // res2 = res;  // not used later in this function
          // cout << "residual " << j << " " << res << endl;
        }
        if (verbosity >=2)
        {
          // compute defect
          CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, res);
          OutPut("smoother " << smoother <<"level 0 res after smoothing: "<< res 
          << " reduction " <<  res/oldres );
          OutPut(" number of Vanka iters: " << j << endl);
        }
        if (slc)
        { 
          alpha = CurrentLevel->
            StepLengthControl(CurrentU, OldU, CurrentOldDefU, 
                              N_Parameters,Parameters);       
        
          for (j=0;j<N_DOF;j++)
            CurrentU[j] = OldU[j] + alpha *( CurrentU[j]-OldU[j]);
        }
        break;

      case 11:
        slc =0;
        if (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE)
          slc = 1;
        else
          if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE)
                 &&(i==N_Levels-1))
            slc = 1;
  
        maxit = TDatabase::ParamDB->SC_COARSE_MAXIT_SADDLE; // max # of iterations
        j = 0;
        // compute defect 
        CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, oldres);
        // cout << "residual " << j << " " << res << endl;
        // res2 = res = oldres;   // not used later in this function
        if (slc)
        {
          memcpy(OldU, CurrentU, N_DOF*SizeOfDouble);
          memcpy(CurrentOldDefU, CurrentDefectU, N_DOF*SizeOfDouble);
        }

        // iterate
        while (((res>TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SADDLE*oldres)
                && (res>0.1*TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE))
               || (j==0))
        {
          // apply smother
          CurrentLevel->BraessSarazin(CurrentU, CurrentDefectU, CurrentAux,
                                      N_Parameters, Parameters,N_Levels);
          j++;
          // maxit reached             
          if(j>=maxit) break;
          // compute defect
          CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, res);
          if (res > divfactor)
          {
            OutPut("Braess-Sarazin : coarse grid solver diverged " <<  res  << endl);
            exit(4711);
          }
          // cout << "residual " << j << " " << res << endl;
        }
        if (verbosity >=2)
        {
          // compute defect
          CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, res);
          OutPut("level 0: number of Braess Sarazin iters: " << j << " " << res/oldres << endl);
        }
        if (slc)
        { 
          alpha = CurrentLevel->StepLengthControl(CurrentU, OldU, CurrentOldDefU,
                                                  N_Parameters,Parameters);
      
          for (j=0;j<N_DOF;j++)
            CurrentU[j] = OldU[j] + alpha *( CurrentU[j]-OldU[j]);
        }
        break;
      case 20:
      case 21:
      case 22:
      case 30:
      case 31:
      case 32: //new vanka implementations, code is copied from case 1, case 2.
        slc = 0;
        if (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE)
          slc = 1;
        else
          if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE)
            &&(i==N_Levels-1))
	        slc = 1;
        maxit = TDatabase::ParamDB->SC_COARSE_MAXIT_SADDLE; // max # of iterations
        j = 0;
        // compute defect
        CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, oldres);
        
        if (verbosity >= 2)
        {
          OutPut("smoother " <<smoother <<" , level 0 res before smoothing: " << oldres << endl);
        }
        divfactor *= res;
        if (slc)
        {
          ii = N_DOF*SizeOfDouble;
          memcpy(OldU, CurrentU, ii);
          memcpy(CurrentOldDefU, CurrentDefectU, ii);
        }
        // iterate, do at least one iteration
        while (((res>TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SADDLE*oldres)
          || (j==0)))
        {
          // apply smother
          CurrentLevel->applySmoother(CurrentU, CurrentRhsU, CurrentAux);
          j++;
          // maxit reached
          if(j>=maxit) break;
          // compute defect
          CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, res);
          if (res > divfactor)
          {
            OutPut("mesh cell Vanka : coarse grid solver diverged " <<  res  << endl);
            exit(4711);
          }          
        }
        if (verbosity >=2)
        {
          // compute defect
          CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, res);
          OutPut("smoother " << smoother <<", level 0 res after smoothing: "<< res << " reduction " <<  res/oldres );
          OutPut(" number of Vanka iters: " << j << endl);
        }
        if (slc)
        {
          alpha = CurrentLevel->StepLengthControl(CurrentU, OldU, CurrentOldDefU,
                                                  N_Parameters,Parameters);
          
          for (j=0;j<N_DOF;j++)
            CurrentU[j] = OldU[j] + alpha *( CurrentU[j]-OldU[j]);
        }
        break;
        
      default:
        OutPut("coarse smoother not found !! Set SC_COARSE_SMOOTHER properly !! " << smoother);
        OutPut(endl);
        exit(4711);
    } // end iterative solver
  } // end coarsest grid
  // not the coarsest grid *******************************************************
  else                        
  {                            // pre smoothing
                               // check if step length control is required
    slc = 0;
    if (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE)
      slc = 1;
    else if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE)
             &&(i==N_Levels-1))
      slc = 1;
    defect_calc = 0;

    CoarserLevel = MultiGridLevels[i-1];
    ii = GEO_DIM*CoarserLevel->GetN_UDOF();
    CoarserLevel->GetSolutionVector(CoarserU);
    CoarserP = CoarserU + ii;
    CoarserLevel->GetRhsVector(CoarserRhsU);
    CoarserRhsP = CoarserRhsU + ii;

    // if step length control
    if (slc)
    {
      // compute defect
      CurrentLevel->CorrectNodes(CurrentU);
      CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, oldres);
      ii = N_DOF*SizeOfDouble;
      memcpy(OldU, CurrentU, ii);
      memcpy(CurrentOldDefU, CurrentDefectU, ii);
    }
    if (verbosity>=2)
    {
      if (!defect_calc)
      {
        CurrentLevel->CorrectNodes(CurrentU);
        CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, oldres);        
        defect_calc = 1;
      }
      OutPut("level " << i << " ");
      OutPut("res before presmoothing: " <<  
             sqrt(Ddot(GEO_DIM*N_UDOF, CurrentDefectU,CurrentDefectU))
             << " " <<  sqrt(Ddot(N_PDOF, CurrentDefectP,
                                  CurrentDefectP)) << " " << 
             oldres << endl);       
    }

    // apply the smoother 
    smoother = TDatabase::ParamDB->SC_SMOOTHER_SADDLE;
    switch(smoother)
    {
      case 1 :
      case 2 :
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SADDLE;j++)
          CurrentLevel->CellVanka(CurrentU, CurrentRhsU, CurrentAux,
                                  N_Parameters, Parameters,smoother,N_Levels);  
        break;
      case 3 :
      case 4 :  
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SADDLE;j++)
          CurrentLevel->NodalVanka(CurrentU, CurrentRhsU, CurrentAux, 
                                   N_Parameters, Parameters,smoother,N_Levels);
         break;
      case 11 :         
        if (!defect_calc)
          CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, oldres);        
        oldres = sqrt(Ddot(2*N_UDOF, CurrentU,CurrentU));
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SADDLE;j++)
        {
          CurrentLevel->BraessSarazin(CurrentU, CurrentDefectU, CurrentAux, 
                                      N_Parameters, Parameters, N_Levels);
          
          // compute defect 
          CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, res);
        }
        break;
      case 20:
      case 21:
      case 22:
      case 30:
      case 31:
      case 32:
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SADDLE;j++)
          CurrentLevel->applySmoother(CurrentU, CurrentRhsU, CurrentAux);
        break;
      default :
        OutPut("smoother not found !! Set SC_SMOOTHER properly !!"<< endl);
        Error("smoother not found !! Set SC_SMOOTHER properly !!"<< endl);
        exit(4711);
    }
     
    // calculate defect
    CurrentLevel->CorrectNodes(CurrentU);
    if (!slc)
      CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, oldres);
    else
    {
      alpha = CurrentLevel->StepLengthControl(CurrentU, OldU, 
                                              CurrentOldDefU, 
                                              N_Parameters,Parameters);       
      
      for (j=0;j< N_DOF;j++)
        CurrentU[j] = OldU[j] + alpha *( CurrentU[j]-OldU[j]);
      CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, oldres);
      ii = N_DOF*SizeOfDouble;
      memcpy(OldU, CurrentU, ii);
      memcpy(CurrentOldDefU, CurrentDefectU, ii);
    }

    if (verbosity>=2)
    {
        OutPut("level " << i << " ");
        OutPut("res after presmoothing: " <<  
               sqrt(Ddot(GEO_DIM*N_UDOF, CurrentDefectU,CurrentDefectU))
               << " " <<  sqrt(Ddot(N_PDOF, CurrentDefectP,
                                    CurrentDefectP)) << " " << 
               oldres << endl);
    }

    // calculate u-defect from u*-defect
    if(CurrentLevel->GetType() == 5)
    {
      // if you want this to work, implement and document (!) it.
      ErrThrow("unsupported NSTYPE 5");
    }

    // restrict defect
    DefectRestriction(USpaces[i-1], USpaces[i], GEO_DIM,
                   CoarserRhsU, CurrentDefectU, CurrentAux);
    DefectRestriction(PSpaces[i-1], PSpaces[i], 
                   CoarserRhsP, CurrentDefectP, CurrentAux);

    // calculate u*-representation from u-representation
    if(CoarserLevel->GetType() == 5)
    {
      // if you want this to work, implement and document (!) it.
      ErrThrow("unsupported NSTYPE 5");
    }

    CoarserLevel->CorrectDefect(CoarserRhsU);

    CoarserLevel->Reset(CoarserU);

    // coarse grid correction, apply mg recursively
    for(j=0;j<mg_recursions[i];j++)
      Cycle(i-1, res);
    /* F--cycle */
    if (TDatabase::ParamDB->SC_MG_CYCLE_SADDLE<1) mg_recursions[i] = 1;

    // post smoothing

     // calculate u from u*
    if(CoarserLevel->GetType() == 5)
    {
      // if you want this to work, implement and document (!) it.
      ErrThrow("unsupported NSTYPE 5");
    }

    // prolongate correction
    Prolongate(USpaces[i-1], USpaces[i], GEO_DIM,
                   CoarserU, CurrentDefectU, CurrentAux);
    Prolongate(PSpaces[i-1], PSpaces[i],  
                   CoarserP, CurrentDefectP, CurrentAux);

    // calculate u*-representation from u-representation
    if(CurrentLevel->GetType() == 5)
    {
      // if you want this to work, implement and document (!) it.
      ErrThrow("unsupported NSTYPE 5");
    }

    // update dofs
    CurrentLevel->CorrectNodes(CurrentDefectU);
    CurrentLevel->Update(CurrentU, CurrentDefectU);
    // apply step length control for prolongation
    if (slc)
    {
      alpha = CurrentLevel->StepLengthControl(CurrentU, OldU, 
                                              CurrentOldDefU, 
                                              N_Parameters,Parameters);       
      //OutPut("slc " <<alpha << endl);
      //alpha = 1;
      for (j=0;j<N_DOF;j++)
        CurrentU[j] = OldU[j] + alpha *( CurrentU[j]-OldU[j]);
      CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, oldres);
      memcpy(OldU, CurrentU, N_DOF*SizeOfDouble);
      memcpy(CurrentOldDefU, CurrentDefectU, N_DOF*SizeOfDouble);
    }

    defect_calc = 0;
 
    if (verbosity>=2)
    {
      // compute defect
      CurrentLevel->CorrectNodes(CurrentU);
      CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, oldres);
      defect_calc = 1;
      OutPut("level " << i << " ");
      OutPut("res before postsmoothing: " <<  
             sqrt(Ddot(GEO_DIM*N_UDOF, CurrentDefectU,CurrentDefectU))
             //sqrt(Ddot(N_UDOF, CurrentDefectU+2*N_UDOF,CurrentDefectU+2*N_UDOF))
             << " " <<  sqrt(Ddot(N_PDOF, CurrentDefectP,CurrentDefectP)) << " " << 
             oldres << endl);       
    }

    // apply the smoother 
    switch(smoother)
    {
      case 1 :
      case 2 :
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SADDLE;j++)
          CurrentLevel->CellVanka(CurrentU, CurrentRhsU, CurrentAux, 
                N_Parameters, Parameters, smoother, N_Levels);         
        break;
      case 3 :
      case 4 :  
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SADDLE;j++)
           CurrentLevel->NodalVanka(CurrentU, CurrentRhsU, CurrentAux, 
                                    N_Parameters, Parameters,smoother,N_Levels);         
        break;
      case 11 :  
        // compute defect
        if (!defect_calc) 
          CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, oldres);
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SADDLE;j++)
        {
          CurrentLevel->BraessSarazin(CurrentU, CurrentDefectU, CurrentAux, 
                                      N_Parameters, Parameters, N_Levels);
          
          // compute defect 
          CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, res);
        }
        break;
      case 20:
      case 21:
      case 22:
      case 30:
      case 31:
      case 32:
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SADDLE;j++)
          CurrentLevel->applySmoother(CurrentU, CurrentRhsU, CurrentAux);
        break;
      default :
        OutPut("smoother not found !! Set SC_SMOOTHER properly !!"<< endl);
        Error("smoother not found !! Set SC_SMOOTHER properly !!"<< endl);
        exit(4711);
    }
 
    // apply step length control
    if (slc)
    {
      alpha = CurrentLevel->StepLengthControl(CurrentU, OldU, 
                                              CurrentOldDefU, 
                                              N_Parameters,Parameters);       
      
      for (j=0;j<N_DOF;j++)
        CurrentU[j] = OldU[j] + alpha *( CurrentU[j]-OldU[j]);
    }

    if (verbosity>=2)
    {
      // compute defect
      CurrentLevel->CorrectNodes(CurrentU);
      CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, res);
      OutPut("level " << i << " ");
      OutPut("res after postsmoothing: " <<  
             sqrt(Ddot(GEO_DIM*N_UDOF, CurrentDefectU,CurrentDefectU))
             << " " <<  sqrt(Ddot(N_PDOF, CurrentDefectP,CurrentDefectP)) << " " << 
             res << endl);       
    }
  }
}

void TNSE_MultiGrid::SetDirichletNodes(int i)
{
  TNSE_MGLevel *CurrentLevel;
  double *CurrentU;
  double *CurrentRhsU;
  int HangingNodeBound, N_Dirichlet, N_UDOF, j;

  if(i>=N_Levels) return;

  CurrentLevel = MultiGridLevels[i];
  N_UDOF = CurrentLevel->GetN_UDOF();
  CurrentLevel->GetSolutionVector(CurrentU);
  CurrentLevel->GetRhsVector(CurrentRhsU);

  HangingNodeBound = CurrentLevel->GetHangingNodeBound();
  N_Dirichlet = CurrentLevel->GetN_Dirichlet();

  for (j=0;j<GEO_DIM;j++)
    memcpy(CurrentU+j*N_UDOF+HangingNodeBound, 
           CurrentRhsU+j*N_UDOF+HangingNodeBound,
           SizeOfDouble*N_Dirichlet);
}

/** return residual on grid i */
double TNSE_MultiGrid::GetResidual(int i)
{
  TNSE_MGLevel *CurrentLevel;
  double *CurrentU;
  double *CurrentRhsU;
  double *CurrentDefectU;
  double oldres;
  
  CurrentLevel = MultiGridLevels[i];

  CurrentLevel->GetSolutionVector(CurrentU);
  CurrentLevel->GetRhsVector(CurrentRhsU);
  CurrentDefectU = CurrentLevel->GetAuxVector(0);

  CurrentLevel->Defect(CurrentU, CurrentRhsU, CurrentDefectU, oldres);

  return oldres;
}

/** set recursion for multigrid */ 
void TNSE_MultiGrid::SetRecursion(int levels)
{
  int gam = TDatabase::ParamDB->SC_MG_CYCLE_SADDLE,k;

  // coarsest grid 
  mg_recursions[1] = 1;
  if (gam>0)
    for (k=2;k<=levels;k++)        
      mg_recursions[k] = gam;
  else                /* F -- cycle */
    for (k=2;k<=levels;k++)        
      mg_recursions[k] = 2;
}
