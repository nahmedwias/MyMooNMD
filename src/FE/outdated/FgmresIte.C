// =======================================================================
// @(#)FgmresIte.C        1.24 06/27/00
//
// Class:       TFgmresIte
// Purpose:     flexible gmres
//
// Author:      Volker John 24.10.2000
//
// History:     24.10.2000 start of implementation
//
// =======================================================================
#include <ItMethod.h>
#include <FgmresIte.h>
#include <MooNMD_Io.h>
#include <Database.h>
#include <LinAlg.h>
#include <NSE_MultiGrid.h>

#ifdef _MPI 
#include <ParFECommunicator3D.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/** constructor with initialization */
TFgmresIte::TFgmresIte(
MatVecProc *MatVec,                               // a function to perform Mat*Vec
DefectProc *Defect,                               // a function to compute residuals
TItMethod *Prec,                                  // a preconditioner
int n_aux,                                        // ?
int n_dof,                                        // dof
int scalar                                        // flag for scalar problems
)
: TItMethod(MatVec, Defect, Prec, n_aux, n_dof)
{
  int i;
  double *aux;

  matvec = MatVec;
  matvecdefect = Defect;
  prec = Prec;
  N_DOF = n_dof;

  if (scalar)
  {
    res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
    red_factor= TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
    maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
  }
  else
  {
    res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE;
    red_factor= TDatabase::ParamDB->SC_LIN_RED_FACTOR_SADDLE;
    maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SADDLE;
  }

  prec_maxit = TDatabase::ParamDB->SC_AMG_PREC_IT;// max. no of iteration of prec.
  div_factor = TDatabase::ParamDB->SC_DIV_FACTOR;
  minit = TDatabase::ParamDB->SC_MINIT;
  restart = TDatabase::ParamDB->SC_GMRES_RESTART; // restart gmres

  if (restart > maxit)
    restart = maxit;
  if (restart <=0 )
  {
    OutPut("WARNING: restart too small: " << restart);
    OutPut("         Set restart to 10 !!!");
    restart = 10;
  }
  /* if (restart > max_gmres_restart)
  {
    restart = MAX_GMRES_RESTART;
    TDatabase::ParamDB->SC_GMRES_RESTART =  MAX_GMRES_RESTART;
    OutPut("WARNING: restart greater than MAX_GMRES_RESTART !!!" << endl);
    OutPut("WARNING: set restart to " << MAX_GMRES_RESTART << endl);
    }*/

  if (prec_maxit<1)
  {
    OutPut("WARNING: Number of preconditioner iterations too small: "
      << prec_maxit << endl);
    OutPut("         Set number of preconditioner iterations to 1 !!!");
    prec_maxit = TDatabase::ParamDB->SC_AMG_PREC_IT = 1;
  }

  s = new double[restart+1];
  cosi = new double[restart+1];
  sn =  new double[restart+1];

  H = new double*[restart+1];
  H[0] = new double[(restart+1)*restart];
  memset(H[0],0,((restart+1)*restart)*SizeOfDouble);
  for(i=1;i<=restart;i++)
    H[i] = H[0] + i*restart;

  v = new double*[restart+1];
  v[0] = new double[(restart+1)*N_DOF];
  memset(v[0],0,((restart+1)*N_DOF)*SizeOfDouble);
  for(i=1;i<=restart;i++)
    v[i] = v[0] + i*N_DOF;

  zv = new double*[restart+1];
  zv[0] = new double[(restart+1)*N_DOF];
  memset(zv[0],0,((restart+1)*N_DOF)*SizeOfDouble);
  for(i=1;i<=restart;i++)
    zv[i] = zv[0] + i*N_DOF;

  defect =  new double[N_DOF];

  N_Aux = n_aux;
  if (n_aux>0)
  {
    AuxArray = new double* [n_aux];
    aux = new double[n_aux*N_DOF];
    for(i=0;i<n_aux;i++)
      AuxArray[i] = aux+i*N_DOF;
  }
}


#ifdef _MPI
TFgmresIte::TFgmresIte(MatVecProc *MatVec,
			DefectProc *Defect,
			TItMethod *Prec,
			int n_aux, int n_dof,
			int scalar,
  	#ifdef  __3D__
			       TParFECommunicator3D *parcomm
	#else		       
				TParFECommunicator2D *parcomm
	#endif 
                      )
: TItMethod(MatVec, Defect, Prec, n_aux, n_dof)
{
  int i;
  double *aux;

  matvec = MatVec;
  matvecdefect = Defect;
  prec = Prec;
  N_DOF = n_dof;

  if (scalar)
  {
    res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
    red_factor= TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
    maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
  }
  else
  {
    res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE;
    red_factor= TDatabase::ParamDB->SC_LIN_RED_FACTOR_SADDLE;
    maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SADDLE;
  }

  prec_maxit = TDatabase::ParamDB->SC_AMG_PREC_IT;
  div_factor = TDatabase::ParamDB->SC_DIV_FACTOR;
  minit = TDatabase::ParamDB->SC_MINIT;
  restart = TDatabase::ParamDB->SC_GMRES_RESTART;

  if (restart > maxit)
    restart = maxit;
  if (restart <=0 )
  {
    OutPut("WARNING: restart too small: " << restart);
    OutPut("         Set restart to 10 !!!");
    restart = 10;
  }
  /* if (restart > max_gmres_restart)
  {
    restart = MAX_GMRES_RESTART;
    TDatabase::ParamDB->SC_GMRES_RESTART =  MAX_GMRES_RESTART;
    OutPut("WARNING: restart greater than MAX_GMRES_RESTART !!!" << endl);
    OutPut("WARNING: set restart to " << MAX_GMRES_RESTART << endl);
    }*/

  if (prec_maxit<1)
  {
    OutPut("WARNING: Number of preconditioner iterations too small: "
      << prec_maxit << endl);
    OutPut("         Set number of preconditioner iterations to 1 !!!");
    prec_maxit = TDatabase::ParamDB->SC_AMG_PREC_IT = 1;
  }

  s = new double[restart+1];
  cosi = new double[restart+1];
  sn =  new double[restart+1];

  H = new double*[restart+1];
  H[0] = new double[(restart+1)*restart];
  memset(H[0],0,((restart+1)*restart)*SizeOfDouble);
  for(i=1;i<=restart;i++)
    H[i] = H[0] + i*restart;

  v = new double*[restart+1];
  v[0] = new double[(restart+1)*N_DOF];
  memset(v[0],0,((restart+1)*N_DOF)*SizeOfDouble);
  for(i=1;i<=restart;i++)
    v[i] = v[0] + i*N_DOF;

  zv = new double*[restart+1];
  zv[0] = new double[(restart+1)*N_DOF];
  memset(zv[0],0,((restart+1)*N_DOF)*SizeOfDouble);
  for(i=1;i<=restart;i++)
    zv[i] = zv[0] + i*N_DOF;

  defect =  new double[N_DOF];

  N_Aux = n_aux;
  if (n_aux>0)
  {
    AuxArray = new double* [n_aux];
    aux = new double[n_aux*N_DOF];
    for(i=0;i<n_aux;i++)
      AuxArray[i] = aux+i*N_DOF;
  }
  
   ParComm = parcomm;
}
#endif


TFgmresIte::~TFgmresIte()
{
  delete s;
  delete cosi;
  delete sn;
  delete H[0];
  delete H;
  delete v[0];
  delete v;
  delete zv[0];
  delete zv;
  delete defect;

  if (N_Aux>0)
  {
    delete AuxArray[0];
    delete AuxArray;
  }
}


void GeneratePlaneRotation(double dx, double dy, double *cs, double *sn)
{
  if (dy == 0.0)
  {
    cs[0] = 1.0;
    sn[0] = 0.0;
  }
  else if (fabs(dy) > fabs(dx))
  {
    double temp = dx / dy;
    sn[0] = 1.0 / sqrt( 1.0 + temp*temp );
    cs[0] = temp * sn[0];
  }
  else
  {
    double temp = dy / dx;
    cs[0] = 1.0 / sqrt( 1.0 + temp*temp );
    sn[0] = temp * cs[0];
  }
}

// rotate a vector by theta (cs = cos(theta), sn= sin(theta)
void ApplyPlaneRotation(double *dx, double *dy, double cs, double sn)
{
  double temp  =  cs * dx[0] + sn * dy[0];
  dy[0] = -sn * dx[0] + cs * dy[0];
  dx[0] = temp;
}


void UpdateGmresIterate(double *x, int Len_x, int k, double **h,
double *s, double **v)
{
  int i,j;
  /* Backsolve: */

  for (i = k; i >= 0; i--)
  {
    s[i] /= h[i][i];
    for (j = i - 1; j >= 0; j--)
      s[j] -= h[j][i] * s[i];
  }

  for (j = 0; j <= k; j++)
    for(i=0; i<Len_x; i++)
      x[i] += v[j][i] * s[j];
}


int TFgmresIte::Iterate (TSquareMatrix **sqmat,
TMatrix **mat, double *sol,
double *rhs)
{
  Output::warn("TFgmresIte::Iterate", "Use new FGMRES implementation instead.");
  int verbose = 1;
  int i=0, j,k,l;
  int maxite;
  double res, res0, reslast, t1, t2; //temp,tempGlobal;
  double beta,residlast,dnorm0; //end_residual;
//  double eps = 1e-14;
//#ifdef _MPI
//  const int* MasterOfDof;
//  int ii, rank;
//  double  resglobal,temp,tempGlobal;
//
//   MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
//   MasterOfDof = ParComm->GetMaster();
//   ParComm->CommUpdate(sol, rhs);
//#endif
  int flexible = TDatabase::ParamDB->SC_FLEXIBLE_KRYLOV_SPACE_SOLVER;
  if(flexible)
  {
    t1 = GetTime();
    if (verbose>1
//#ifdef _MPI
//    && rank==0
//#endif
    )
      OutPut("Entering fgmres" << endl);
    
    if ((TDatabase::ParamDB->INTERNAL_GMRES_INFO)&&(verbose==-1))
    {
      OutPut("fgmres maxite " << maxit << " restart " << restart << endl);
      TDatabase::ParamDB->INTERNAL_GMRES_INFO = 0;
    }

    // compute defect, result in v[0]
    matvecdefect(sqmat, mat, sol, rhs, v[0]);
    
//#ifdef _MPI
//    ParComm->CommUpdate(v[0], v[0]);
//    res=0.;
//    for(ii=0; ii<N_DOF; ii++)
//      if(MasterOfDof[ii] == rank)
//        res += v[0][ii]*v[0][ii];
//
//    MPI_Allreduce(&res, &resglobal, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
//    res=res0=reslast = sqrt(resglobal);
//#else
    
    res=(Ddot(N_DOF,v[0],v[0]));    
    res=res0=reslast=sqrt(Ddot(N_DOF,v[0],v[0]));     
     
//#endif
                                 /* norm of residual */
    beta=res;
 

    //res = res0 = beta;
    maxite = maxit;
    // residuals[0]=resid;
    // residual_cnt++;
    //OutPut("beta " << beta << endl);
    
    if (beta <= res_norm_min )     /* stopping criterion fulfilled */
    {
        //if ((minit==0)||(beta<eps))
      if ((minit==0)||(beta<1e-20))
      {
        if (verbose>0)
          OutPut("(no) fgmres iteration " << 0 << " " << beta << endl);
        // iteration_cnt=0;
//        end_residual=res;
        return(0);
      }
      else
      {
        maxite = minit;
      }
    }
    
    if (verbose>0
//#ifdef _MPI
//    && rank==0
//#endif
    )
      OutPut("fgmres Iteration " << 0 << " " << beta << endl);
    
    j=1;
    // iterate
    while (j<=maxite)
    {
      Dscal(N_DOF, 1.0/beta, v[0]);
      memset(s,0,(restart+1)*SizeOfDouble);
      s[0]=beta;

      // iterate until converged or restart
      for(i=0; i<restart && j<=maxite; i++, j++)
      {
        if (prec==nullptr)
          memcpy(zv[i],v[i],N_DOF*SizeOfDouble);
        else
        {
          memset(zv[i],0,N_DOF*SizeOfDouble);
          if (prec_maxit>1)
            dnorm0 = sqrt(Ddot(N_DOF,v[i],v[i]))*TDatabase::ParamDB->SC_AMG_PREC_RED_FACTOR;
          // preconditioner iteration
          for (k=0;k<prec_maxit;k++)
          {
            // solution of prec ite in zv[i]
            // rhs for prec in v[i]
            prec->Iterate(sqmat, mat, zv[i], v[i]);
      
//#ifdef _MPI
//            ParComm->CommUpdate(zv[i], v[i]);
//#endif
            if (k<prec_maxit-1)
            {
              // compute defect
              matvecdefect(sqmat, mat, zv[i], v[i], defect);
//#ifdef _MPI
//              res=0.;
//              for(ii=0; ii<N_DOF; ii++)
//                if(MasterOfDof[ii] == rank)
//                  res += defect[ii]*defect[ii];
//
//              MPI_Allreduce(&res, &resglobal, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
//              res=res0=reslast = sqrt(resglobal);
//#else
   
              res=(Ddot(N_DOF,defect,defect));    
              res=res0=reslast=sqrt(Ddot(N_DOF,defect,defect));     
              
//#endif
// 	      printf("res is %lf ---------\n",res);
// #ifdef _MPI
//   MPI_Finalize();
// #endif
//   exit(0);
              if (res <= dnorm0)
                break;
            }
          }
        }

                                 /* d=A*defect; */
        matvec(sqmat, mat, zv[i], defect);
        
        //OutPut(sqrt(Ddot(N_DOF,defect,defect))<< endl);
        for(k=0; k<=i; k++)
        {
//#ifdef _MPI
//          temp=0.;
//          for(ii=0; ii<N_DOF; ii++)
//            if(MasterOfDof[ii] == rank)
//              temp += v[k][ii]*defect[ii];
//
//          MPI_Allreduce(&temp, &H[k][i], 1, MPI_DOUBLE, MPI_SUM,
//                        TDatabase::ParamDB->Comm);
//#else
          H[k][i]= Ddot(N_DOF, defect,v[k]);    
//#endif

          Daxpy(N_DOF, -H[k][i], v[k], defect);
//#ifdef _MPI
//          ParComm->CommUpdate(defect, defect);
//#endif
        }

//#ifdef _MPI
//        temp=0.;
//        for(ii=0; ii<N_DOF; ii++)
//          if(MasterOfDof[ii] == rank)
//            temp += defect[ii]*defect[ii];
//
//        MPI_Allreduce(&temp, &tempGlobal, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
//        H[i+1][i]= sqrt(tempGlobal);
//#else
        H[i+1][i]=sqrt(Ddot(N_DOF,defect,defect));
//#endif
      
     
        Dcopy(N_DOF, defect, v[i+1]);
        if(H[i+1][i]!=0)
          Dscal(N_DOF, 1/H[i+1][i], v[i+1]);
        else
        {
          Error("Fehler!" << endl);
          OutPut("Fehler in FgmresIte.C !!!" << endl);
          // exit(-1);
        }

        for(k=0;k<i;k++)
          ApplyPlaneRotation(&H[k][i], &H[k+1][i], cosi[k], sn[k]);

        GeneratePlaneRotation(H[i][i], H[i+1][i], &cosi[i], &sn[i]);
        ApplyPlaneRotation(&H[i][i], &H[i+1][i], cosi[i], sn[i]);
        ApplyPlaneRotation(&s[i], &s[i+1], cosi[i], sn[i]);

        residlast=res;

        // stopping criterion fulfilled
        if ((j >= minit) && (
          (((res=fabs(s[i+1])) < res_norm_min)
          || ((res=fabs(s[i+1]))<red_factor*res0))))
        {
          memset(defect,0,N_DOF*SizeOfDouble);
          UpdateGmresIterate(defect, N_DOF, i, H, s, zv);
                                 /* new iterate */
          Daxpy(N_DOF,1.0,defect,sol);
          t2 = GetTime();
          if (verbose>=1
//#ifdef _MPI
//          && rank==0
//#endif
          )
          {
            OutPut("FGMRES ITE: " << setw(4) << j);
            OutPut(" t: " << setw(6) << t2-t1);
            OutPut(" t/cyc: " << setw(6) << (t2-t1)/j);
            OutPut(" res : "  << setw(8) << res);
            OutPut(" rate: " << pow(res/res0,1.0/j) << endl);
            //iteration_cnt=j;
            //end_residual=res;
          }
          if (verbose==-2
//#ifdef _MPI
//            && rank==0
//#endif
          )
          {
            OutPut("vanka fgmres ite: " << setw(4) << j);
            OutPut(" t: " << setw(6) << t2-t1);
            OutPut(" t/cyc: " << setw(6) << (t2-t1)/j);
            OutPut(" res : "  << setw(8) << res);
            OutPut(" rate: " << pow(res/res0,1.0/j) << endl);
            //iteration_cnt=j;
            //end_residual=res;
          }
          return(j);
        }
        if (verbose>0
//#ifdef _MPI
//          && rank==0
//#endif
        )
        {
          OutPut("fgmres iteration " << j << " " << res << " " << res/residlast);
          OutPut(endl);
        }
        if (res>div_factor*res0)
        {
          OutPut("fgmres iteration diverges !!!" << endl);
          exit(4711);
        }
        //if (j<maxit)
        // {
        // residuals[residual_cnt%AMG_CONV_RATE_BACK]=res;
        //  residual_cnt++;
        //}
        //  finish = clock();
        //  if (finish>=start)
        //    elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
        //  else
        //    elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;
        //  start = clock();
      }                            /* endfor i */

      /* compute the update=rhs for prec. */
                                 /* set solution to zero */
      memset(v[0],0, N_DOF*SizeOfDouble);
      UpdateGmresIterate(v[0], N_DOF, restart-1, H, s, zv);

      Daxpy(N_DOF,1.0,v[0],sol);   /* new iterate */
                                  /* new defect v[0] = b - A*x */

      matvecdefect(sqmat, mat, sol, rhs, v[0]);
//#ifdef _MPI
//      ParComm->CommUpdate(v[0], rhs);
//      // ParComm->CommUpdate(sol, rhs);
//#endif                                 /* norm of residual */
 
//#ifdef _MPI
//      res=0.;
//      for(ii=0; ii<N_DOF; ii++)
//        if(MasterOfDof[ii] == rank)
//          res += v[0][ii]*v[0][ii];
//
//      MPI_Allreduce(&res, &resglobal, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
//      res=res0=reslast = sqrt(resglobal);
//#else
   
      res=(Ddot(N_DOF,v[0],v[0]));    
      res=res0=reslast=sqrt(Ddot(N_DOF,v[0],v[0]));     
     
//#endif
      beta=res;

      if (verbose>1)
      {
        if(fabs(beta-res)>0.01*beta)
        {
          OutPut("WARNING: restart residual changed " <<  beta << " " << res);
          OutPut(endl);
        }
      }
    }                              /* endwhile */

    //iteration_cnt=j-1;
    //end_residual=res;

    t2 = GetTime();
    // iteration stopped because of reaching maximal iteration number
    if (j>=maxite && res>res_norm_min &&res>res0*red_factor )
    {
      if (verbose>-1)
        OutPut("FGMRES not converged !!!" << endl);
      //OutPut("FGMRES Iteration: (maximal) iterations " << i-1 <<
      //     " residual " << res << endl);
      //    iteration_cnt=i;
      //    end_residual=res;
      //return(0);
    }

    if (j>maxite)
      j=maxite;
    if (verbose>-1)
    {
      OutPut("FGMRES ITE: " << setw(4) << j);
      OutPut(" t: " << setw(6) << t2-t1);
      OutPut(" t/cyc: " << setw(6) << (t2-t1)/j);
      OutPut(" res : "  << setw(8) << res);
      OutPut(" rate: " << pow(res/res0,1.0/j) << endl);
    }
    if (verbose==-2)
    {
      OutPut("vanka fgmres ite: " << setw(4) << j);
      OutPut(" t: " << setw(6) << t2-t1);
      OutPut(" t/cyc: " << setw(6) << (t2-t1)/j);
      OutPut(" res : "  << setw(8) << res);
      OutPut(" rate: " << pow(res/res0,1.0/j) << endl);
    }
    // iteration_cnt=i;
    // end_residual=res;
    return(j);
  }
  else // not flexible
  {
    double *r = new double[N_DOF];
    if (verbose>1)
      OutPut("Entering GMRES" <<endl);
    //AMG_dcopy(r[0],b);                     /* copy rhs (b) into r */
    Dcopy(N_DOF, rhs, r);
    //dmatminus(r[0],A[0],B,x);               /* r = rhs - A*x, r is overwritten */
    matvecdefect(sqmat, mat, sol, rhs, r);

    //beta=sqrt(AMG_ddot(r[0],r[0]));           /* norm of residual */
    beta=sqrt(Ddot(N_DOF, r,r));
    double resid = beta;
    double start_residual=resid;
    if ((beta  <= res_norm_min )&&(minit==0))     /* stopping criterion fulfilled */
    {
      if (verbose>0)
      {
        OutPut("GMRES Iteration 0 " << beta << endl);
//        end_residual=resid;
      }
      return(0);
    }
    if (verbose>0)
    {
      OutPut("GMRES Iteration 0 " << beta << endl);
    }
    if (verbose>1)
    {
      OutPut("Entering GMRES iteration cycle" << endl);
    }
    j=1;
    while (j<=maxit)
    {
      Dcopy(N_DOF, r, v[0]);
      Dscal(N_DOF, 1.0/beta, v[0]);
      memset(s,0,(restart+1)*SizeOfDouble);
      s[0]=beta;
      for (i=0;i<restart && j<=maxit;i++, j++)
      {
        if (prec==nullptr)
          /* d=A*v; */
          matvec(sqmat, mat, v[i], defect);
        else
        {
          Dcopy(N_DOF, v[i], zv[0]);
          memset(r,0.0,N_DOF*SizeOfDouble);
          // preconditioner iteration
          for (k=0;k<prec_maxit;k++)
          {
            // solution of prec ite in zv[i]
            // rhs for prec in v[i]
            prec->Iterate(sqmat, mat, r, v[i]);
          }
          matvec(sqmat, mat, r, defect);
        }

        //TODO teste norm defect
        //     OutPut("norm v[i] " << sqrt(Ddot(N_DOF, v[i], v[i])) << endl);
        //     OutPut("norm r " << sqrt(Ddot(N_DOF, r, r)) << endl);
        //     OutPut("norm d[0] " << sqrt(Ddot(N_DOF, defect, defect)) << endl);
        for(k=0; k<=i; k++)
        {
          H[k][i]=Ddot(N_DOF, defect, v[k]);
          Daxpy(N_DOF, -H[k][i], v[k], defect);
        }
        H[i+1][i]=sqrt(Ddot(N_DOF, defect,defect));
        Dcopy(N_DOF, defect, v[i+1]);
        Dscal(N_DOF, 1.0/H[i+1][i], v[i+1]);
        for(k=0;k<i;k++)
          ApplyPlaneRotation(&H[k][i], &H[k+1][i], cosi[k], sn[k]);
        GeneratePlaneRotation(H[i][i], H[i+1][i], &cosi[i], &sn[i]);
        ApplyPlaneRotation(&H[i][i], &H[i+1][i], cosi[i], sn[i]);
        ApplyPlaneRotation(&s[i], &s[i+1], cosi[i], sn[i]);
        residlast=resid;
        if (((resid=fabs(s[i+1])) < res_norm_min)
          || ((resid=fabs(s[i+1]))<red_factor*start_residual))
        {
                                                  /* set solution to zero */
          memset(zv[0],0,((restart+1)*N_DOF)*SizeOfDouble);
          UpdateGmresIterate(zv[0],N_DOF, i, H, s, v);
          //TODO
          //  OutPut("norm z[0] " << sqrt(Ddot(N_DOF, zv[0], zv[0])) << endl);
          if (prec==nullptr)                         /* preconditioner should be applied */
          {
            Daxpy(N_DOF, 1.0, zv[0], sol);
          }
          else
          {
            memset(r,0.0,N_DOF*SizeOfDouble);
            // preconditioner iteration
            for (k=0;k<prec_maxit;k++)
            {
              // solution of prec ite in zv[i]
              // rhs for prec in v[i]
              prec->Iterate(sqmat, mat, r, zv[0]);
            }
            //TODO
            //    OutPut("norm sol " << sqrt(Ddot(N_DOF, sol, sol)) << endl);
            //    OutPut("norm r " << sqrt(Ddot(N_DOF, r, r)) << endl);

            Daxpy(N_DOF, 1.0, r, sol);            /* new iterate */
          }
          //TODO
          //    OutPut("norm sol " << sqrt(Ddot(N_DOF, sol, sol)) << endl);
          if (verbose>0)
          {
            OutPut("GMRES (right) : iterations " << j << " residual " <<resid << endl);
          }
//          end_residual=resid;
          return(j);
        }
        if (verbose>0)
        {
          OutPut("GMRES (right) Iteration " << j << " residuum " << resid << " " << resid/residlast << endl);
        }
        if (resid>div_factor*start_residual)
        {
          OutPut("GMRES (right) iteration diverges !!!" << endl);
          exit(4711);
        }
        //            if (j<maxit)
        //            {
        //                residuals[residual_cnt%AMG_CONV_RATE_BACK]=resid;
        //                residual_cnt++;
        //            }
        //            finish = clock();
        //            if (finish>=start)
        //              elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
        //            else
        //              elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;
        //            start = clock();
      }                                           /* endfor i */
      /* compute the update=rhs for prec. */
      memset(zv[0],0,N_DOF*SizeOfDouble);
      //TODO                       /* set solution to zero */
      UpdateGmresIterate(zv[0], N_DOF, restart-1, H, s, v);
      if (prec==nullptr)                             /* preconditioner should be applied */
      {
        Dcopy(N_DOF, zv[0],r);
      }
      else
      {
        memset(r,0,N_DOF*SizeOfDouble);           /* set solution to zero */
        for (l=0;l<prec_maxit;l++)                /* number of preconditioner iterations */
        {                                         /* apply preconditioner */
          prec->Iterate(sqmat, mat, r, zv[0]);
        }
      }
      Daxpy(N_DOF, 1.0, r, sol);                  /* new iterate */
      Dcopy(N_DOF,rhs,r);                         /* copy rhs (b) into r */
      /* new defect r = r - A*x, r is overwritten */
      matvecdefect(sqmat, mat, sol, rhs, r);
      beta=sqrt(Ddot(N_DOF, r, r));               /* norm of residual */

      if (verbose>1)
      {
        if(fabs(beta-resid)>0.01*beta)
          OutPut("restart residual changed " << beta << " " << resid << endl);
      }
    }
    /* endwhile */
    if (verbose>0)
    {
      OutPut("GMRES (right) : (maximal) iterations " << maxit << " residual " << resid << endl);
    }
//    end_residual=resid;
    return(j);
  }
}
