#include <Database.h>
#include <LinAlg.h>
#ifdef __2D__
#include <FEFunction2D.h>
#elif __3D__
#include <FEFunction3D.h>
#endif
#include <AlgebraicFluxCorrection.h>
#include <IsoBoundEdge.h>
#include <BoundComp.h>
#include <BaseCell.h>

#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <MooNMD_Io.h>

#ifdef _MPI
#include <mpi.h>
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

//! Anonymous namespace for helper methods which do not have to be assessed
//! from the outside.
namespace
{
  /** Check whether the setup fulfils the cfl-like condition needed in a partly
   * explicit afc scheme (Kuzmin 2009, eq (11) ).
   * Print a warning if not so.
   * @param tau Current timestep length.
   * @param theta2 Is expliciteness parameter.
   */
  void check_cfl_condition(double mass_diag_entry,
    double system_diag_entry, double tau, double theta2)
  {
    //check CFL
    double cfl = mass_diag_entry/system_diag_entry;
    if (theta2 >0)
    {
      cfl /= theta2;
      if (tau > cfl)
      {
        Output::print("WARNING in FEM-FCT system matrix: ",
          mass_diag_entry, " cfl violated: cfl ",
          cfl, " delta t: ", tau);
      }
    }
  }
/** ************************************************************************ */
  
  /**
   * Compute the entire system matrix of a FEM-FCT corrected system, after
   * its mass matrix got lumped and its stiffness matrix flux corrected.
   *
   * Also performs checking of the cfl condition in case an explicit
   * method has been used as predictor ("intermediate solution").
   *
   * @note Makes sense only for (at least partly) implicit schemes.
   * @param[in, out] system_matrix The systems stiffness matrix A, with additional
   * diffusion already added to it. Must be square. Will be turned to
   * S = M_(\text{lumped}) + \Delta_t * \theta_1 * S.
   * @param[in] lumped_mass_matrix The row wise lumped mass matrix of the system,
   * as a vector of its diagonal entries. Must have length equal to the stiffness
   * matrix' dimension.
   * @param[in] tau The current time step length.
   * @param theta1 The impliciteness parameter. Must be 0.5 so far
   * (and defaults to that value).
   * @param theta2 The expliciteness parameter. Must be 0.5 so far
   * (and defaults to that value).
   */
  void fem_fct_compute_system_matrix(
    FEMatrix& system_matrix,
    const std::vector<double>& lumped_mass_matrix,
    double tau, double theta1 = 0.5, double theta2 = 0.5)
  {
    if (theta1 != 0.5 || theta2 != 0.5  )
    {
      ErrThrow("fem_fct_compute_system_matrix for Crank Nicolson only."
        "TODO Implement other schemes!");
    }
    // check if stiffness matrix is square
    if(!system_matrix.is_square())
    {
      ErrThrow("system_matrix must be square for FEM-FCT!");
    }

    int nDof = system_matrix.GetN_Rows();

    if((int) lumped_mass_matrix.size() != nDof)
    {
      ErrThrow("Lumped mass matrix vector does not fit the stiffness matrix!");
    }

#ifdef _MPI
    const TParFECommunicator3D& comm = system_matrix.GetFESpace3D()->get_communicator();
    const int* masters = comm.GetMaster();
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    // get pointers to columns, rows and entries of system_matrix
    const int* ColInd = system_matrix.GetKCol();
    const int* RowPtr = system_matrix.GetRowPtr();
    double*  Entries = system_matrix.GetEntries();

    for(int i=0;i<nDof;i++)
    {
#ifdef _MPI
      // do this only for master rows
      if(masters[i] != rank)
        continue;
#endif
      // i-th row of sqmatrix
      int j0 = RowPtr[i];
      int j1 = RowPtr[i+1];

      for( int j = j0 ; j < j1 ; j++ )
      {
        int index = ColInd[j];
        if( i == index )
        {                                         //diagonal entry

          check_cfl_condition(lumped_mass_matrix[i], Entries[j], tau, theta2);

          // calculate the new system matrix entry
          Entries[j] = lumped_mass_matrix[i] + tau * theta1 * Entries[j];

        }
        else
        {                                         //non-diagonal entry
          // calculate the system matrix entry
          Entries[j]= tau * theta1 * Entries[j];

          //print warning if S violates M-matrix necessary condition (a lot!)
          if ( Entries[j]>1e-10)
          {
            Output::print("WARNING in FEM-FCT system matrix:  ",
              Entries[j], " is a positive off-diagonal entry!");
          }
        }
      }
    }
  }
/** ************************************************************************ */
  
  /*!
   * Compute the artificial diffusion matrix to a given stiffness matrix.
   *
   * @param[in] A The stiffness matrix which needs additional
   * diffusion. Must be square.
   * @param[out] matrix_D A vector to write the entries
   * of D, the artificial diffusion matrix to. Must contain only 0.
   * MPI: D will leave this method with matrix-consistency level 0, meaning that
   * only its master rows will be correct. Everything else is left at 0.
   *
   */
  void compute_artificial_diffusion_matrix(
    const FEMatrix& A, FEMatrix& D)
  {
    // catch non-square matrix
    if (!A.is_square() )
    {
      ErrThrow("Matrix must be square!");
    }
    if(A.GetStructure()!=D.GetStructure())
    {
      ErrThrow("A and D should have the same structure!");
    }
    // store number of dofs
    int nDof = A.GetN_Rows();
    double* D_Entries=D.GetEntries();

#ifdef _MPI
    const TParFECommunicator3D& comm = A.GetFESpace3D()->get_communicator();
    const int* masters = comm.GetMaster();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    // get pointers to columns, rows and entries of matrix A
    const int* ColInd = A.GetKCol();
    const int* RowPtr = A.GetRowPtr();
    const double* Entries = A.GetEntries();

    // compute off-diagonal entries of the matrix D
    for(int i=0;i<nDof;i++)                       //row loop
    {
#ifdef _MPI
      // do this only for master rows - they couple only to other masters,
      // slaves and halo1 d.o.f., therefore one can be sure, that if A_ij
      // is on this rank, so is the transposed entry, A_ji
      if(masters[i] != rank)
        continue;
#endif
      for(int l=RowPtr[i];l<RowPtr[i+1];l++)
      {
        int j = ColInd[l];                        //column index
        double k_ij = 0;
        double k_ji = 0;
        if (j!=i)
        {                                         // consider only off-diagonals
          k_ij = Entries[l];
          // now get the transposed entry
          for (int ll=RowPtr[j];ll<RowPtr[j+1];ll++)
          {
            if (ColInd[ll]==i)
            {
              k_ji = Entries[ll];
              break;
            }
          }
          //determine entry D_ij
          D_Entries[l] = std::min({-k_ij , 0.0 , -k_ji}
          );
        }
      }
    }

    // compute diagonal entries of the matrix D
    for(int i=0;i<nDof;i++)                       //row loop
    {
#ifdef _MPI
      // do this only for master rows
      if(masters[i] != rank)
        continue;
#endif
      double val = 0.0;
      // add all entries of i-th row
      int ll = -1;
      for(int l=RowPtr[i];l<RowPtr[i+1];l++)
      {
        val +=  D_Entries[l];
        int j = ColInd[l];
        if (j==i)                                 //diagonal found
          ll = l;                                 //place of the diagonal entry in the matrix_D entries array
      }
      D_Entries[ll] = -val;
    }
  }

/** ************************************************************************ */
  /**
   * @brief Lump a mass matrix row wise.
   *
   * @param[in] M The mass matrix to be lumped. Must be square.
   *
   * @param[out] lump_mass A vector containing the diagonal entries
   *  of the lumped mass matrix. Must be of according length.
   */
  void lump_mass_matrix_to_vector(
    const FEMatrix& M, std::vector<double>& lump_mass)
  {
    // catch non-square matrix
    if (!M.is_square() )
    {
      ErrThrow("Mass matrix must be square!");
    }
    // catch not fitting vector size
    if ( (int) lump_mass.size() != M.GetN_Columns())
    {
      ErrThrow("Vector to write lumped mass matrix to has incorrect size! ",
        lump_mass.size(), " != ", M.GetN_Columns());
    }

    const int * RowPtr = M.GetRowPtr();
    const double * Entries = M.GetEntries();
    int rows = M.GetN_Rows();

#ifdef _MPI
    const TParFECommunicator3D& comm = M.GetFESpace3D()->get_communicator();
    const int* masters = comm.GetMaster();
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    for (int i=0; i<rows; i++)
    {
#ifdef _MPI
      //skip the non-masters in MPI case
      if(masters[i] != rank)
        continue;
#endif
      lump_mass[i]=0.0;
      int j0 = RowPtr[i];
      int j1 = RowPtr[i+1];

      for (int j=j0;j<j1;j++)
      {
        lump_mass[i] += Entries[j];
      }
      if(lump_mass[i]==0)
      {
        ErrThrow("zero entry in lumped matrix ", i, " ", lump_mass[i]);
      }
    }
  }

  /**
   * @brief The MinMod Flux Prelimiter.
   * @param a first value
   * @param b second value
   * @return The result of the MinMod Prelimiter.
   */
  double MinMod(double a, double b)
  {
    if (a*b < 0)
    {
      return 0.0;
    }
    //CB 2015/07/27 should not this return a negative value, when a<0 and b<0?
    // int sign = a < 0 ? -1 : 1;

    if (fabs(a) <  fabs(b))
    {
      return a;                                   //return sign * a;
    }
    else
    {
      return b;                                   //return sign * b;
    }
  }
  
/** ************************************************************************ */
  /**
   * Apply Zalesaks flux limiter as described in Kuzmin 2009 pp.6-7
   * to compute flux correction factors. Used with FEM-FCT for time-
   * dependent CDR problems.
   *
   * @param[out] alphas The flux correction factors, row by row (sorted like
   * the entries array of the mass_matrix!).
   * MPI: alphas will be ROWWISE LEVEL-0-CONSISTENT IN MPI on output
   *
   * @param[in] mass_matrix The complete mass matrix.
   * MPI: mass_matrix MUST BE ROWWISE LEVEL-0-CONSISTENT IN MPI.
   *
   * @param[in] lumped_mass The row-wise lumped mass matrix as vector
   * of its diagonal values.
   * MPI: lumped_mass MUST BE ROWWISE LEVEL-0-CONSISTENT IN MPI.
   *
   * @param[in] raw_fluxes The raw fluxes, sorted like alphas. Must have been
   * multiplied with current timesteplength!
   * MPI: raw_fluxes MUST BE ROWWISE LEVEL-0-CONSISTENT IN MPI
   *
   * @param[in] sol_approx The approximate solution used to calculate
   * the flux correction.
   * MPI: sol_approx MUST BE LEVEL-2-CONSISTENT IN MPI.
   *
   * @param[in] neum_to_diri See FEM-FCT, same story.
   * @note I'm pretty sure we can get rid of the whole "neum_to_diri"
   * thing when passing an FEMatrix here which got assembled with Neumann
   * treatment of dirichlet rows.
   */
  void ZalesaksFluxLimiter(
    std::vector<double>& alphas,
    const FEMatrix& mass_matrix,
    const std::vector<double>& lumped_mass,
    const std::vector<double>& raw_fluxes,
    const std::vector<double>& sol_approx,
    const std::vector<int>& neum_to_diri
    )
  {

#ifdef _MPI
    const TParFECommunicator3D& comm = mass_matrix.GetFESpace3D()->get_communicator();
    const int* masters = comm.GetMaster();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    //check if all dimensions fit
    if (!mass_matrix.is_square() )
    {
      ErrThrow("Matrix must be square!");
    }

    int nDofs = mass_matrix.GetN_Rows();
    int nEntries = mass_matrix.GetN_Entries();

    if((int) lumped_mass.size() != nDofs)
    {
      ErrThrow("lumped_mass has incorrect size ", lumped_mass.size(),
        " != ", nDofs);
    }

    if((int) sol_approx.size() != nDofs)
    {
      ErrThrow("sol_approx has incorrect size ", sol_approx.size(),
        " != ", nDofs);
    }

    if((int) alphas.size() != nEntries)
    {
      ErrThrow("alphas has incorrect size ", alphas.size(), " != ", nEntries);
    }

    if((int) raw_fluxes.size() != nEntries)
    {
      ErrThrow("raw_fluxes has incorrect size ", raw_fluxes.size(),
        " != ", nEntries);
    }
    //end checking input

    // get pointers to columns, rows and entries of mass_matrix
    const int* ColInd_M = mass_matrix.GetKCol();
    const int* RowPtr_M = mass_matrix.GetRowPtr();

    // auxiliary space
    std::vector<double> P_plus (nDofs, 0.0);
    std::vector<double> P_minus (nDofs, 0.0);
    std::vector<double> Q_plus (nDofs, 0.0);
    std::vector<double> Q_minus (nDofs, 0.0);

    std::vector<double> R_plus (nDofs, 0.0);
    std::vector<double> R_minus (nDofs, 0.0);

    // Compute Ps and Qs (negative/positive diffusive/antidiffusive fluxes),
    // where the Qs are "distances to local extrema of the auxiliary solution"
    // MPI: Needs raw_fluxes matrix rowwise level-0-consistent and
    //      sol_approx in level-2-consistency
    for(int i=0;i<nDofs;i++)
    {
#ifdef _MPI
      //skip the non-master rows in MPI case
      if(masters[i] != rank)
        continue;
#endif

      // i-th row of mass matrix
      int j0 = RowPtr_M[i];
      int j1 = RowPtr_M[i+1];
      for(int j=j0;j<j1;j++)
      {
        int index = ColInd_M[j];

        if(i != index)                            //diagonal entries are skipped
        {
          if (raw_fluxes[j] > 0.0)
            P_plus[i] += raw_fluxes[j];           //P_plus was initialized with 0
          if (raw_fluxes[j] <= 0.0)
            P_minus[i] += raw_fluxes[j];          //P_minus was initialized with 0

          double help = sol_approx[index]-sol_approx[i];

          if (help > Q_plus[i])                   //Q_plus was initialized with 0
            Q_plus[i]= help;

          if (help < Q_minus[i])                  //Q_minus was initialized with 0
            Q_minus[i]= help;
        }
      }                                           // end loop j
    }

    // Compute R_plus and R_minus (nodal correction factors)
    // MPI: Needs lumped_mass, Q_plus, Q_minus, P_plus, P_minus in Level 0 consistency
    for(int i=0;i<nDofs;i++)
    {

#ifdef _MPI
      //skip the non-master rows in MPI case
      if(masters[i] != rank)
        continue;
#endif

      //determine R_plus
      if(fabs(P_plus[i]) == 0.0)
        R_plus[i] = 1.0;
      else
      {
        //OutPut(Q_plus[i] << " "  << P_plus[i] << " " << lump_mass[i] << " ");
        R_plus[i] = std::min(1.0 , (lumped_mass[i] * Q_plus[i])/P_plus[i] );
      }

      //determine R_minus
      if(fabs(P_minus[i]) == 0.0)
        R_minus[i] = 1.0;
      else
      {
        //OutPut(Q_minus[i] << " "  << P_minus[i] << " " << lump_mass[i] << " ");
        R_minus[i] = std::min(1.0, (lumped_mass[i] * Q_minus[i])/P_minus[i]);
      }

    }

    // treat Dirichlet nodes (e.g. Kuzmin & Moeller 2005 S. 27)
    // TODO actually there it is about inlet and outlet!
    for (int j=0;j< (int) neum_to_diri.size();j++)
    {
      int i = neum_to_diri[j];
      R_plus[i] = 1;
      R_minus[i] = 1;
    }

#ifdef _MPI
    // put R_plus and R_minus into level 2 consistency
    comm.consistency_update(R_plus.data(), 2);
    comm.consistency_update(R_minus.data(), 2);
#endif

    // Compute alphas (final correction factors)
    // MPI: needs R_plus and R_minus in level 2 consistency
    for( int i=0 ; i<nDofs ; i++ )
    {

#ifdef _MPI
      //skip the non-master rows in MPI case
      if(masters[i] != rank)
        continue;
#endif

      // i-th row of sqmatrix
      int j0 = RowPtr_M[i];
      int j1 = RowPtr_M[i+1];

      for(int j=j0;j<j1;j++)
      {
        int index = ColInd_M[j];
        if(raw_fluxes[j] > 0.0)
        {
          alphas[j] = std::min(R_plus[i], R_minus[index]);
        }
        if(raw_fluxes[j]<=0.0)
        {
          //initialisation
          alphas[j] = std::min(R_minus[i], R_plus[index]);
        }
        // clipping, see Kuzmin (2009), end of Section 5
        //if ((fabs(Q_plus[i])< eps)&&(fabs(Q_minus[i])< eps))
        //  alpha[j] = 0;
      }                                           //end loop j
    }

    //MPI: At output, correction factors alpha_ij are rowwise level-0-consistent
  }

/** ************************************************************************ */  
  /** Computation of the Jacobian times the fluxes for the Zalesak limiter
   *            used in Newton's method
   * @param[in] A: The system matrix A
   * @param[in] F: The Matrix where entries f_ij=d_ij-(uj_-u_i)
   * @param[in] P_plus, Q_plus, P_minus, Q_minus
   *            R_plus, R_minus: The fluxes used for computation of the 
   *            limiters alpha_ij
   * @param[in] index_i: Denotes the ith row of matrix DF for which summation 
   *                     is required
   * @param[in] index_j: Column number for entry a_ij
   * @param[in] entries_pointer_ij: Denotes the pointer to the array Entries
   * @param[in] D: The artificial diffusion matrix
   *
   * @param[out] The summation required in the entry of the Jacobian matrix.
   */
  double Compute_Jacobian_times_flux_Kuzmin(const FEMatrix& A,
    const double * F ,
    const double * P_plus,
    const double * Q_plus,
    const double * Q_minus,
    const double * P_minus,
    const double * R_minus,
    const double * R_plus,
    const int index_i,
    const int index_j,
    const int entries_pointer_ij,
    const FEMatrix& D)
  {
    double sum=0, epsinv;
    //DQ_minus, DQ_plus, DP_minus, DP_plus: 
    //             Denotes the Derivative of the raw fluxes.
    double DQ_plus, DQ_minus,DP_plus,DP_minus, sumD=0;
    const double eps = 1e-14;
    const int* ColInd = A.GetKCol();
    const int* RowPtr = A.GetRowPtr();
    const double* Entries = A.GetEntries();
    int row_i_start, row_i_end, temp, index_a_ji, index_k, index_a_ki;
    int row_k_start, row_k_end, index_a_kj, index_a_jk, a_kj;
    
    const double* afc_matrix_D_entries=D.GetEntries();

    epsinv = 1.0/eps;
    // computation of the index of entry a_ji
    for(int mm=RowPtr[index_j];mm<RowPtr[index_j+1];mm++)
    {
      if(ColInd[mm]==index_i)
      {
        index_a_ji=mm;
        break;
      }
    }

    row_i_start=RowPtr[index_i];
    row_i_end=RowPtr[index_i+1];
    // loop over all dofs that are connected with the matrix entry a_{ij}
    for(int k=row_i_start; k<row_i_end;k++)
    {
      index_k=ColInd[k];
      row_k_start=RowPtr[index_k];
      row_k_end=RowPtr[index_k+1];
      index_a_kj=-1;
      //to find the index of the entry a_ki
      for(int jj=row_k_start;jj<row_k_end;jj++)
      {
        // index of a_kj
        if(ColInd[jj]==index_j)
          index_a_kj=jj;
        // index of transposed entry a_ki
        if(ColInd[jj]==index_i)
          index_a_ki=jj;
      }                                           //end of loop jj
      //case when a_ik>=a_ki
      if (Entries[k]>=Entries[index_a_ki])
      {
        // cases without contribution to the sum
        if(F[k]>0 && fabs(R_plus[index_i]-1)<eps)
          continue;
        if(F[k]<0 && fabs(R_minus[index_i]-1)<eps)
          continue;
        if(fabs(F[k])<eps)
          continue;
	if ((index_i!=index_j)&&fabs(afc_matrix_D_entries[entries_pointer_ij])<1e-14)
	  continue;
        // *** first case with contribution ***
        if(F[k]>0)
        {
          //Calculations for finding derivative of Q_plus
          if (index_i!=index_j)
          {
            if(F[entries_pointer_ij]>=0)
              DQ_plus=0;
            else
              DQ_plus=-afc_matrix_D_entries[entries_pointer_ij];
          }
          else
          {
            sumD=0;
            for(int m=row_i_start;m<row_i_end;m++)
            {
              if(F[m]<0)
                sumD+=afc_matrix_D_entries[m];
            }
            DQ_plus=sumD;
          }

          // calculations for finding derivative of P_plus
          if (index_i!=index_j)
          {
            if(F[entries_pointer_ij]<=0)
            {
              DP_plus=0;
            }
            else
            {
              if (Entries[entries_pointer_ij]>=Entries[index_a_ji])
                DP_plus=afc_matrix_D_entries[entries_pointer_ij];
              else
                DP_plus=0;
            }
          }
          else
          {
            sumD=0;
            for(int m=row_i_start;m<row_i_end;m++)
            {
              int transposed_entry=ColInd[m];
              for(int ll=RowPtr[transposed_entry];ll<RowPtr[transposed_entry+1];
                  ll++)
              {
                if(ColInd[ll]==index_i)
                {
                  temp=ll;
                  break;
                }
              }
              if(F[m]>0 && Entries[m]>=Entries[temp])
                sumD+=afc_matrix_D_entries[m];
            }                                     //End of loop m
            DP_plus=-sumD;
          }
          double val =  DQ_plus/P_plus[index_i]-Q_plus[index_i]*
                                DP_plus/(P_plus[index_i]*P_plus[index_i]);
          if (fabs(val) < epsinv)
            sum += val*F[k];
          continue;
        }
        // *** second case with contribution ***
        if(F[k]<0)
        {
          // calculations for finding derivative of Q_minus
          if (index_i!=index_j)
          {
            if(F[entries_pointer_ij]<=0)
              DQ_minus=0;
            else
              DQ_minus=-afc_matrix_D_entries[entries_pointer_ij];
          }
          else
          {
            sumD=0;
            for(int m=row_i_start;m<row_i_end;m++)
            {
              if(F[m]>0)
                sumD+=afc_matrix_D_entries[m];
            }
            DQ_minus=sumD;
          }
          //Calculations for finding derivative of P_minus
          if (index_i!=index_j)
          {
            if(F[entries_pointer_ij]>=0)
              DP_minus=0;
            else
            {
              if(Entries[entries_pointer_ij]>=Entries[index_a_ji])
                DP_minus=afc_matrix_D_entries[entries_pointer_ij];
              else
                DP_minus=0;
            }
          }
          else
          {
            sumD=0;
            for(int m=row_i_start;m<row_i_end;m++)
            {
              int transposed_entry=ColInd[m];
              for(int ll=RowPtr[transposed_entry];ll<RowPtr[transposed_entry+1];
                  ll++)
              {
                if(ColInd[ll]==index_i)
                {
                  temp=ll;
                  break;
                }
              }
              if(F[m]<0 && Entries[m]>=Entries[temp])
                sumD+=afc_matrix_D_entries[m];
            }                                     //End of loop m
            DP_minus=-sumD;
          }
          double val = DQ_minus/P_minus[index_i]-Q_minus[index_i]*
                          DP_minus/(P_minus[index_i]*P_minus[index_i]);
          if (fabs(val) < epsinv)
            sum += val*F[k];
          continue;
        }
        ErrThrow("Compute_Jacobian_times_flux_Kuzmin: index pair without case");
      }                                        //end of if case when a_ik>=a_ki
      //case when a_ki>a_ik
      else
      {
        // cases without contribution to the sum
        if(F[k]>0 && fabs(R_minus[index_k]-1)<eps)
          continue;
        if(F[k]<0 && fabs(R_plus[index_k]-1)<eps)
          continue;
        if(fabs(F[k])<eps)
          continue;
        double f_kj, d_kj;
        if(index_a_kj != -1)
        {
          f_kj=F[index_a_kj];
          d_kj=afc_matrix_D_entries[index_a_kj];
          a_kj=Entries[index_a_kj];
        }
        else
        {
          // next k
          continue;
        }     
        if ((index_k != index_j)&& (fabs(d_kj) < 1e-14))
	  continue;
        // computation of the index of entry a_jk
        for(int mm=RowPtr[index_j];mm<RowPtr[index_j+1];mm++)
        {
          if(ColInd[mm]==index_k)
          {
            index_a_jk=mm;
            break;
          }
        }
        // *** first case with contribution ***
        if(F[k]>0)
        {
          //Calculations for finding derivative of Q_minus
          if (index_k!=index_j)
          {
            if(f_kj<=0)
              DQ_minus=0;
            else
              DQ_minus=-d_kj;
          }
          else
          {
            sumD=0;
            for(int m=row_k_start;m<row_k_end;m++)
            {
              if(F[m]>0)
                sumD+=afc_matrix_D_entries[m];
            }
            DQ_minus=sumD;
          }
          // calculations for finding derivative of P_minus
          if (index_k!=index_j)
          {
            if(f_kj>=0)
            {
              DP_minus=0;
            }
            else
            {
              if (a_kj>=Entries[index_a_jk])
                DP_minus=d_kj;
              else
                DP_minus=0;
            }
          }
          else
          {
            sumD=0;
            for(int m=row_k_start;m<row_k_end;m++)
            {
              int transposed_entry=ColInd[m];
              for(int ll=RowPtr[transposed_entry];ll<RowPtr[transposed_entry+1];
                   ll++)
              {
                if(ColInd[ll]==index_k)
                {
                  temp=ll;
                  break;
                }
              }
              if(F[m]<0 && Entries[m]>=Entries[temp])
                sumD+=afc_matrix_D_entries[m];
            }                                     //End of loop m
            DP_minus=-sumD;
          }
          double val =  DQ_minus/P_minus[index_k]-Q_minus[index_k]*
                                 DP_minus/(P_minus[index_k]*P_minus[index_k]);
          if (fabs(val) < epsinv)
            sum += val*F[k];
          continue;
        }
        // *** second case with contribution ***
        if(F[k]<0)
        {
          // calculations for finding derivative of Q_plus
          if (index_k!=index_j)
          {
            if(f_kj>=0)
              DQ_plus=0;
            else
              DQ_plus=-d_kj;
          }
          else
          {
            sumD=0;
            for(int m=row_k_start;m<row_k_end;m++)
            {
              if(F[m]<0)
                sumD+=afc_matrix_D_entries[m];
            }
            DQ_plus=sumD;
            //Output::print<2>(DQ_plus);
          }
          //Calculations for finding derivative of P_plus
          if (index_k!=index_j)
          {
            if(f_kj<=0)
              DP_plus=0;
            else
            {
              if(a_kj>=Entries[index_a_jk])
                DP_plus=d_kj;
              else
                DP_plus=0;
            }
          }
          else
          {
            sumD=0;
            for(int m=row_k_start;m<row_k_end;m++)
            {
              int transposed_entry=ColInd[m];
              for(int ll=RowPtr[transposed_entry];ll<RowPtr[transposed_entry+1];
                   ll++)
              {
                if(ColInd[ll]==index_k)
                {
                  temp=ll;
                  break;
                }
              }
              if(F[m]>0 && Entries[m]>=Entries[temp])
                sumD+=afc_matrix_D_entries[m];
            }                                     //End of loop m
            DP_plus=-sumD;
          }
          double val = DQ_plus/P_plus[index_k]-Q_plus[index_k]*
                                   DP_plus/(P_plus[index_k]*P_plus[index_k]);
          if (fabs(val) < epsinv)
            sum += val*F[k];
          continue;
        }
        ErrThrow("Compute_Jacobian_times_flux_Kuzmin: index pair without case");
      }                                      //end of else case when a_ki>a_ik
    }
    return sum;
  }
  
/** ************************************************************************ */  

  /* regularized maximum function given as 
   * max{x,y}=(x+y+sqrt((x+y)^2+sigma))/2*/
  double smooth_max(const double x, const double y, const double sigma)
  {
    double value;
    value=(x+y+sqrt((x-y)*(x-y)+sigma))*0.5;
    return value;
  }

  /* regularixed minimum function given as 
   * min{x,y}=(x+y-sqrt((x-y)^2+sigma))/2*/
  double smooth_min(const double x, const double y, const double sigma)
  {
    double value;
    value=(x+y-sqrt((x-y)*(x-y)+sigma))*0.5;
    return value;
  }

  /* derivative of regularized minimum function*/
  double der_smooth_min(const double x, const double y, const double sigma)
  {
    double value;
    value=1+(x-y)/sqrt((x-y)*(x-y)+sigma);
    return value;
  }

/** ************************************************************************ */  
  /** Computation of the Jacobian times the fluxes for the regularized Kuzmin 
   *            limiter used in Regularized Newton's method
   * @param[in] A: The system matrix A
   * @param[in] F: The Matrix where entries f_ij=d_ij-(uj_-u_i)
   * @param[in] index_i: Denotes the ith row of matrix DF for which summation 
   *                     is required
   * @param[in] index_j: Column number for entry a_ij
   * @param[in] entries_pointer_ij: Denotes the pointer to the array Entries
   * @param[in] afc_matrix_D_entries: The artificial diffusion matrix
   *
   * @param[out] The summation required in the entry of the Jacobian matrix.
   */
  double Compute_Jacobian_times_flux_Kuzmin_Regularized(const FEMatrix& A,
    const double * F ,
    const int index_i,
    const int index_j,
    const int entries_pointer_ij,
    const FEMatrix& D, 
    const double sigma, 
    const double omega_matrix, const double omega_derivative)
  {
    //Regularized P_plus, P_minus, Q_plus, Q_minus
    double P_plus_i_reg=0.0;
    double P_minus_i_reg=0.0;
    double Q_plus_i_reg=0.0;
    double Q_minus_i_reg=0.0;
    double P_plus_k_reg, P_minus_k_reg, Q_plus_k_reg, Q_minus_k_reg;
    double a_kj, a_jk = 0, a_ji = 0;
    double sum=0;
    
    //DQ_minus, DQ_plus, DP_minus, DP_plus: Denotes the Derivative of the raw fluxes.
    double DQ_plus_i, DQ_minus_i, DP_plus_i,DP_minus_i, DF_ik, 
           DR_plus_i, DR_minus_i;
    double DQ_plus_k, DQ_minus_k, DP_plus_k,DP_minus_k, 
           DR_plus_k, DR_minus_k;
    double R_plus_i, R_minus_i, R_plus_k, R_minus_k;

    const int* ColInd = A.GetKCol();
    const int* RowPtr = A.GetRowPtr();
    const double* Entries = A.GetEntries();
    int row_i_start, row_i_end, temp = 0, index_a_ji = 0;
    
    // indicies used in the loop
    int k, m, index_m, m2, m3, mm, row_k_start, row_k_end, index_k;
    int index_a_ki = -1, index_a_kj, index_a_jk = 0;
    
    const double* afc_matrix_D_entries=D.GetEntries();

    // i-th row of sqmatrix
    row_i_start = RowPtr[index_i];
    row_i_end = RowPtr[index_i+1];

    // computation of the index of entry a_ji
    for(int mm=RowPtr[index_j];mm<RowPtr[index_j+1];mm++)
    {
      if(ColInd[mm]==index_i)
      {
        index_a_ji=mm;
        a_ji=Entries[index_a_ji];
        break;
      }
    }

    // computation of regularized Q_plus_i, Q_minus_i, P_plus_i, P_minus_i
    for(m=row_i_start;m<row_i_end;m++)
    {
      Q_minus_i_reg -= smooth_max(0,F[m],sigma);
      Q_plus_i_reg -= smooth_min(0,F[m],sigma);
      index_m=ColInd[m];
      m2=RowPtr[index_m];
      m3=RowPtr[index_m+1];
      //finding index of transposed entry a_im used in P_plus_i and P_minus_i
      for(mm=m2;mm<m3;mm++)
      {
        if(ColInd[mm]==index_i)
          break;
      }                                           //end of loop mm
      if(Entries[mm]>Entries[m])
        continue;
      P_plus_i_reg += smooth_max(0,F[m],sigma);
      P_minus_i_reg += smooth_min(0,F[m],sigma);
    }                                             //end of loop m
    R_plus_i = smooth_min(1,Q_plus_i_reg/P_plus_i_reg,sigma)/2.0;
    R_minus_i = smooth_min(1,Q_minus_i_reg/P_minus_i_reg,sigma)/2.0;

    // computation of derivatives in node i
    if (index_i!=index_j)
    {
      double d_ij = afc_matrix_D_entries[entries_pointer_ij];
      double f_ij = F[entries_pointer_ij];
      double f_ij_f_ij2_sigma = f_ij/sqrt(f_ij * f_ij + sigma);

      DQ_plus_i = -(1-f_ij_f_ij2_sigma)*d_ij/2.0;
      DQ_minus_i = -(1+f_ij_f_ij2_sigma)*d_ij/2.0;
      if (Entries[entries_pointer_ij]>=a_ji)
      {
        DP_plus_i = (1+f_ij_f_ij2_sigma)*d_ij/2.0;
        DP_minus_i = (1-f_ij_f_ij2_sigma)*d_ij/2.0;
      }
      else
      {
        DP_plus_i = DP_minus_i = 0;
      }
    }
    else                                          // index_i = index_j
    {
      // DQ_plus, DQ_minus
      DQ_plus_i=DQ_minus_i=0;
      for(int m=row_i_start;m<row_i_end;m++)
      {
        //the case k!=i
        if(ColInd[m]==index_i)
          continue;
        DQ_plus_i += (1-F[m]/sqrt(F[m]*F[m]+sigma))*afc_matrix_D_entries[m];
        DQ_minus_i += (1+F[m]/sqrt(F[m]*F[m]+sigma))*afc_matrix_D_entries[m];
      }
      DQ_plus_i /= 2.0;
      DQ_minus_i /= 2.0;
      // DP_plus, DP_minus
      DP_plus_i=DP_minus_i=0;
      for(int m=row_i_start;m<row_i_end;m++)
      {
        //the case m!=i
        if(ColInd[m]==index_i)
          continue;
        int transposed_entry=ColInd[m];
        for(int ll=RowPtr[transposed_entry];ll<RowPtr[transposed_entry+1];ll++)
        {
          if(ColInd[ll]==index_i)
          {
            temp=ll;
            break;
          }
        }
        if(Entries[m]>=Entries[temp])
        {
          DP_plus_i -= (1+F[m]/sqrt(F[m]*F[m]+sigma)) * afc_matrix_D_entries[m];
          DP_minus_i -= (1-F[m]/sqrt(F[m]*F[m]+sigma)) * afc_matrix_D_entries[m];
        }
      }                                           //End of loop m
      DP_plus_i /= 2.0;
      DP_minus_i /= 2.0;
    }
    DR_plus_i = ((DQ_plus_i*P_plus_i_reg-Q_plus_i_reg*DP_plus_i)/(P_plus_i_reg*P_plus_i_reg))
      *der_smooth_min(1,Q_plus_i_reg/(P_plus_i_reg),sigma)/2.0;
    DR_minus_i = ((DQ_minus_i*P_minus_i_reg-Q_minus_i_reg*DP_minus_i)/(P_minus_i_reg*P_minus_i_reg))
      *der_smooth_min(1,Q_minus_i_reg/(P_minus_i_reg),sigma)/2.0;

    //loop to find sum over beta_ik
    for(k=row_i_start;k<row_i_end;k++)
    {
      index_k=ColInd[k];
      row_k_start=RowPtr[index_k];
      row_k_end=RowPtr[index_k+1];
      index_a_kj=-1;
      //to find the index of the entry a_ki
      for(int jj=row_k_start;jj<row_k_end;jj++)
      {
        // index of a_kj
        if(ColInd[jj]==index_j)
          index_a_kj=jj;
        // index of transposed entry a_ki
        if(ColInd[jj]==index_i)
          index_a_ki=jj;
      }                                           //end of loop jj

      //case when a_ik>=a_ki
      if (Entries[k]>=Entries[index_a_ki])
      {

        // derivative of F_ik
        if (index_i==index_j)
        {
          if(index_k==index_j)
            DF_ik=0;
          else
            DF_ik=-afc_matrix_D_entries[k];
        }
        else
        {
          if(index_k==index_j)
            DF_ik=afc_matrix_D_entries[entries_pointer_ij];
          else
            DF_ik=0;
        }
        sum += omega_derivative * DR_plus_i*smooth_max(0,F[k],sigma) 
              + omega_matrix * R_plus_i*der_smooth_min(F[k],0,sigma)*DF_ik
              + omega_derivative * DR_minus_i*smooth_min(0,F[k],sigma) 
              + omega_matrix * R_minus_i*der_smooth_min(0,F[k],sigma)*DF_ik;
      }
      else                                        //case when a_ki>a_ik
      {
        double f_kj, d_kj;
        if(index_a_kj != -1)
        {
          f_kj=F[index_a_kj];
          d_kj=afc_matrix_D_entries[index_a_kj];
          a_kj=Entries[index_a_kj];
        }
        else
        {
          // next k
          continue;
        }
        // computation of the index of entry a_jk
        for(int mm=RowPtr[index_j];mm<RowPtr[index_j+1];mm++)
        {
          if(ColInd[mm]==index_k)
          {
            index_a_jk = mm;
            a_jk = Entries[index_a_jk];
            break;
          }
        }

        // computation of regularized Q_plus_k, Q_minus_k, P_plus_k, P_minus_k
        P_plus_k_reg = P_minus_k_reg = Q_plus_k_reg =  Q_minus_k_reg = 0.0;
        for(m=row_k_start;m<row_k_end;m++)
        {
          Q_minus_k_reg -= smooth_max(0,F[m],sigma);
          Q_plus_k_reg -= smooth_min(0,F[m],sigma);
          index_m=ColInd[m];
          m2=RowPtr[index_m];
          m3=RowPtr[index_m+1];
          //finding index of transposed entry a_km used in P_plus_k and P_minus_k
          for(mm=m2;mm<m3;mm++)
          {
            if(ColInd[mm]==index_k)
              break;
          }                                       //end of loop mm
          if(Entries[mm]>Entries[m])
            continue;
          P_plus_k_reg += smooth_max(0,F[m],sigma);
          P_minus_k_reg += smooth_min(0,F[m],sigma);
        }
        R_plus_k = smooth_min(1,Q_plus_k_reg/P_plus_k_reg,sigma)/2.0;
        R_minus_k = smooth_min(1,Q_minus_k_reg/P_minus_k_reg,sigma)/2.0;
        
        // computation of derivatives in node k
        if (index_k!=index_j)
        {
          //double d_kj = afc_matrix_D_entries[entries_pointer_ij];
          //double f_kj = F[entries_pointer_ij];
     
          double f_kj_f_kj2_sigma = f_kj/sqrt(f_kj * f_kj + sigma);
          DQ_plus_k = -(1-f_kj_f_kj2_sigma)*d_kj/2.0;
          DQ_minus_k = -(1+f_kj_f_kj2_sigma)*d_kj/2.0;
          
          if (a_kj >= a_jk)
          {
            DP_plus_k = (1+f_kj_f_kj2_sigma)*d_kj/2.0;
            DP_minus_k = (1-f_kj_f_kj2_sigma)*d_kj/2.0;
          }
          else
          {
            DP_plus_k = DP_minus_k = 0;
          }
        }
        else                                          // index_k = index_j
        {
          // DQ_plus, DQ_minus
          DQ_plus_k=DQ_minus_k=0;
          for(int m=row_k_start;m<row_k_end;m++)
          {
            //the case k!=m
            if (ColInd[m]==index_k)
              continue;
            DQ_plus_k += (1-F[m]/sqrt(F[m]*F[m]+sigma))*afc_matrix_D_entries[m];
            DQ_minus_k += (1+F[m]/sqrt(F[m]*F[m]+sigma))*afc_matrix_D_entries[m];
          }
          DQ_plus_k /= 2.0;
          DQ_minus_k /= 2.0;
          // DP_plus, DP_minus
          DP_plus_k=DP_minus_k=0;
          for(int m=row_k_start;m<row_k_end;m++)
          {
            //the case m!=k
            if(ColInd[m]==index_k)
              continue;
            int transposed_entry=ColInd[m];
            for(int ll=RowPtr[transposed_entry];ll<RowPtr[transposed_entry+1];ll++)
            {
              if(ColInd[ll]==index_k)
              {
                temp=ll;
                break;
              }
            }
            if(Entries[m]>=Entries[temp])
            {
              DP_plus_k -= (1+F[m]/sqrt(F[m]*F[m]+sigma))
                          *afc_matrix_D_entries[m];
              DP_minus_k -= (1-F[m]/sqrt(F[m]*F[m]+sigma))
                          *afc_matrix_D_entries[m];
            }
          }                                           //End of loop m
          DP_plus_k /= 2.0;
          DP_minus_k /= 2.0;
        }
        DR_plus_k = ((DQ_plus_k*P_plus_k_reg-Q_plus_k_reg*DP_plus_k)/(P_plus_k_reg*P_plus_k_reg))
          *der_smooth_min(1,Q_plus_k_reg/P_plus_k_reg,sigma)/2.0;
        DR_minus_k = ((DQ_minus_k*P_minus_k_reg-Q_minus_k_reg*DP_minus_k)/(P_minus_k_reg*P_minus_k_reg))
          *der_smooth_min(1,Q_minus_k_reg/P_minus_k_reg,sigma)/2.0;

        // derivative of F_ik
        if (index_k==index_j)
        {
          if(index_i==index_j)
            DF_ik=0;
          else
            DF_ik=-afc_matrix_D_entries[k];
        }
        else
        {
          if(index_i==index_j)
            DF_ik= afc_matrix_D_entries[k];
          else
            DF_ik=0;
        }
        sum += -omega_derivative * DR_plus_k*smooth_max(0,-F[k],sigma) 
               - omega_matrix * R_plus_k*der_smooth_min(-F[k],0,sigma)*DF_ik
               -omega_derivative * DR_minus_k*smooth_min(0,-F[k],sigma) 
               - omega_matrix * R_minus_k*der_smooth_min(0,-F[k],sigma)*DF_ik;
      }
    }
    return sum;
  }

/** ************************************************************************ */  
  
  /** Computation of the Jacobian times the fluxes for the BJK17 limiter
   *            used in regularixed Newton's method
   * @param[in] A: The system matrix A
   * @param[in] F: The Matrix where entries f_ij=d_ij-(uj_-u_i)
   * @param[in] P_plus, Q_plus, P_minus, Q_minus, R_plus, R_minus: 
   *            The fluxes used for computation of the limiters alpha_ij
   * @param[in] umax: Vector where umax_i=max{u_ij, 1<=j<=Ndofs}
   * @param[in] umin: Vector where umin_i=min{u_ij, 1<=j<=Ndofs}
   * @param[in] q: Vector where q_i=gamma_i\sum d_ij
   * @param[in] sol_col_j: entry of sol[col_j]
   * @param[in] row_i: Denotes the ith row of matrix DF for which summation 
   *                   is required
   * @param[in] col_j: Column number for entry a_ij
   * @param[in] entries_pointer: Denotes the pointer to the array Entries
   * @param[in] D: The artificial diffusion matrix
   *
   * @param[out] The summation required in the entry of the Jacobian matrix.
   */

  double Compute_Jacobian_times_flux_BJK17(const FEMatrix& A,
    const double * F ,
    const double * P_plus,
    const double * P_minus,
    const double * Q_plus,
    const double * Q_minus,
    const double * R_plus,
    const double * R_minus,
    const std::vector<double>& umax,
    const std::vector<double>& umin,
    const std::vector<double>& q,
    const double sol_col_j,
    const int row_i,
    const int col_j,
    const int entries_pointer,
    const FEMatrix& D)

  {
    double sum=0, epsinv;
    double DQ_plus, DQ_minus,DP_plus,DP_minus, sumD=0;
    const double eps = 1e-14;
    const int* ColInd = A.GetKCol();
    const int* RowPtr = A.GetRowPtr();
    int k0, k1, m0, m1;

    epsinv = 1.0/eps;

    k0=RowPtr[row_i];
    k1=RowPtr[row_i+1];
    m0=RowPtr[col_j];
    m1=RowPtr[col_j+1];
    const double* afc_matrix_D_entries=D.GetEntries();

    for(int k=k0;k<k1;k++)
    {
      // situations that do not lead to contributions
      if(F[k]>0 && fabs(R_plus[row_i]-1.0) < eps)
        continue;
      if(F[k]<0 && fabs(R_minus[row_i]-1.0) < eps)
        continue;
      if(fabs(F[k]) < eps)
        continue;

      int col_k = ColInd[k];
      // *** first condition with a contribution ***
      // if(F[k]>0 && R_plus[row_i]<1)
      if(F[k]>0)
      {
        if(R_plus[row_i]<=R_minus[col_k])
        {
          // computation of derivatives of Q_plus and P_plus
          // NOTE: the index corresponding to umax need not to be uniquely 
          //       defined
          DQ_plus = DP_plus = 0.0;
          if(row_i!=col_j)
          {
            if(fabs(umax[row_i]-sol_col_j) < eps)
              DQ_plus=-q[row_i];
            if(F[entries_pointer]>0)
              DP_plus=afc_matrix_D_entries[entries_pointer];
          }
          else
          {
            if(fabs(umax[row_i]-sol_col_j) >= eps)
              DQ_plus=q[row_i];
            sumD=0.0;
            for(int m=k0;m<k1;m++)
            {
              if(F[m]>0)
                sumD+=afc_matrix_D_entries[m];
            }
            DP_plus=-sumD;
          }
          if (P_plus[row_i] > 0.0)
          {
            double val = DQ_plus/P_plus[row_i]
                       - Q_plus[row_i]*DP_plus/(P_plus[row_i]*P_plus[row_i]);
            if (fabs(val) < epsinv)
              sum+=val*F[k];
          }
          continue;
        }
        else
        {
          // computation of derivatives of Q_minus and P_minus
          DQ_minus = DP_minus = 0.0;
          if(col_k!=col_j)
          {
            if(fabs(umin[col_k]-sol_col_j) < eps)
              DQ_minus=-q[col_k];
            if(F[entries_pointer]<0)
              DP_minus=afc_matrix_D_entries[entries_pointer];
          }
          // col_k = col_j
          else
          {
            if(fabs(umin[col_j]-sol_col_j) >= eps)
              DQ_minus=q[col_j];
            sumD=0;
            for(int m=m0;m<m1;m++)
            {
              if(F[m]<0)
                sumD+=afc_matrix_D_entries[m];
            }
            DP_minus=-sumD;
          }
          if (P_minus[col_k] < 0.0)
          {
            double val = DQ_minus/P_minus[col_k] 
                       - Q_minus[col_k]*DP_minus/(P_minus[col_k]*P_minus[col_k]);
            if (fabs(val) < epsinv)
              sum+=val*F[k];
          }
          continue;
        }
      }
      // *** second condition with a contribution ***
      //if(F[k]<0 && R_minus[row_i]<1)
      if (F[k] < 0)
      {
        if(R_minus[row_i]<=R_plus[col_k])
        {
          DQ_minus = DP_minus = 0.0;
          // computation of derivatives of Q_minus and P_minus
          if(row_i!=col_j)
          {
            if(fabs(umin[row_i]-sol_col_j) < eps)
              DQ_minus=-q[row_i];
            if(F[entries_pointer]<0)
              DP_minus=afc_matrix_D_entries[entries_pointer];
          }
          else
          {
            if(fabs(umin[row_i]-sol_col_j) >= eps)
              DQ_minus=q[row_i];
            sumD=0;
            for(int m=k0;m<k1;m++)
            {
              if(F[m]<0)
                sumD+=afc_matrix_D_entries[m];
            }
            DP_minus=-sumD;
          }
          if (P_minus[row_i] < 0.0)
          {
            double val = DQ_minus/P_minus[row_i]
                       - Q_minus[row_i]*DP_minus/(P_minus[row_i]*P_minus[row_i]);
            if (fabs(val) < epsinv)
              sum+=val*F[k];
          }
          continue;
        }
        // R_minus[row_i] > R_plus[ColInd[k]]
        else
        {
          DQ_plus = DP_plus = 0.0;
          // computation of derivatives of Q_plus and P_plus
          if(col_k!=col_j)
          {
            if(fabs(umax[col_k]-sol_col_j)<eps)
              DQ_plus=-q[col_k];
            if(F[entries_pointer]>0)
              DP_plus=afc_matrix_D_entries[entries_pointer];
          }
          // col_k = col_j
          else
          {
            if(fabs(umax[col_j]-sol_col_j) >= eps)
              DQ_plus=q[col_j];
            sumD=0;
            for(int m=m0;m<m1;m++)
            {
              if(F[m]>0)
                sumD+=afc_matrix_D_entries[m];
            }
            DP_plus=-sumD;
          }
          if (P_plus[col_k] > 0.0)
          {
            double val = DQ_plus/P_plus[col_k] 
                       - Q_plus[col_k]*DP_plus/(P_plus[col_k]*P_plus[col_k]);
            if (fabs(val) < epsinv)
              sum+=val*F[k];
          }
          continue;
        }
      }
      ErrThrow("Compute_Jacobian_times_flux_BJK17: index pair without case",
               k," ",F[k]);
    }
    return sum;
  }

/** ************************************************************************ */  
  
  /**
   * Compute the weights for the linearty preserving limiter from 
   * Barrenechea, John, Knobloch; M3AS 2017.
   *
   * @param[in] FEMatrix system matrix
   * @param[in] current solution vector
   * @param[out] gamma vector with the weights
   *
   * NOTE: The first step of [BJK17], Rem. 6.2 is not performed, i.e., no 
   *       vertex is shifted. For non-convex patches, the denominator of gamma_i
   *       might become too small, i.e. gamma_i too large. This might lead to a 
   *       smearing of layers but does not change the linearity preservation
   */

  void Compute_Parameter_For_Linearity_Preservation(FEMatrix& system_matrix,
                                                    const std::vector<double>& u,
                                                    std::vector<double>& gamma)
  {
#ifdef __2D__

    int i, j, index, N_Cells, N_Unknowns, N_V;
    int *global_numbers, *begin_index, *dof;
    double area, edge, x[3], y[3], dist, dist_max, dist_max2;
    FE2D CurrentElement;

    Output::print<4>("AFC: compute parameter for linearity preservation");
    // get fe space
    auto fespace = system_matrix.GetFESpace2D();
    // get collection
    auto Coll = fespace->GetCollection();
    // number of mesh cells
    N_Cells = Coll->GetN_Cells();
    // array with global numbers of d.o.f.
    global_numbers = fespace->GetGlobalNumbers();
    // array which points to the beginning of the global numbers in
    // global_numbers for each mesh cell
    begin_index = fespace->GetBeginIndex();
    // number of unknowns
    N_Unknowns = u.size();
    // allocate memory for the parameter
    std::vector<double> parameter(N_Unknowns);
    parameter.resize(2*N_Unknowns, 4711.0);

    // loop over the mesh cells
    for(i=0;i<N_Cells;i++)
    {
      auto cell = Coll->GetCell(i);
      N_V = cell->GetN_Vertices();
      // get pointer to local dofs
      dof = global_numbers+begin_index[i];
      // finite element on the mesh cell
      CurrentElement = fespace->GetFE2D(i, cell);
      // only for P_1
      if (CurrentElement != C_P1_2D_T_A)
      {
        ErrThrow("Compute_Parameter_For_Linearity_Preservation only implemented"
                  "for P_1 !!!");
      }
      area = cell->GetMeasure();
      // loop over the vertices
      for (j=0;j<N_V;j++)
      {
        // get coordinates
        cell->GetVertex(j)->GetCoords(x[j], y[j]);
      }
      for (j=0;j<N_V;j++)
      {
        index = dof[j];
        // check for obtuse angles with opposite edge
        // vertex j+1
        if ((x[(j+2)%N_V] - x[(j+1)%N_V])*(x[(j)] - x[(j+1)%N_V])
          + (y[(j+2)%N_V] - y[(j+1)%N_V])*(y[(j)] - y[(j+1)%N_V]) <=0)
        {
          dist = sqrt((x[j] - x[(j+1)%N_V])*(x[(j)] - x[(j+1)%N_V])
            + (y[j] - y[(j+1)%N_V])*(y[(j)] - y[(j+1)%N_V]));
          dist_max = sqrt((x[j] - x[(j+2)%N_V])*(x[(j)] - x[(j+2)%N_V])
            + (y[j] - y[(j+2)%N_V])*(y[(j)] - y[(j+2)%N_V]));
          // OutPut("obt1 ");
        }
        else
        {
          // check for obtuse angles with opposite edge
          // vertex j+2
          if ((x[(j+1)%N_V] - x[(j+2)%N_V])*(x[(j)] - x[(j+2)%N_V])
            + (y[(j+1)%N_V] - y[(j+2)%N_V])*(y[(j)] - y[(j+2)%N_V]) <=0)
          {
            dist = sqrt((x[j] - x[(j+2)%N_V])*(x[(j)] - x[(j+2)%N_V])
              + (y[j] - y[(j+2)%N_V])*(y[(j)] - y[(j+2)%N_V]));
            dist_max = sqrt((x[j] - x[(j+1)%N_V])*(x[(j)] - x[(j+1)%N_V])
              + (y[j] - y[(j+1)%N_V])*(y[(j)] - y[(j+1)%N_V]));
            // OutPut("obt2 ");
          }
          else
            // no obtuse angle
          {
            // length opposite egde
            edge = (x[(j+1)%N_V] - x[(j+2)%N_V]) *  (x[(j+1)%N_V] - x[(j+2)%N_V]);
            edge +=  (y[(j+1)%N_V] - y[(j+2)%N_V]) *  (y[(j+1)%N_V] - y[(j+2)%N_V]);
            edge = sqrt(edge);
            dist = 2.*area/edge;
            //OutPut(edge << " ");
            dist_max = sqrt((x[j] - x[(j+1)%N_V])*(x[(j)] - x[(j+1)%N_V])
              + (y[j] - y[(j+1)%N_V])*(y[(j)] - y[(j+1)%N_V]));
            dist_max2 = sqrt((x[j] - x[(j+2)%N_V])*(x[(j)] - x[(j+2)%N_V])
              + (y[j] - y[(j+2)%N_V])*(y[(j)] - y[(j+2)%N_V]));
            if (dist_max2 > dist_max)
              dist_max = dist_max2;
          }
        }
        // compute numerator
        if (dist_max > parameter[index])
          parameter[index] = dist_max;
        // compute denominator
        if (dist < parameter[index+N_Unknowns])
          parameter[index+N_Unknowns] = dist;
      }
    }

    // loop over the unknowns
    for (i=0;i<N_Unknowns;i++)
    {
      gamma[i] = parameter[i]/ parameter[i+N_Unknowns];
    }
#endif
#ifdef __3D__
    int i, j, index, N_Cells, N_Unknowns, N_V;
    int *global_numbers, *begin_index, *dof;
    //A,B,C,D: Coefficients of the plane oppostie the vertex (x[j],y[j],z[j])
    double dist, dist_max, dist_max1, dist_max2, dist_max3, x[4], y[4], z[4];
    double A, B, C, area;
    FE3D CurrentElement;    
    Output::print<4>("AFC: compute parameter for linearity preservation");
    auto fespace=system_matrix.GetFESpace3D();
    // get collection  
    auto Coll = fespace->GetCollection();
    // number of mesh cells
    N_Cells = Coll->GetN_Cells();
    // array with global numbers of d.o.f.
    global_numbers = fespace->GetGlobalNumbers();
    // array which points to the beginning of the global numbers in
    // global_numbers for each mesh cell
    begin_index = fespace->GetBeginIndex();  
    // number of unknowns
    N_Unknowns = u.size();
    // allocate memory for the parameter
    std::vector<double> parameter(N_Unknowns);
    parameter.resize(2*N_Unknowns, 4711.0);
	
    // loop over the mesh cells
    for(i=0;i<N_Cells;i++)
    {
      auto cell = Coll->GetCell(i);
      N_V = cell->GetN_Vertices();
      // get pointer to local dofs
      dof = global_numbers+begin_index[i];
      // finite element on the mesh cell
      CurrentElement = fespace->GetFE3D(i, cell);
      // only for P_1
      if (CurrentElement != C_P1_3D_T_A)
      {
	ErrThrow("Compute_Parameter_For_Linearity_Preservation only implemented"
	         "for P_1 !!!");
      }
      area=cell->GetMeasure();
      // loop over the vertices
      for (j=0;j<N_V;j++)
      {
	// get coordinates
	cell->GetVertex(j)->GetCoords(x[j], y[j],z[j]);  
      }    
      for (j=0;j<N_V;j++)
      {
	index = dof[j];
	dist_max1= std::sqrt((x[j] - x[(j+1)%N_V])*(x[(j)] - x[(j+1)%N_V])
	+ (y[j] - y[(j+1)%N_V])*(y[(j)] - y[(j+1)%N_V])
	+ (z[j] - z[(j+1)%N_V])*(z[(j)]-z[(j+1)%N_V]));
	dist_max2= std::sqrt((x[j] - x[(j+2)%N_V])*(x[(j)] - x[(j+2)%N_V])
	+ (y[j] - y[(j+2)%N_V])*(y[(j)] - y[(j+2)%N_V])
	+ (z[j] - z[(j+2)%N_V])*(z[(j)]-z[(j+2)%N_V]));
	dist_max3= std::sqrt((x[j] - x[(j+3)%N_V])*(x[(j)] - x[(j+3)%N_V])
	+ (y[j] - y[(j+3)%N_V])*(y[(j)] - y[(j+3)%N_V])
	+ (z[j] - z[(j+3)%N_V])*(z[(j)]-z[(j+3)%N_V]));
	dist_max=std::max(dist_max1,dist_max2); 
	dist_max=std::max(dist_max, dist_max3);
	A=(z[(j+3)%N_V]-z[(j+1)%N_V])*(y[(j+2)%N_V]-y[(j+1)%N_V])
   -(z[(j+2)%N_V]-z[(j+1)%N_V])*(y[(j+3)%N_V]-y[(j+1)%N_V]);
	B=(z[(j+2)%N_V]-z[(j+1)%N_V])*(x[(j+3)%N_V]-x[(j+1)%N_V])
   -(z[(j+3)%N_V]-z[(j+1)%N_V])*(x[(j+2)%N_V]-x[(j+1)%N_V]);
	C=(y[(j+3)%N_V]-y[(j+1)%N_V])*(x[(j+2)%N_V]-x[(j+1)%N_V])
   -(x[(j+3)%N_V]-x[(j+1)%N_V])*(y[(j+2)%N_V]-y[(j+1)%N_V]);
	
	// OTHER IMPLEMENTATION
	//Equation of a plane with three points (x[j+1],y[j+1],z[j+1]),
  // (x[j+2],y[j+2],z[j+2]) and (x[j+3],y[j+3],z[j+3])
	//D=-[A*x[(j+1)%N_V]+B*y[(j+1)%N_V]+C*z[(j+1)%N_V]];
	//Minimum distance between a point (x[j],y[j],z[j]) and a plane Ax+By+Cz+D=0
	//dist=abs(A*x[j]+B*y[j]+C*z[j]+D)/sqrt(A*A+B*B+C*C);
	
	//Idea from BJK17
	dist=(6.0*area)/std::sqrt(A*A+B*B+C*C);
	// compute numerator
	if (dist_max > parameter[index])
	  parameter[index] = dist_max;
	// compute denominator
	if (dist < parameter[index+N_Unknowns])
	  parameter[index+N_Unknowns] = dist;  
      }
    }
    for (i=0;i<N_Unknowns;i++)
    {
      gamma[i] =(parameter[i]/ parameter[i+N_Unknowns]);  
    }
#endif
  }
}

/** ************************************************************************ */

ParameterDatabase AlgebraicFluxCorrection::default_afc_database()
{
  ParameterDatabase db("default algebraic flux correction database");

  // Type of AFC to be applied.
  db.add("algebraic_flux_correction", "none", " Chose which type of afc to use.",
    {"none", "afc", "fem-fct-cn"}
  );

  db.add("afc_prelimiter", 0, "Choose an afc flux prelimiting scheme. Options "
    "are 0 (none), 1 (min-mod), 2 (grad-direction), 3 (both)", {0,1,2,3}
  );

  // type of the limiter to be applied for steady-state case
  db.add("afc_limiter", "kuzmin", "Choose an afc limiter. Options are"
    "kuzmin, BJK17", {"kuzmin", "BJK17"}
  );

  // type of the limiter to be applied for steady-state case
  db.add("afc_initial_iterate", "afc_zero", "Choose the initial iterate"
    "afc_zero, galerkin, supg, upwind", {"afc_zero", "galerkin", 
       "supg", "upwind"}
  );

  // iteration scheme for the afc methods
  db.add("afc_iteration_scheme", "fixed_point_matrix", 
         "Choose an iteration scheme for the afc methods. Options are"
    "fixed_point_rhs, fixed_point_matrix, newton, newton_regu",
    {"fixed_point_rhs", "fixed_point_matrix", "newton", 
      "newton_regu", "newton_no_damp"}
  );
  
  //maximum number of iterations in the non linear loop
  db.add("afc_nonlinloop_maxit", (size_t) 1 ,
         "Maximal number of iterations for the nonlinear loop in AFC."
         "Must be a value between 0 and 100000.", (size_t)  0, (size_t)100000);
  
  //tolerance for non linear loop in AFC scheme
  db.add("afc_nonlinloop_epsilon", 1e-10, 
         "Stopping criterion for the nonlinear loop in AFC."
         "Must be a value between 1e-20 an 1e20.", 1e-20, 1e20);
  

  //Constants related to Dynamic Damping from [JK08]
  db.add("afc_nonlinloop_damping_factor", 1.0, "A damping parameter"
          "for the nonlinear loop in AFC. Must be a value between 1" 
          "(no damping) and 0 (no update).", 0.0,1.0);

  db.add("afc_nonlinloop_damping_factor_max", 1.0, 
         "Maximal number for afc_nonlinloop_damping_factor."
         "Only changed internally", 0.0,1.0);

  db.add("afc_nonlinloop_damping_factor_min", 0.01, 
         "Minimal number for afc_nonlinloop_damping_factor."
         "Intended to stay constant", 0.0,1.0);
  
  db.add("afc_nonlinloop_damping_factor_max_global", 1.0, 
         "Maximal number for afc_nonlinloop_damping_factor."
         "for complete iteration", 0.0,1.0);

  db.add("afc_nonlinloop_damping_factor_increase", 1.1, 
         "Increase factor for afc_nonlinloop_damping_factor"
         "Intended to stay constant", 1.0,2.0);        // c_2 in [JK08]

  db.add("afc_nonlinloop_damping_factor_decrease", 0.5, 
         "Decrease factor for afc_nonlinloop_damping_factor"
         "Intended to stay constant", 0.0,1.0);

  db.add("afc_nonlinloop_damping_factor_max_increase", 1.001, 
         "Increase factor for afc_nonlinloop_damping_factor_max."
         "Intended to stay constant", 1.0,2.0);        // c_3 in [JK08]

  db.add("afc_nonlinloop_damping_factor_max_decrease", 0.9, 
         "Decrease factor for afc_nonlinloop_damping_factor_max."
         "Intended to stay constant", 0.0,1.0);        // c_4 in [JK08]

  db.add("afc_nonlinloop_damping_factor_min_tol", 1.001, 
         "Tolerance for afc_nonlinloop_damping_factor_min."
         "Intended to stay constant", 1.0,2.0);        // c_1 in [JK08]

  db.add("afc_nonlinloop_damping_factor_constant", "no", 
         "Whether or not afc_nonlinloop_damping_factor"
         "should be constant",                         //
  {
    "no", "yes"
  });
  
  //Anderson acceleration parameters
  db.add("afc_nonlinloop_anderson_acc", "yes", 
         "Whether or not to use Anderson acceleration", {"no", "yes"}
  );
  
  db.add("afc_nonlinloop_anderson_acc_vec", (size_t) 10, 
         "Number of vectors in Anderson acceleration",  
         (size_t)  1, (size_t)  1000);
  
  db.add("afc_nonlinloop_anderson_acc_start", (size_t) 0, 
         "Starting iterate for Anderson acceleration",  
         (size_t)  0, (size_t)  1000);
  
  db.add("afc_anderson_damping", "yes", "Anderson Damping from WN11 choice", 
         {"yes", "no"});

  //parameters related to fixed point matrix and Newton method
  db.add("afc_newton_regu_sigma", 1e-8, "Penalty for regularized Newton method"
    "is scaled with fourth power of mesh width", 0.0,1.0);

  db.add("afc_fixed_point_matrix_weight", 0.75, 
         "weight for implicit fixed point method"
         "fixed_point_matrix", 0.0,1.0);

  db.add("afc_fixed_point_derivative_weight", 0.1, 
         "weight for contribution of the derivative"
         "in Newton's method", 0.0,1.0);
  
  db.add("afc_fixed_point_derivative_weight_factor", 0.0, 
         "factor for weight for contribution of the derivative"
         "in Newton's method", 0.0,1.0);
   
  db.add("afc_damping_bound_newton", 0.02, 
         "bound for switching from fixed_point_rhs to Newton", 0.02, 1.0);
 
  db.add("afc_change_method_threshold", 1e-5, 
         "threshold for changing to Newton's method", 0.0,100.0);

  /*db.add("afc_nonlinloop_switch_to_newton_scheme", 1, 
   * "scheme for switching to formal Newton's method"
   *    "1 - standard, 10 - first fpr (for BAIL proceedings), 
   *    11 - first fpm (for BAIL proceedings)"
   *    "20 - first fpr with reswitch (for BAIL proceedings), 
   *    21 - first fpm with reswitch (for BAIL proceedings)"
   *    ,{1,10,11,20,21});*/
  
  //projection to admissible values. Idea from [BB17.CMAME] 
  db.add("afc_project_to_admissible_values", "no", 
         "Whether or not to project intermediate iterates to admissible values",
         {"no", "yes"});
   
  return db;
}


// ////////////////////////////////////////////////////////////////////////////
// Implementation of the methods in the namespace AlgebraicFluxCorrection    //
// ////////////////////////////////////////////////////////////////////////////

void AlgebraicFluxCorrection::steady_state_algorithm(
FEMatrix& system_matrix,
const std::vector<double>& sol,
std::vector<double>& rhs,
const std::vector<int>& neum_to_diri,
FEMatrix& D,
std::vector<double>& gamma,
std::vector<double>& alphas,
bool compute_D_and_gamma,
const ParameterDatabase& db,
//std::vector<double>& exact_interpolant,
//double& d_h_error,
Limiter limiter,
Iteration_Scheme it_scheme,
const int is_not_afc_fixed_point_rhs)
{
  Output::print<4>("AFC: enter steady_state_algorithm");
  //catch non-square matrix
  if (!system_matrix.is_square())
  {
    ErrThrow("System matrix must be square for AFC!");
  }

  // store the total number of dofs
  int nDofs = system_matrix.GetN_Rows();

  // heritage style index declaration
  int i,j,j0,j1,j2,j3,jj,index;
  double alpha_ij;

  // get pointers to columns, rows and entries of matrix A
  const int * ColInd = system_matrix.GetKCol();
  const int * RowPtr = system_matrix.GetRowPtr();
  // non-const, matrix entries get modified!
  double* Entries = system_matrix.GetEntries();
  double* D_Entries=D.GetEntries();

  int N_Entries = system_matrix.GetN_Entries();
  // allocate memory for matrix F and flux limiters
  double* F = new double[N_Entries+6*nDofs];
  memset(F, 0, (N_Entries+6*nDofs)*SizeOfDouble);
  double* P_plus = F + N_Entries;
  double* P_minus = P_plus + nDofs;
  double* Q_plus = P_minus + nDofs;
  double* Q_minus = Q_plus + nDofs;
  double* R_plus = Q_minus + nDofs;
  double* R_minus = R_plus + nDofs;
  
  std::vector<double> umin, umax, q;
  
 
  if ((it_scheme == Iteration_Scheme::NEWTON)||
    (it_scheme == Iteration_Scheme::FIXEDPOINT_MATRIX)
    ||(it_scheme == Iteration_Scheme::NEWTON_REGU)
    ||(it_scheme == Iteration_Scheme::FIXEDPOINT_RHS))
    alphas.resize(N_Entries+1,0.0);
  // compute entries of the artificial diffusion matrix D
  // TODO make matrix D an actual TMatrix and not only an entries vector
  if (compute_D_and_gamma)
  {
    Output::print<4>("AFC: compute matrix D");
    compute_artificial_diffusion_matrix(system_matrix, D);
    if (limiter == Limiter::BJK17)
    {
      Output::print<4>("AFC: compute vector gamma");
      Compute_Parameter_For_Linearity_Preservation(system_matrix, sol, gamma);
    }
  }
  // add this matrix to A giving \tilde A (Entries)
  // this is the matrix with the properties of an M matrix
  //for fixed_point_rhs and iteration>1 we don't need to add matrix D again.
  if (is_not_afc_fixed_point_rhs)
    system_matrix += D;
  /*Previous Implementation
    Daxpy(N_Entries, 1.0, &afc_matrix_D_entries[0], Entries);*/
  
  // allocate and fill arrays for linearity preserving limiter
  // from Barrenechea, John, Knobloch M3AS (2017)

  if (limiter == Limiter::BJK17)
  {
    umin.resize(nDofs,0.0);
    umax.resize(nDofs,0.0);
    q.resize(nDofs,0.0);

    for (i=0;i<nDofs;i++)
    {
      umin[i] = umax[i] = sol[i];
      q[i] = 0;
      // i-th row of sqmatrix
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      // check values of the neighbor dofs
      for(j=j0;j<j1;j++)
      {
        // column
        index = ColInd[j];
        if (sol[index] < umin[i])
          umin[i] = sol[index];
        if (sol[index] > umax[i])
          umax[i] = sol[index];
        if (i != index)
          q[i] += D_Entries[j];
      }
      q[i] *= gamma[i];
    }
  }

  // compute matrix F
  // loop over all rows
  for(i=0;i<nDofs;i++)
  {
    // i-th row of sqmatrix
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];

    for(j=j0;j<j1;j++)
    {
      // column
      index = ColInd[j];
      // d_ij (u_j - u_i)
      F[j] = D_Entries[j] * (sol[index]-sol[i]);
    }
  }
  // matrix F is computed

  // compute flux limiters
  // loop over all rows
  // linearity preserving limiter from [BJK17]
  if (limiter == Limiter::BJK17)
  {
    for(i=0;i<nDofs;i++)
    {
      // i-th row of sqmatrix
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      // loop over the neighbors
      for(j=j0;j<j1;j++)
      {
        // column
        index = ColInd[j];
        // diagonal
        if (index==i)
          continue;
        if (F[j] > 0)
          P_plus[i] += F[j];
        if (F[j] < 0)
          P_minus[i] += F[j];
      }
      Q_plus[i] = q[i]*(sol[i]-umax[i]);
      Q_minus[i] = q[i]*(sol[i]-umin[i]);
    }
  }

  // Kuzmin limiter
  if (limiter == Limiter::KUZMIN)
  {
    for(i=0;i<nDofs;i++)
    {
      // i-th row of sqmatrix
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      for(j=j0;j<j1;j++)
      {
        // VJ: removed these lines 18/05/09 since they affect 
        // just diagonal entries
        // if ((it_scheme == Iteration_Scheme::FIXEDPOINT_RHS) 
        //      && (Entries[j] > 0))
        // continue;
        // column
        index = ColInd[j];
        // check transposed entry -> jj
        // diagonal
        if (index==i)
          continue;
        j2 = RowPtr[index];
        j3 = RowPtr[index+1];
        for (jj=j2;jj<j3;jj++)
        {
          if (ColInd[jj]==i)
          {
            break;
          }
        }
        // check upwind condition
        // this ensures that the 'link' between i and index is treated only once
        // note that a[jj] > a[j] <==> a[jj]+d[jj] > a[j]+d[j] 
        //  since d is a symmetric matrix
        if (Entries[jj] > Entries[j])
          continue;
        // only the active part of the matrix
        if (F[j] > 0)
        {
          P_plus[i] += F[j];
          if (index<nDofs)
            Q_plus[index] += F[j];
          Q_minus[i] -= F[j];
        }
        if (F[j] < 0)
        {
          P_minus[i] += F[j];
          Q_plus[i] -= F[j];
          if (index<nDofs)
            Q_minus[index] +=  F[j];
        }
      }                                           // end loop j
    }
  }

  // apply the nodal correction factor evaluated at the upwind node i
  // loop over all nodes

  // original but discontinuous proposal
  // other propsals in MooNMD
  for(i=0;i<nDofs;i++)
  {
    // initialization
    R_plus[i] = R_minus[i] = 1.0;
    if (fabs(P_plus[i])>0)
    {
      R_plus[i] = Q_plus[i]/P_plus[i];
      if (R_plus[i] >1)
        R_plus[i] = 1;
    }
    if (fabs(P_minus[i])>0)
    {
      R_minus[i] = Q_minus[i]/P_minus[i];
      if (R_minus[i] >1)
        R_minus[i] = 1;
    }
  }

  // treat Dirichlet nodes
  for (j=0;j < (int) neum_to_diri.size();j++)
  {
    i=neum_to_diri[j];
    R_plus[i] = 1;
    R_minus[i] = 1;
  }

  // apply the flux limiters
  // loop over all rows
  for(i=0;i<nDofs;i++)
  {
    int diag_index=i;
    // i-th row of sqmatrix
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];
    // loop over the columns
    for(j=j0;j<j1;j++)
    {
      // column
      index = ColInd[j];
      // diagonal entry
      if (index==i)
      {
        diag_index = j;
        continue;
      }
      // this should not be happen
      if (Entries[j] > 0)
      {
        ErrThrow("positive non-diagonal entry in FEM-TVD ", i, " ", j, " ",
          Entries[j]);
      }
      if ((it_scheme == Iteration_Scheme::FIXEDPOINT_RHS) && (Entries[j] > 0))
      {
        OutPut("positive entry in FEMTVD "<<i<< " " <<j<<" "<<Entries[j]<<endl);
        exit(4711);
      }

      // original, symmetric application
      if (limiter == Limiter::BJK17)
      {
        alpha_ij = 1;
        if (F[j] > 0)
        {
          alpha_ij = R_plus[i];
          if (R_minus[index] < alpha_ij)
            alpha_ij = R_minus[index];
        }
        if (F[j] < 0)
        {
          alpha_ij = R_minus[i];
          if (R_plus[index] < alpha_ij)
            alpha_ij = R_plus[index];
        }

        // explicit treatment
        if (it_scheme == Iteration_Scheme::FIXEDPOINT_RHS)
        {
          // update rhs
          alphas[j]=alpha_ij;
          rhs[i] += alpha_ij*F[j];          
        }
        if ((it_scheme==Iteration_Scheme::NEWTON)
            || (it_scheme == Iteration_Scheme::FIXEDPOINT_MATRIX)
            || (it_scheme==Iteration_Scheme::NEWTON_REGU))
          // store limiter
          alphas[j]=alpha_ij;
      }

      if (limiter == Limiter::KUZMIN)
      {
        // check transposed entry
        j2 = RowPtr[index];
        j3 = RowPtr[index+1];
        for (jj=j2;jj<j3;jj++)
        {
          // index of transposed entry
          if (ColInd[jj]==i)
          {
            break;
          }
        }
        // check upwind condition
        // this ensures that the 'link' between i and index is treated only once
        if (Entries[jj] > Entries[j])
          continue;

        alpha_ij = 1;
        if (F[j] > 0)
        {
          alpha_ij = R_plus[i];
        }
        if (F[j] < 0)
        {
          alpha_ij = R_minus[i];
        }

        // explicit treatment
        if (it_scheme == Iteration_Scheme::FIXEDPOINT_RHS)
        {
          alphas[j]=alpha_ij;
          alphas[jj]=alpha_ij;
          F[j] = alpha_ij *F[j];
          // update rhs of current row
          rhs[i] += F[j];          
          // update rhs wrt to current column
          // note that F[j] = -F[jj] and 
          // alpha_j = alpha_jj (symmetry of alpha matrix)
          if (index<nDofs)
            rhs[index] -= F[j];
        }
        if ((it_scheme == Iteration_Scheme::NEWTON)
            ||(it_scheme == Iteration_Scheme::FIXEDPOINT_MATRIX)
            || (it_scheme == Iteration_Scheme::NEWTON_REGU))
        {
          // store limiters
          alphas[j] = alpha_ij;
          alphas[jj] = alpha_ij;
        }
      }                                           //End of Kuzmin limiter
    }                                             //End of For loop j
    if((it_scheme==Iteration_Scheme::NEWTON)
        ||(it_scheme == Iteration_Scheme::FIXEDPOINT_MATRIX)
        || (it_scheme == Iteration_Scheme::NEWTON_REGU)
          ||(it_scheme == Iteration_Scheme::FIXEDPOINT_RHS))
      alphas[diag_index]=1.0;
  }                                               //End of loop i

                                               

  if (it_scheme == Iteration_Scheme::FIXEDPOINT_MATRIX)
  {
    int j3, j4, index_j;
    double tau0 =  (double)db["afc_fixed_point_matrix_weight"];
    double tau1 = 1.0-tau0;
    // update matrix
    for(int i=0;i<nDofs;i++)
    {
      j3=RowPtr[i];
      j4=RowPtr[i+1];
      // loop over the columns of the matrix
      for(int j=j3;j<j4;j++)
      {
        index_j=ColInd[j];
        //Non-Diagonal Entries
        if(i!=index_j)
        {
          Entries[j]-=alphas[j]*D_Entries[j]*tau0;
          // if (Entries[j] > 0)
          //  Output::print<2>("non-diag pos ", Entries[j]);
          rhs[i] += alphas[j]*D_Entries[j]*sol[index_j]*tau1;
        }
        else
        {
          double sum_flux=0.0;
          for(int jj=j3;jj<j4;jj++)
          {
            if(ColInd[jj]!=i)
              sum_flux += alphas[jj]*D_Entries[jj];
          }
          Entries[j]+=sum_flux*tau0;
          //if (Entries[j] < 0)
          //  Output::print<2>("diag neg ", Entries[j]);

          rhs[i] -= sum_flux*sol[i]*tau1;
        }
      }
    }                                             //Formation of matrix complete
  }

  if ((it_scheme == Iteration_Scheme::NEWTON)
      || (it_scheme == Iteration_Scheme::NEWTON_REGU))
  {
    int j3, j4, index_j;
    std::vector<double> df(N_Entries,0.0);
    double tau0 = (double)db["afc_fixed_point_matrix_weight"];
    double tau1 = (double)db["afc_fixed_point_derivative_weight"]
                 *(double)db["afc_fixed_point_derivative_weight_factor"];
    if (db["afc_iteration_scheme"].is("newton_no_damp"))
    {
      if ((double)db["afc_fixed_point_derivative_weight_factor"]>1e-3)
      {
        tau0 =  tau1 = 1.0;
      }
    }
    Output::print<4>("tau0/tau1 ", tau0 , " " , tau1);
    
    // compute first part of rhs
    // with old matrix
    for(int i=0;i<nDofs;i++)
    {
      j3=RowPtr[i];
      j4=RowPtr[i+1];
      for(int j=j3;j<j4;j++)
      {
        index_j = ColInd[j];
        rhs[i] += -(Entries[j])*sol[index_j]+alphas[j]*F[j];
      }
    }

    Output::print<4>("AFC: computing Jacobian");
    // compute Jacobian, store on matrix entries
    if (limiter == Limiter::KUZMIN)
    {
      if (it_scheme == Iteration_Scheme::NEWTON)
      {
        // Computation of Matrix DF
        // loop over the degrees of freedom
        for(int i=0;i<nDofs;i++)
        {
          j3=RowPtr[i];
          j4=RowPtr[i+1];
          // loop over the columns of the matrix
          for(int j=j3;j<j4;j++)
          {
            index_j=ColInd[j];
            //Non-Diagonal Entries
            if(i!=index_j)
            {
              //Derivative Product is a function which returns the summation 
              //inside the DF matrix
              df[j]=Entries[j]-tau0*alphas[j]*D_Entries[j]
                -tau1*Compute_Jacobian_times_flux_Kuzmin(system_matrix,
                                                         F,P_plus,Q_plus,
                                                         Q_minus,P_minus,
                                                         R_minus,R_plus,i,
                                                         index_j,j,
                                                         D);
            }
            else
            {
              double sum_flux=0.0;
              for(int jj=j3;jj<j4;jj++)
              {
                if(ColInd[jj]!=i)
                  sum_flux += alphas[jj]*D_Entries[jj];
              }
              df[j]=Entries[j]+tau0*sum_flux
                -tau1*Compute_Jacobian_times_flux_Kuzmin(system_matrix,
                                                         F,P_plus,Q_plus,
                                                         Q_minus,P_minus,
                                                         R_minus,R_plus,i,
                                                         index_j,j,
                                                         D);
            }
          }               //end of loop j
        }                 //end of loop i     //Formation of matrix DF complete
      }                   //end of Newton case

      if (it_scheme == Iteration_Scheme::NEWTON_REGU)
      {
        // Computation of Matrix DF
        // loop over the degrees of freedom
        double sigma = (double)db["afc_newton_regu_sigma"];
        
        for(int i=0;i<nDofs;i++)
        {
          j3=RowPtr[i];
          j4=RowPtr[i+1];
          // loop over the columns of the matrix
          for(int j=j3;j<j4;j++)
          {
            index_j=ColInd[j];
            // Entries contains already a+d
            df[j]=Entries[j]
             -Compute_Jacobian_times_flux_Kuzmin_Regularized(system_matrix,
                                                             F,i,index_j,j,
                                                             D,
                                                             sigma, tau0, tau1);
          }
        }                //end of loop i      //Formation of matrix DF complete
      }                  //end of Newton_regu
    }                    //end of KUZMIN limiter case

    if(limiter == Limiter::BJK17)
    {
      if (it_scheme == Iteration_Scheme::NEWTON)
      {
        //Computation of matrix DF
        //Loop over the degrees of freedom
        for(int i=0;i<nDofs;i++)
        {
          j3=RowPtr[i];
          j4=RowPtr[i+1];
          //Loop over the columns
          for(int j=j3;j<j4;j++)
          {
            index_j=ColInd[j];
            //Non-Diagonal Entries
            if(i!=index_j)
            {
              //Derivative Product is a function which returns the summation inside the DF matrix
              df[j]=Entries[j]
                   -tau0*alphas[j]*D_Entries[j]
                   -tau1*Compute_Jacobian_times_flux_BJK17(system_matrix, 
                                                           F, P_plus,P_minus,
                                                           Q_plus,Q_minus,
                                                           R_plus,R_minus,
                                                           umax,umin,q,
                                                           sol[index_j],i,
                                                           index_j,j,
                                                           D);
            }
            else
            {
              double sum_flux=0.0;
              for(int jj=j3;jj<j4;jj++)
              {
                if(ColInd[jj]!=i)
                  sum_flux +=alphas[jj]*D_Entries[jj];
              }
              df[j]=Entries[j]+tau0*sum_flux
                   -tau1*Compute_Jacobian_times_flux_BJK17(system_matrix, 
                                                           F, P_plus,P_minus,
                                                           Q_plus,Q_minus,
                                                           R_plus,R_minus,
                                                           umax,umin,q,
                                                           sol[index_j],i,
                                                           index_j,j,
                                                           D);
            }
          }                //end of loop j
        }                  //end of loop i  //Formation of matrix DF complete
      }                    //end of Newton case
    }                      //end of BJK17 limiter case

    // update rhs for system to be solved
    for(int i=0;i<nDofs;i++)
    {
      j3=RowPtr[i];
      j4=RowPtr[i+1];
      for(int j=j3;j<j4;j++)
      {
        index_j = ColInd[j];
        rhs[i] += df[j]*sol[index_j];           //Using DF instead of matrix A
        Entries[j]=df[j];
      }
    }
  }

  delete [] F;                                    
}                                               //End of function definition

/** ************************************************************************ */

void AlgebraicFluxCorrection::crank_nicolson_fct(
const FEMatrix& M_C, FEMatrix& K,
const std::vector<double>& oldsol,
std::vector<double>& rhs, std::vector<double>& rhs_old,
double delta_t,
const std::vector<int>& neum_to_diri,
AlgebraicFluxCorrection::Prelimiter prelim
)
{

  //make sure Crank-Nicolson is used
  if( TDatabase::TimeDB->TIME_DISC != 2 )
  {
    ErrThrow("crank_nicolson_fct performs Crank-Nicolson time discretization."
      "Change TDatabase::TimeDB->TIME_DISC to 2.")
  }

  //check if both matrices are square and dimensions match
  if(!M_C.is_square() || !K.is_square())
  {
    ErrThrow("M_C and K must be square!");
  }
  if(M_C.GetN_Rows() != K.GetN_Rows())
  {
    ErrThrow("M_C and K dimension mismatch!")
  }

#ifdef _MPI
  const TParFECommunicator3D& comm = K.GetFESpace3D()->get_communicator();
  const int* masters = comm.GetMaster();
  int size, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  int nDofs = K.GetN_Rows();

  // get pointers to columns, rows and entries of matrix M_C
  const int* ColInd_M = M_C.GetKCol();
  const int* RowPtr_M = M_C.GetRowPtr();
  const double* Entries_M = M_C.GetEntries();

  // get pointers to columns, rows and entries of matrix K; plus number of
  // entries.
  const int* ColInd = K.GetKCol();
  const int* RowPtr = K.GetRowPtr();
  int N_Entries = K.GetN_Entries();
  
  FEMatrix D(K);
  
  D.reset();
  
  double* afc_matrix_D_Entries=D.GetEntries();

  //STEP 0.1: Compute artificial diffusion matrix D and add to K, K := K + D ("L").
  compute_artificial_diffusion_matrix(K, D);
  
  K+=D;

  //MPI: Matrix K is in Level-0-consistency now.

  //STEP 0.2: Lump mass matrix.
  std::vector<double> lump_mass (nDofs, 0.0);
  lump_mass_matrix_to_vector(M_C, lump_mass);

  //STEP 1: Computation of the intermediate solution
  // follows Kuzmin(2009), S.10 (44)

  double theta = 0.5;          //hard-coded implicitness. 1/2 for Crank-Nicolson
  double one_minus_theta = 0.5;                  //hard-coded expliciteness

  //approximation to u(t + /frac{1}{2} /delta t)
  std::vector<double> u_interm(nDofs,0.0);
  //approximation to u'(t + /frac{1}{2} /delta t)
  std::vector<double> u_dot_interm(nDofs,0.0);

  // step-by-step build up u_interm...
  // ...first L u_{k-1}
  // MPI: oldsol is normally given in Level2-Consistency,
  // when this method is called
  K.multiply(&oldsol[0], &u_interm[0]);

  // ...then M_L^(-1)(f_{k-1} - L u_{k-1})
  for(int i=0;i<nDofs;i++)
  {
#ifdef _MPI
    // do this only for master rows
    // u_interm gets Level-0-consistency now
    if(masters[i] != rank)
      continue;
#endif
    u_interm[i] = (rhs_old[i] - u_interm[i])/lump_mass[i];
  }

  // ...then u_{k-1} + halftimestep * M_L^(-1)(f_{k-1} - L u_{k-1})
  double halftimestep = delta_t/2.0;
  for(int i=0;i<nDofs;i++)
  {
    u_interm[i] = oldsol[i] + halftimestep * u_interm[i];
  }

  // set non-actives for intermediate solution (to those of interm. solution)
  for ( int j=0 ; j < (int) neum_to_diri.size() ; j++)
  {
    int i=neum_to_diri[j];
    // just copy the values from old solution - this should
    // have the non-actives set correctly!
    u_interm[i] = oldsol[i];
  }

  // compute u_dot_interm (Kuzmin 2009 S.9 (37))
  for (int i = 0; i < nDofs ; ++i)
  {
    u_dot_interm[i] = 2*(u_interm[i] - oldsol[i])/delta_t;
  }
  //intermediate solution is complete

  // MPI: u_interm and u_dot_interm are both in Level-0-consistency
  // and should get level 2 for the next loop
#ifdef _MPI
  comm.consistency_update(u_interm.data(),2);
  comm.consistency_update(u_dot_interm.data(),2);
#endif

  //STEP 2: compute raw antidiffusive fluxes * delta_t, and apply pre-limiting,
  // if required to do so (Kuzmin 2009 S.9 (36) )
  std::vector<double> raw_fluxes(N_Entries, 0.0);
  for(int i=0;i<nDofs;i++)
  {
#ifdef _MPI
    // do this only for master rows
    if(masters[i] != rank)
      continue;
#endif
    // i-th row of mass matrix
    int j0 = RowPtr_M[i];
    int j1 = RowPtr_M[i+1];

    for(int j=j0;j<j1;j++)
    {
      // column
      int index = ColInd_M[j];

      if (index==i)                 //nothing to do if this is a diagonal entry
        continue;

      double val = - afc_matrix_D_Entries[j] * (u_interm[i]-u_interm[index]);
      raw_fluxes[j]= Entries_M[j] * (u_dot_interm[i]-u_dot_interm[index]) + val;

      // pre-limit fluxes
      switch(prelim)
      {
        case Prelimiter::NONE:
          // no prelimiting
          break;
        case Prelimiter::BOTH:
          // no break statement, will execute both following prelimiters!
        case Prelimiter::MIN_MOD:
          raw_fluxes[j] = MinMod(raw_fluxes[j]/delta_t, val);
          break;
        case Prelimiter::GRAD_DIRECTION:
          if (raw_fluxes[j]*(u_interm[index]-u_interm[i])>0)
          {
            raw_fluxes[j] = 0;
          }
          break;
      }
      //multiply with delta, required by our Zalesak implementation
      raw_fluxes[j] *= delta_t;
    }
  }

  //STEP 3: apply Zalesaks flux correction to gain correction factors alpha
  std::vector<double> alphas( M_C.GetN_Entries(), 0.0 );
  ZalesaksFluxLimiter( alphas, M_C, lump_mass, raw_fluxes,
    u_interm, neum_to_diri);

  // STEP 4: Calculate the new right hand side and system matrix.
  for(int i=0;i<nDofs;i++)
  {
#ifdef _MPI
    // do this only for master rows
    if(masters[i] != rank)
      continue;
#endif
    double corrected_flux = 0.0;
    int j0 = RowPtr[i];
    int j1 = RowPtr[i+1];
    for(int j=j0;j<j1;j++)
    {
      int index = ColInd[j];
      if(i != index)
      {
        corrected_flux +=  alphas[j] * raw_fluxes[j];
      }
    }
    //do the right hand side updates
    rhs_old[i]=rhs[i];
    rhs[i] = u_interm[i] * lump_mass[i] + theta*delta_t*rhs[i] +
      corrected_flux;
  }

  //and finally update K to be the final system matrix
  fem_fct_compute_system_matrix(K, lump_mass, delta_t,
    theta, one_minus_theta);

}

/** ************************************************************************ */

void AlgebraicFluxCorrection::correct_dirichlet_rows(FEMatrix& MatrixA)
{
  //hold pointers to row, kcol, entries array
  const int* RowPtr_A      = MatrixA.GetRowPtr();
  const int* KCol_A        = MatrixA.GetKCol();
  double* Entries_A  = MatrixA.GetEntries();

  //determine first and one-after-last dirichlet rows
  size_t diriHighBound;
  size_t diriLowBound;
#ifdef __3D__
  diriHighBound = MatrixA.GetFESpace3D()->GetDirichletBound();
  diriLowBound = diriHighBound - MatrixA.GetFESpace3D()->GetN_Dirichlet();
#elif __2D__
  diriHighBound = MatrixA.GetFESpace2D()->GetDirichletBound();
  diriLowBound = diriHighBound - MatrixA.GetFESpace2D()->GetN_Dirichlet();
#endif

  // loop over rows and set them to unity-vectors
  for (size_t rowIndex = diriLowBound;
    rowIndex < diriHighBound ;++rowIndex)
  {
    int l0 = RowPtr_A[rowIndex];
    int l1 = RowPtr_A[rowIndex+1];
    for (int l=l0;l<l1;l++)
    {
      // diagonal entry
      if (KCol_A[l]== (int) rowIndex)
        Entries_A[l] = 1;
      else
        Entries_A[l] = 0;
    }
  }
}

/** ************************************************************************ */

void AlgebraicFluxCorrection::AFC_Compute_New_Iterate(
  const BlockVector& old_solution,
  BlockVector& new_solution,
  const ParameterDatabase& db)
{
  new_solution.add_scaled(old_solution,-1.0);
  new_solution.scale(db["afc_nonlinloop_damping_factor"]);
  new_solution.add_scaled(old_solution,1.0);
  new_solution.copy_nonactive(old_solution);

  //Note: Projection only possible when the maximum and minimum are already known.
  if (db["afc_project_to_admissible_values"].is("yes"))
  {
    Output::print<2>("  projection to admissible values");

    int len = new_solution.length();
    for (int i=0;i<len;i++)
    {
      if (new_solution[i]<0)
        new_solution[i] = 0;
      if (new_solution[i]>1)
        new_solution[i] = 1;
    }
  }
}
