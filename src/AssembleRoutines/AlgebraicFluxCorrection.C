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
namespace {
 

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

  if(!(int) lumped_mass_matrix.size() == nDof)
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
      {//diagonal entry

        check_cfl_condition(lumped_mass_matrix[i], Entries[j], tau, theta2);

        // calculate the new system matrix entry
        Entries[j] = lumped_mass_matrix[i] + tau * theta1 * Entries[j];

      }
      else
      {//non-diagonal entry
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
    const FEMatrix& A, std::vector<double>& matrix_D)
{
  // catch non-square matrix
  if (!A.is_square() )
  {
    ErrThrow("Matrix must be square!");
  }
  // store number of dofs
  int nDof = A.GetN_Rows();

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
  for(int i=0;i<nDof;i++) //row loop
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
      int j = ColInd[l]; //column index
      double k_ij = 0;
      double k_ji = 0;
      if (j!=i)
      {// consider only off-diagonals
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
        matrix_D[l] = std::min({-k_ij , 0.0 , -k_ji});
      }
    }
  }

  // compute diagonal entries of the matrix D
  for(int i=0;i<nDof;i++)//row loop
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
      val +=  matrix_D[l];
      int j = ColInd[l];
      if (j==i) //diagonal found
        ll = l; //place of the diagonal entry in the matrix_D entries array
    }
    matrix_D[ll] = -val;
  }
}

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
    return a; //return sign * a;
  }
  else
  {
    return b; //return sign * b;
  }
}
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

  if(! (int) lumped_mass.size() == nDofs)
  {
    ErrThrow("lumped_mass has incorrect size ", lumped_mass.size(),
             " != ", nDofs);
  }


  if(! (int) sol_approx.size() == nDofs)
  {
    ErrThrow("sol_approx has incorrect size ", sol_approx.size(),
             " != ", nDofs);
  }

  if(! (int) alphas.size() == nEntries)
  {
    ErrThrow("alphas has incorrect size ", alphas.size(), " != ", nEntries);
  }

  if(! (int) raw_fluxes.size() == nEntries)
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

      if(i != index) //diagonal entries are skipped
      {
        if (raw_fluxes[j] > 0.0)
          P_plus[i] += raw_fluxes[j]; //P_plus was initialized with 0
        if (raw_fluxes[j] <= 0.0)
          P_minus[i] += raw_fluxes[j]; //P_minus was initialized with 0

        double help = sol_approx[index]-sol_approx[i];

        if (help > Q_plus[i]) //Q_plus was initialized with 0
          Q_plus[i]= help;

        if (help < Q_minus[i]) //Q_minus was initialized with 0
          Q_minus[i]= help;
      }
    } // end loop j
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
    }                                             //end loop j
  }

  // MPI: At output, correction factors alpha_ij are rowwise level-0-consistent
}

/* Computation of Derivative Product for the Matrix DF(Jacobian Matrix for Newton iteration)
 *            and the Zalesak limiter
 * @param[in] A: The system matrix A
 * @param[in] F: The Matrix where entries f_ij=d_ij-(uj_-u_i)
 * @param[in] P_plus, Q_plus
 *            P_minus, Q_minus
 *            R_plus, R_minus: The fluxes used for computation of the limiters alpha_ij
 * @param[in] row_i: Denotes the ith row of matrix DF for which summation is required
 * @param[in] afc_matrix_D_entries: The artificial diffusion matrix
 * @param[in] entries_pointer: Denotes the pointer to the array Entries 
 * @param[in] col_j: Column number for entry a_ij
 * 
 * @param[out] The summation required in the entry of the Jacobian matrix.
 */ 
double Derivative_Product_Zalesak(const FEMatrix& A,
                          const double * F , 
			  const double * P_plus,  
			  const double * Q_plus,
			  const double * Q_minus,
			  const double * P_minus,
			  const double * R_minus,
			  const double * R_plus,
			  int row_i,
			  const std::vector<double>& afc_matrix_D_entries,
			  int entries_pointer, 
			  int col_j)
{
  double sum=0;
  //DQ_minus, DQ_plus, DP_minus, DP_plus: Denotes the Derivative of the raw fluxes.
  double DQ_plus, DQ_minus,DP_plus,DP_minus, sumD=0;
  const double eps = 1e-10;
  const int* ColInd = A.GetKCol();
  const int* RowPtr = A.GetRowPtr(); 
  const double* Entries = A.GetEntries();
  int k0, k1, temp, temp1;
  
  k0=RowPtr[row_i];
  k1=RowPtr[row_i+1];  
  //Computation of the entry a_ji
  for(int mm=RowPtr[col_j];mm<RowPtr[col_j+1];mm++)
  {
    if(ColInd[mm]==row_i)
    {
      temp1=mm;
      break;
      
    }
  }
  // loop over all dofs that are connected with the matrix entry a_{ij}
  for(int k=k0; k<k1;k++)
  {
    // cases without contribution to the sum
    if(F[k]>0 && fabs(R_plus[row_i]-1)<eps)
      continue;
    else if(fabs(F[k])<eps)
      continue;
    else if(F[k]<0 && fabs(R_minus[row_i]-1)<eps)
      continue;    
    // first case with contribution
    else if(F[k]>0 && R_plus[row_i]<1)
    {
      //Calculations for finding derivative of Q_plus
      if(F[entries_pointer]>=0 && row_i!=col_j)
	DQ_plus=0;
      else if(F[entries_pointer]<0 && row_i!=col_j)
        DQ_plus=-afc_matrix_D_entries[entries_pointer];
      else if(row_i==col_j)
      {
	sumD=0;
	int m0,m1;
	m0=RowPtr[row_i];
	m1=RowPtr[row_i+1];
	for(int m=m0;m<m1;m++)
	{
	  if(F[m]<0)
	   sumD+=afc_matrix_D_entries[m];
        }//End of loop m
	DQ_plus=sumD;
      }
      //Calculations for finding derivative of P_plus
      if(F[entries_pointer]<=0 && row_i!=col_j)
	DP_plus=0;
      else if(F[entries_pointer]>0 && row_i!=col_j && Entries[entries_pointer]>=Entries[temp1])
        DP_plus=afc_matrix_D_entries[entries_pointer];
      else if(row_i==col_j)
      {
	sumD=0;
	int m0,m1;
	m0=RowPtr[row_i];
	m1=RowPtr[row_i+1];
	for(int m=m0;m<m1;m++)
	{
	 int transposed_entry=ColInd[m];
	 for(int ll=RowPtr[transposed_entry];ll<RowPtr[transposed_entry+1];ll++)
	 {
	   if(ColInd[ll]==row_i)
	   {
	    temp=ll;
	    break;
	   }
         }
	 if(F[m]>0 && Entries[m]>=Entries[temp])
	  sumD+=afc_matrix_D_entries[m];
        }//End of loop m
	DP_plus=-sumD;
      }
      sum+=((DQ_plus*P_plus[row_i]-Q_plus[row_i]*DP_plus)/(P_plus[row_i]*P_plus[row_i]))*F[k];
      continue;
     }
     // second case with contribution
      else if(F[k]<0 && R_minus[row_i]<1)
     {
      //Calculations for finding derivative of Q_minus
       if(F[entries_pointer]<=0 && row_i!=col_j)
	DQ_minus=0;
      else if(F[entries_pointer]>0 && row_i!=col_j)
        DQ_minus=-afc_matrix_D_entries[entries_pointer];
      else if(row_i==col_j)
      {
	sumD=0;
	int m0,m1;
	m0=RowPtr[row_i];
	m1=RowPtr[row_i+1];
	for(int m=m0;m<m1;m++)
	{
	  if(F[m]>0)
	  sumD+=afc_matrix_D_entries[m];
	}//End of loop m
	DQ_minus=sumD;
      }
      //Calculations for finding derivative of P_minus
      if(F[entries_pointer]>=0 && row_i!=col_j)
	DP_minus=0;
      else if(F[entries_pointer]<0 && row_i!=col_j && Entries[entries_pointer]>=Entries[temp1])
        DP_minus=afc_matrix_D_entries[entries_pointer];
      else if(row_i==col_j)
      {
	sumD=0;
        int m0,m1;
	m0=RowPtr[row_i];
	m1=RowPtr[row_i+1];
	for(int m=m0;m<m1;m++)
	{
	 int transposed_entry=ColInd[m];
	 for(int ll=RowPtr[transposed_entry];ll<RowPtr[transposed_entry+1];ll++)
	 {
	   if(ColInd[ll]==row_i)
	   {
	    temp=ll;
	    break;
	   }
          }
          if(F[m]<0 && Entries[m]>=Entries[temp])
	   sumD+=afc_matrix_D_entries[m];
	  }//End of loop m
	DP_minus=-sumD;
      }
      sum+=((DQ_minus*P_minus[row_i]-Q_minus[row_i]*DP_minus)/(P_minus[row_i]*P_minus[row_i]))*F[k];
      continue;
     } 
     else
       Output::print<4>("ERROR ERROR");
  }
return sum;  
}

/* Computation of Derivative Product for the Matrix DF(Jacobian Matrix for Newton iteration)
 *            and the BJK17 Linear preserving limiter
 * @param[in] A: The system matrix A
 * @param[in] F: The Matrix where entries f_ij=d_ij-(uj_-u_i)
 * @param[in] P_plus, Q_plus
 *            P_minus, Q_minus
 *            R_plus, R_minus: The fluxes used for computation of the limiters alpha_ij
 * @param[in] umax: Vector where umax_i=max{u_ij, 1<=j<=Ndofs}
 * @parma[in] umin: Vector where umin_i=min{u_ij, 1<=j<=Ndofs}
 * @param[in] sol: The solution vector of the problem
 * @parma[in] q: Vector where q_i=gamma_i\sum d_ij
 * @param[in] row_i: Denotes the ith row of matrix DF for which summation is required
 * @param[in] afc_matrix_D_entries: The artificial diffusion matrix
 * @param[in] entries_pointer: Denotes the pointer to the array Entries 
 * @param[in] col_j: Column number for entry a_ij
 * 
 * @param[out] The summation required in the entry of the Jacobian matrix.
 */ 

double Derivative_Product_Linear(const FEMatrix A,
			  const double * F , 
			  const double * P_plus,  
			  const double * P_minus,
			  const double * Q_plus,
			  const double * Q_minus,
			  const double * R_plus,
			  const double * R_minus,
			  const std::vector<double> umax,
			  const std::vector<double> umin,
			  const std::vector<double> q,
			  const std::vector<double>& sol,
			  const int row_i,
			  const int entries_pointer, 
			  const int col_j,
			  const std::vector<double> afc_matrix_D_entries)
                                 
{
  double sum=0;
  //DQ_minus, DQ_plus, DP_minus, DP_plus: Denotes the Derivative of the raw fluxes.
  double DQ_plus, DQ_minus,DP_plus,DP_minus, sumD=0;
  const int* ColInd = A.GetKCol();
  const int* RowPtr = A.GetRowPtr(); 
  int k0, k1;  
  k0=RowPtr[row_i];
  k1=RowPtr[row_i+1];  
  for(int k=k0;k<k1;k++)
  {
    //Conditions for which we don't have any contributions
    if(F[k]==0)
      continue;
    else if(F[k]>0 && R_plus[row_i]==1)
      continue;
    else if(F[k]<0 && R_minus[row_i]==1)
      continue;
    //First conditon where we have contribution
    else if(F[k]>0 && R_plus[row_i]<1)
    {
      if(R_plus[row_i]<=R_minus[ColInd[k]])
      {
	//Computation of derivatives of Q_plus
	if(row_i!=col_j)
	{
	  if(umax[row_i]==sol[col_j])
	    DQ_plus=-q[row_i];
	  else
	    DQ_plus=0;
	}
	else
	{
	  if(umax[row_i]==sol[col_j])
	    DQ_plus=0;
	  else
	    DQ_plus=q[row_i];
	}
	//Computation of derivatives of P_plus
	if(F[entries_pointer]<=0 && row_i!=col_j)
	DP_plus=0;
        else if(F[entries_pointer]>0 && row_i!=col_j)
        DP_plus=afc_matrix_D_entries[entries_pointer];
        else if(row_i==col_j)
        {
	 sumD=0;
	 int m0,m1;
	 m0=RowPtr[row_i];
	 m1=RowPtr[row_i+1];
	 for(int m=m0;m<m1;m++)
	 {
	   if(F[m]>0)
	    sumD+=afc_matrix_D_entries[m];
         }//End of loop m
	 DP_plus=-sumD;
        }
       sum+=((DQ_plus*P_plus[row_i]-Q_plus[row_i]*DP_plus)/(P_plus[row_i]*P_plus[row_i]))*F[k];
       continue;
      }
      //R_minus[k]<R_plus[i] then R_minus[k]<1
      else
      {
	 if(ColInd[k]!=col_j)
	  {
	    if(umin[ColInd[k]]==sol[col_j])
	    DQ_minus=-q[ColInd[k]];
	    else
	    DQ_minus=0;
	  }
	  else
	  {
	    if(umin[ColInd[k]]==sol[col_j])
	     DQ_minus=0;
	    else
	     DQ_minus=q[ColInd[k]];
	  }
	  if(F[entries_pointer]>=0 && ColInd[k]!=col_j)
	   DP_minus=0;
          else if(F[entries_pointer]<0 && ColInd[k]!=col_j)
           DP_minus=afc_matrix_D_entries[entries_pointer];
          else if(ColInd[k]==col_j)
          {
	    sumD=0;
	    int m0,m1;
	    m0=RowPtr[ColInd[k]];
	    m1=RowPtr[ColInd[k]+1];
	    for(int m=m0;m<m1;m++)
	    {
	      if(F[m]<0)
	       sumD+=afc_matrix_D_entries[m];
            }//End of loop m
	    DP_plus=-sumD;
          }
          sum+=((DQ_minus*P_minus[ColInd[k]]-Q_minus[ColInd[k]]*DP_minus)/(P_minus[ColInd[k]]*P_minus[ColInd[k]]))*F[k];
          continue;
      }
    }
    else if(F[k]<0 && R_minus[row_i]<1)
    {
      if(R_minus[row_i]<=R_plus[ColInd[k]])
      {
	if(row_i!=col_j)
	{
	  if(umin[row_i]==sol[col_j])
	    DQ_minus=-q[row_i];
	  else
	    DQ_minus=0;
	}
	else
	{
	  if(umin[row_i]==sol[col_j])
	    DQ_minus=0;
	  else
	    DQ_minus=q[row_i];
	}
	if(F[entries_pointer]>=0 && row_i!=col_j)
	DP_minus=0;
        else if(F[entries_pointer]<0 && row_i!=col_j)
        DP_minus=afc_matrix_D_entries[entries_pointer];
        else if(row_i==col_j)
        {
	 sumD=0;
	 int m0,m1;
	 m0=RowPtr[row_i];
	 m1=RowPtr[row_i+1];
	 for(int m=m0;m<m1;m++)
	 {
	   if(F[m]<0)
	    sumD+=afc_matrix_D_entries[m];
         }//End of loop m
	 DP_minus=-sumD;
        }
       sum+=((DQ_minus*P_minus[row_i]-Q_minus[row_i]*DP_minus)/(P_minus[row_i]*P_minus[row_i]))*F[k];
       continue;
      }
      else
      {
	  if(ColInd[k]!=col_j)
	  {
	    if(umax[ColInd[k]]==sol[col_j])
	    DQ_plus=-q[ColInd[k]];
	    else
	    DQ_plus=0;
	  }
	  else
	  {
	    if(umax[ColInd[k]]==sol[col_j])
	     DQ_plus=0;
	    else
	     DQ_plus=q[ColInd[k]];
	  }
	  if(F[entries_pointer]<=0 && ColInd[k]!=col_j)
	   DP_plus=0;
          else if(F[entries_pointer]>0 && ColInd[k]!=col_j)
           DP_plus=afc_matrix_D_entries[entries_pointer];
          else if(ColInd[k]==col_j)
          {
	    sumD=0;
	    int m0,m1;
	    m0=RowPtr[ColInd[k]];
	    m1=RowPtr[ColInd[k]+1];
	    for(int m=m0;m<m1;m++)
	    {
	      if(F[m]>0)
	       sumD+=afc_matrix_D_entries[m];
            }//End of loop m
	    DP_plus=-sumD;
          }
          sum+=((DQ_plus*P_plus[ColInd[k]]-Q_plus[ColInd[k]]*DP_plus)/(P_plus[ColInd[k]]*P_plus[ColInd[k]]))*F[k];
          continue;
      }
    }
    else
      Output::print<4>("ERROR");
  }
  return sum;
}

/**
 * Compute the weights for the linearty preserving limiter from Barrenechea, John, Knobloch; M3AS 2017.
 *
 * @param[in] FEMatrix system matrix
 * @param[in] current solution vector
 * @param[out] gamma vector with the weights 
 * 
 */

void Compute_Parameter_For_Linearity_Preservation(FEMatrix& system_matrix, const std::vector<double>& u,
  std::vector<double>& gamma)
{
 
  int i, j, index, N_Cells, N_Unknowns, N_V;
  int *global_numbers, *begin_index, *dof;
  double area, edge, x[3], y[3], dist, dist_max, dist_max2;
  TCollection *Coll;
  const TFESpace2D *fespace;
  TBaseCell *cell;
  FE2D CurrentElement;

  Output::print<4>("AFC: compute parameter for linearity preservation");
  // get fe space  
  fespace = system_matrix.GetFESpace2D();  
  // get collection
  Coll = fespace->GetCollection();
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
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    // get pointer to local dofs
    dof = global_numbers+begin_index[i];
    // finite element on the mesh cell
    CurrentElement = fespace->GetFE2D(i, cell);
    // only for P_1
    if (CurrentElement != C_P1_2D_T_A)
    {
      ErrThrow("Compute_Parameter_For_Linearity_Preservation only implemented for P_1 !!!");
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
}
 
}

ParameterDatabase AlgebraicFluxCorrection::default_afc_database()
{
  ParameterDatabase db("default algebraic flux correction database");

  // Type of AFC to be applied.
  db.add("algebraic_flux_correction", "none", " Chose which type of afc to use.",
         {"none", "afc", "fem-fct-cn"});

  db.add("afc_prelimiter", 0, "Choose an afc flux prelimiting scheme. Options "
      "are 0 (none), 1 (min-mod), 2 (grad-direction), 3 (both)", {0,1,2,3});
  
  // type of the limiter to be applied for steady-state case
  db.add("afc_limiter", "zalesak", "Choose an afc limiter. Options are"
      "zalesak, linearity_preserving_BJK17", {"zalesak", "linearity_preserving_BJK17"});
  
  // iteration scheme for the afc methods 
  db.add("afc_iteration_scheme", "fixed_point_matrix", "Choose an iteration scheme for the afc methods. Options are"
      "fixed_point_rhs, fixed_point_matrix, newton", {"fixed_point_rhs", "fixed_point_matrix", "newton"});
  
  db.add("afc_nonlinloop_damping_factor", 1.0, "A damping parameter for the nonlinear loop in AFC."
    "Must be a value between 1 (no damping) and 0 (no update).", 0.0,1.0);

  db.add("afc_nonlinloop_maxit", (size_t) 1 , "Maximal number of iterations for the nonlinear loop in AFC."
    "Must be a value between 0 and 100000.", (size_t)  0, (size_t)  100000);

  db.add("afc_nonlinloop_epsilon", 1e-10, "Stopping criterion for the nonlinear loop in AFC."
    "Must be a value between 1e-20 an 1e20.", 1e-20, 1e20);
  
  return db;
}

// ////////////////////////////////////////////////////////////////////////////
// Implementation of the methods in the namespace AlgebraicFluxCorrection    //
// ////////////////////////////////////////////////////////////////////////////


void AlgebraicFluxCorrection::steady_state_algorithm(
    FEMatrix& system_matrix,
    const std::vector<double>& sol,
    std::vector<double>& rhs,
    //this argument is used in outcommented code block only
    const std::vector<int>& neum_to_diri,
    std::vector<double>& afc_matrix_D_entries,
    std::vector<double>& gamma,
    bool compute_D_and_gamma,
    Limiter limiter,
    Iteration_Scheme it_scheme)
{
  Output::print<4>("AFC: enter steady_state_algorithm");
  //catch non-square matrix
  if (!system_matrix.is_square())
  {
    ErrThrow("System matrix must be square for FEM-TVD!");
  }

  // store the total number of dofs
  int nDofs = system_matrix.GetN_Rows();

  // heritage style index declaration
  int i,j,j0,j1,j2,j3,jj,kk,index;
  double alpha_ij;

  // get pointers to columns, rows and entries of matrix A
  const int * ColInd = system_matrix.GetKCol();
  const int * RowPtr = system_matrix.GetRowPtr();
  // non-const, matrix entries get modified!
  double* Entries = system_matrix.GetEntries();

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
  std::vector<double> alphas;
  std::vector<double> umin, umax, q;

  if (it_scheme == Iteration_Scheme::NEWTON)
     alphas.resize(N_Entries,0.0);
  // compute entries of the artificial diffusion matrix D
  // TODO make matrix D an actual TMatrix and not only an entries vector  
  if (compute_D_and_gamma)
  {
    Output::print<4>("AFC: compute matrix D");
    compute_artificial_diffusion_matrix(system_matrix, afc_matrix_D_entries);
    if (limiter == Limiter::LINEAR_PRESERVE_BJK17)
    {
      Output::print<4>("AFC: compute vector gamma");
      Compute_Parameter_For_Linearity_Preservation(system_matrix, sol, gamma);
    }
  }

  // add this matrix to A giving \tilde A (Entries)
  // this is the matrix with the properties of an M matrix
  Daxpy(N_Entries, 1.0, &afc_matrix_D_entries[0], Entries);
  
  // allocate and fill arrays for linearity preserving limiter
  // from Barrenechea, John, Knobloch M3AS (2017)
  //if ((TDatabase::ParamDB->FEM_FCT_LINEAR_TYPE==4)||(TDatabase::ParamDB->FEM_FCT_LINEAR_TYPE==14))
  if (limiter == Limiter::LINEAR_PRESERVE_BJK17)
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
          q[i] += afc_matrix_D_entries[j];
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
      F[j] = afc_matrix_D_entries[j] * (sol[index]-sol[i]);
    }
  }
  // matrix F is computed

  // compute flux limiters
  // loop over all rows
  // linearity preserving limiter from [BJK17]
  if (limiter == Limiter::LINEAR_PRESERVE_BJK17)
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
   
  // Zalesak limiter 
  if (limiter == Limiter::ZALESAK)
  {
    for(i=0;i<nDofs;i++)
    {
      // i-th row of sqmatrix
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      for(j=j0;j<j1;j++)
      {
	if ((it_scheme == Iteration_Scheme::FIXEDPOINT_RHS) && (Entries[j] > 0))
          continue;
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
    double update_diag;
    
    // i-th row of sqmatrix
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];
    update_diag = 0;
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
        OutPut("positive entry in FEMTVD " << i << " " << j << " " << Entries[j] << endl);
        exit(4711);
      }
      
        // original, symmetric application
        if (limiter == Limiter::LINEAR_PRESERVE_BJK17)
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
          if (F[j]==0)
	    alpha_ij = 1;
	  
         // explicit treatment
	  if (it_scheme == Iteration_Scheme::FIXEDPOINT_RHS) 
	 {
	    // update rhs
            rhs[i] += alpha_ij*F[j];
	 }
	  if (it_scheme == Iteration_Scheme::FIXEDPOINT_MATRIX) 
	 { // implicit treatment
	   // off diagonal
	   Entries[j] += (0-alpha_ij)*afc_matrix_D_entries[j];
	   update_diag -= (0-alpha_ij)*afc_matrix_D_entries[j];
	 }
	 if(it_scheme==Iteration_Scheme::NEWTON)
	   alphas[j]=alpha_ij;
        }
        
        if (limiter == Limiter::ZALESAK)
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
	  
	  if (F[j] > 0)
          {
	    alpha_ij = R_plus[i];
          }
          else
          {
	    alpha_ij = R_minus[i];
          }
          if (F[j]==0)
	    alpha_ij = 1;
 
	 // explicit treatment   
	  if (it_scheme == Iteration_Scheme::FIXEDPOINT_RHS) 
	 {
            F[j] = alpha_ij *F[j];
            // update rhs of current row
            rhs[i] += F[j];
            // update rhs wrt to current column
            // note that F[j] = -F[jj] and alpha_j = alpha_jj (symmetry of alpha matrix)
            if (index<nDofs)
              rhs[index] -= F[j];
	 }
	  if (it_scheme == Iteration_Scheme::FIXEDPOINT_MATRIX) 
	 {
	   // implicit treatment
	   // both entries are change with the same magnitude since matrix_D_Entries is symmetric
	   Entries[j] += (0-alpha_ij)*afc_matrix_D_entries[j];
	   update_diag -= (0-alpha_ij)*afc_matrix_D_entries[j];
	   Entries[jj] += (0-alpha_ij)*afc_matrix_D_entries[jj];
  	   // find diagonal entry of the row from the transposed entry
           for (kk=j2;kk<j3;kk++)
           {
             if (ColInd[kk]==index)
             {
               break;
             }
           }
	   Entries[kk] -= (0-alpha_ij)*afc_matrix_D_entries[jj];
	 }
	 if(it_scheme == Iteration_Scheme::NEWTON)
	 {
	   alphas[j] = alpha_ij;
	   alphas[jj] = alpha_ij;
	 }
	}//End of Zalesak limiter
    } //End of For loop j

    // update diagonal 
         if (it_scheme == Iteration_Scheme::FIXEDPOINT_MATRIX) 
         {
            //OutPut(i << " update_diag " << update_diag << endl);
            Entries[diag_index] += update_diag;
	    //Output::print<4>(i," ", j, " ", diag_index, " ",  Entries[diag_index]);
         }
         if(it_scheme==Iteration_Scheme::NEWTON)
	   alphas[diag_index]=1.0;
    
  }//End of loop i
 
 
  if(it_scheme == Iteration_Scheme::NEWTON)
  {
    int j3, j4, index1;
    std::vector<double> df(N_Entries,0.0);    
    // compute first part of rhs
    // with old matrix 
    for(int i=0;i<nDofs;i++)
    {
      j3=RowPtr[i];
      j4=RowPtr[i+1];
      for(int j=j3;j<j4;j++)
      {
	index1 = ColInd[j];
	rhs[i] -= Entries[j]*sol[index1]+(0-alphas[j])*F[j]-sol[i]*afc_matrix_D_entries[j];
        //Output::print<4>(index1, " ", sol[index1],  " " , alphas[j], " " , F[j]);
      } 
    }
    
    Output::print<4>("AFC: computing Jacobian");
    // compute Jacobian, store on matrix entries
    if (limiter == Limiter::ZALESAK)
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
	  index1=ColInd[j];
	  //Non-Diagonal Entries
	  if(i!=index1) 
	  {
	    //Derivative Product is a function which returns the summation inside the DF matrix
	    df[j]=Entries[j]+(0-alphas[j])*afc_matrix_D_entries[j]
	    -Derivative_Product_Zalesak(system_matrix,F,P_plus,Q_plus,Q_minus,P_minus,R_minus,R_plus,i,afc_matrix_D_entries,j, index1);
          }
	  else
	  {
	    double sum_flux=0.0;
	    for(int jj=j3;jj<j4;jj++)
	      {
	        if(ColInd[jj]!=i)
	        sum_flux+=(0-alphas[jj])*afc_matrix_D_entries[jj];
	      }
	    df[j]=Entries[j]-sum_flux-Derivative_Product_Zalesak(system_matrix,F,P_plus,Q_plus,Q_minus,P_minus,R_minus,R_plus,i,afc_matrix_D_entries,j, index1);
	  }
        }
      }//Formation of matrix DF complete
    }
    
    if(limiter == Limiter::LINEAR_PRESERVE_BJK17)
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
	  index1=ColInd[j];
	  //Non-Diagonal Entries
	  if(i!=index1) 
	  {
	    //Derivative Product is a function which returns the summation inside the DF matrix
	    df[j]=Entries[j]+(0-alphas[j])*afc_matrix_D_entries[j]
	    -Derivative_Product_Linear(system_matrix, F, P_plus,P_minus,Q_plus,Q_minus,R_plus,R_minus,umax,umin,q,sol,i,j,index1,afc_matrix_D_entries);
          }
	  else
	  {
	    double sum_flux=0.0;
	    for(int jj=j3;jj<j4;jj++)
	      {
	        if(ColInd[jj]!=i)
	        sum_flux+=(0-alphas[jj])*afc_matrix_D_entries[jj];
	      }
	    df[j]=Entries[j]-sum_flux-Derivative_Product_Linear(system_matrix, F, P_plus,P_minus,Q_plus,Q_minus,R_plus,R_minus,umax,umin,q,sol,i,j,index1,afc_matrix_D_entries);
	  }
        }
      }
    }
        
    // update rhs for system to be solved
    for(int i=0;i<nDofs;i++)
    {
      j3=RowPtr[i];
      j4=RowPtr[i+1];
      for(int j=j3;j<j4;j++)
      {
	index1 = ColInd[j];
	rhs[i] += df[j]*sol[index1]; //Using DF instead of matrix A
	//Output::print<4>(index1, " ", sol[index1],  " " , alphas[j], " " , F[j]);
	Entries[j]=df[j];
      }
    } 
  }

  delete [] F; // here it should be
}//End of function definition

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
  double* Entries = K.GetEntries();
  int N_Entries = K.GetN_Entries();

  //STEP 0.1: Compute artificial diffusion matrix D and add to K, K := K + D ("L").
  std::vector<double> matrix_D_Entries(N_Entries, 0.0);
  compute_artificial_diffusion_matrix(K, matrix_D_Entries);

  Daxpy(N_Entries, 1.0, &matrix_D_Entries[0], Entries);

  //MPI: Matrix K is in Level-0-consistency now.

  //STEP 0.2: Lump mass matrix.
  std::vector<double> lump_mass (nDofs, 0.0);
  lump_mass_matrix_to_vector(M_C, lump_mass);


  //STEP 1: Computation of the intermediate solution
  // follows Kuzmin(2009), S.10 (44)

  double theta = 0.5; //hard-coded implicitness. 1/2 for Crank-Nicolson
  double one_minus_theta = 0.5; //hard-coded expliciteness

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

      if (index==i) //nothing to do if this is a diagonal entry
        continue;

      double val = - matrix_D_Entries[j] * (u_interm[i]-u_interm[index]);
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

void AlgebraicFluxCorrection::AFC_Compute_New_Iterate(const BlockVector& old_solution, BlockVector& new_solution,
  const ParameterDatabase& db)
{
    {
      // update direction
      new_solution.add_scaled(old_solution,-1.0);
      new_solution.scale(db["afc_nonlinloop_damping_factor"]);
      new_solution.add_scaled(old_solution,1.0);
      new_solution.copy_nonactive(old_solution);
    }
}
 