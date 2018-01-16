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

}

ParameterDatabase AlgebraicFluxCorrection::default_afc_database()
{
  ParameterDatabase db("default algebraic flux correction database");

  // Type of AFC to be applied.
  db.add("algebraic_flux_correction", "none", " Chose which type of afc to use.",
         {"none", "default", "fem-tvd", "fem-fct-cn"});

  db.add("afc_prelimiter", 0, "Chose an afc flux prelimiting scheme. Options "
      "are 0 (none), 1 (min-mod), 2 (grad-direction), 3 (both)", {0,1,2,3});

  return db;
}

// ////////////////////////////////////////////////////////////////////////////
// Implementation of the methods in the namespace AlgebraicFluxCorrection    //
// ////////////////////////////////////////////////////////////////////////////


void AlgebraicFluxCorrection::fem_tvd_algorithm(
    FEMatrix& system_matrix,
    const std::vector<double>& sol,
    std::vector<double>& rhs,
    //this argument is used in outcommented code block only
    const std::vector<int>& neum_to_diri,
    bool continuous_proposal,
    bool nonsymmetric_application)
{
  //catch non-square matrix
  if (!system_matrix.is_square())
  {
    ErrThrow("System matrix must be square for FEM-TVD!");
  }

  // store the total number of dofs
  int nDofs = system_matrix.GetN_Rows();

  // heritage style index declaration
  int i,j,j0,j1,j2,j3,jj,index;
  double nenner, zaehler;

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

  // compute entries of the artificial diffusion matrix D
  // TODO make matrix D an actual TMatrix and not only an entries vector
  std::vector<double> matrix_D_Entries(N_Entries, 0.0);
  compute_artificial_diffusion_matrix(system_matrix, matrix_D_Entries);

  // add this matrix to A giving \tilde A (Entries)
  // this is the matrix with the properties of an M matrix
  Daxpy(N_Entries, 1.0, &matrix_D_Entries[0], Entries);

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
      // d_ij (u_i - u_j)
      F[j] = matrix_D_Entries[j] * (sol[index]-sol[i]);
    }
  }
  // matrix F is computed

  // compute flux limiters
  // loop over all rows
  for(i=0;i<nDofs;i++)
  {
    // i-th row of sqmatrix
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];
    for(j=j0;j<j1;j++)
    {
      if (Entries[j] > 0)
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
    }                                             // end loop j
  }

  // apply the nodal correction factor evaluated at the upwind node i
  // loop over all nodes
  if (!continuous_proposal)
  { // original but discontinuous proposal
    for(i=0;i<nDofs;i++)
    {
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
  }
  else
  { // continuous proposal
    for(i=0;i<nDofs;i++)
    {
      zaehler =  Q_plus[i];
      if (-Q_minus[i] < zaehler)
        zaehler = -Q_minus[i];
      nenner = 1e-32;
      if (P_plus[i] > nenner)
        nenner = P_plus[i];
      if (-P_minus[i] > nenner)
        nenner = -P_minus[i];
      R_plus[i] = zaehler/nenner;
      if (R_plus[i] > 1)
        R_plus[i] = 1;
      R_minus[i] = R_plus[i];
    }
  }


//  // treat Dirichlet nodes CB 2016/ I commented this out,
    // because it breaks my example. TODO Needs further investigation.
//  for (j=0;j < (int) neum_to_diri.size();j++)
//  {
//    i=neum_to_diri[j];
//    if (i >= nDofs)
//    {
//      ErrThrow("neum_to_diri index ", i, " is out of scope!");
//    }
//    R_plus[i] = 1;
//    R_minus[i] = 1;
//  }

  // apply the flux limiters
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
      if (index==i)
        continue;
       // this should not happen
      if (Entries[j] > 0)
      {
        ErrThrow("positive non-diagonal entry in FEM-TVD ", i, " ", j, " ",
                 Entries[j]);
      }

      // check transposed entry
      j2 = RowPtr[index];
      j3 = RowPtr[index+1];
      for (jj=j2;jj<j3;jj++)
      {
        if (ColInd[jj]==i)
        {
          break;
        }
      }

      if (!nonsymmetric_application)
      {
        // original, symmetric application
        // check upwind condition
        // this ensures that the 'link' between i and index is treated only once
        if (Entries[jj] > Entries[j])
          continue;
        // compute contribution to rhs
        if (F[j] > 0)
        {
          F[j] = R_plus[i]*F[j];
        }
        else
          F[j] = R_minus[i]*F[j];
        // update rhs of current row
        rhs[i] += F[j];
        // update rhs wrt to current column
        // note that F[j] = -F[jj] and alpha_j = alpha_jj (symmetry)
        if (index<nDofs)
          rhs[index] -= F[j];
      }
      else
      { // nonsymmetric application
        // compute contribution to rhs
        if (F[j] > 0)
        {
          F[j] = R_plus[i]*F[j];
        }
        else
          F[j] = R_minus[i]*F[j];
        // update rhs of current row
        rhs[i] += F[j];
      }
    }
  }
  delete [] F;
}

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
