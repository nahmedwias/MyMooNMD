/**
 * Implementation of class MumpsWrapper, declared in MumpsWrapper.h
 *
 * @author Clemens Bartsch
 * @date 2016/03/14
 */
#ifdef _MPI

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <MumpsWrapper.h>
#include <ParFECommunicator3D.h>

#include <mpi.h>
#include <memory>

//define macros for mumps integer parameters
#define JOB_INIT -1
#define JOB_END -2
#define ICNTL(I) icntl[(I)-1] //macro doing the fortran shift for control params
#define INFOG(I) infog[(I)-1]

MumpsWrapper::MumpsWrapper(
    const BlockFEMatrix& bmatrix, std::vector<const TParFECommunicator3D*> comms)
{
  //input checks
  if (comms.size() != bmatrix.get_n_cell_rows())
  {
    throw std::runtime_error("Number of given communicators does not match"
        " block order.");
  }
  size_t n_comms = comms.size();

  for(size_t index = 0 ; index<n_comms; ++index) //check if fespaces match
  {
    //compare adresses of the FESpaces
    if(comms.at(index)->get_fe_space() != &bmatrix.get_row_space(index))
      throw std::runtime_error("Adresses of spaces do not match.");
  }
  for(size_t index = 0; index<n_comms; ++index)
  // check if all communicators have dimension 1 - this is a restriction which
  // must be removed later!
  {
    if(comms.at(index)->get_n_dim() != 1 )
      throw std::runtime_error("MumpsWrapper can only be used with "
          "ParFeCommunicators of dimension 1");
  }
  //todo these input checks are makeshift and will be adapted later when the
  //class takes clearer shape (and eventually Mapper & Communicator got reworked)

  // initialize the mumps entity
  // set "hard" parameters
  id_.par = 1; // root used in computation
  id_.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
  id_.sym = 0; //non-symmetric matrix
  id_.job=JOB_INIT;
  dmumps_c(&id_);
  // the "softer" parameters must be set after initializing
  set_mumps_parameters();

  // transform and store the matrix distributed in coordinate format
  store_in_distributed_coordinate_form(bmatrix, comms);

}

int MumpsWrapper::solve(
    const BlockVector& rhs, BlockVector& solution,
    std::vector<TParFECommunicator3D*> comms)
{

  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int n_masters_local_comms;
  int n_dofs_global_comms;
  check_input_solve(rhs, solution, comms, n_masters_local_comms, n_dofs_global_comms);

  //1) set up the global right hand side in root
  int root_rank = 0;
  bool i_am_root = (mpi_rank == root_rank);

  std::vector<double> master_values;
  master_values.reserve(n_masters_local_comms); //reserve space for local rhs and local solution
  std::vector<double> rhs_global(n_dofs_global_comms, 0); //used in root only, initialized everywhere FIXME replace by dummy calls

  //three different shifts in this loop, all due to arrangement in blocks
  // - for reading the local rhs "block_read_shift", which
  //   accounts for all dofs (master, slave, halo)
  // - for writing the local rhs "block_write_shift", which
  //   accounts for master dofs only
  // - for writing the global right hand side "global_write_shift", which
  //   accounts for the global number of dofs per block (determined by adding up master dofs)
  size_t loc_master_shift = 0;
  size_t loc_dof_shift = 0;
  size_t glob_dof_shift = 0;

  for(size_t index = 0; index < comms.size() ;++index) //loop over blocks
  {
    //fill the local right hand side with master dof rows only
    const int* masters = comms.at(index)->GetMaster(); //TODO this should be a vector (in ParFECommunicator)!
    size_t n_loc_dofs_block = comms.at(index)->GetNDof();
    for(size_t i = 0; i< n_loc_dofs_block; ++i)
    {//push rhs values for all master dofs on this rank and block into rhs2sol_local_
      if(masters[i] == mpi_rank)
      {
        master_values.push_back(rhs.at(loc_master_shift + i));
      }
    }
    //...and gather those local right hand sides globally, maintaining the local2global dof ordering

//    //mpi calls
//    //determine how many values to send per process - n_masters_per_proc
//    size_t block_n_masters_local = comms.at(index)->GetN_Master();
//    int* n_masters_per_proc = new int[mpi_size];
//    MPI_Allgather(&block_n_masters_local, 1, MPI_INT, //send
//                  n_masters_per_proc, 1, MPI_INT,     //receive
//                  MPI_COMM_WORLD);                    //control
//    //determine the shifts in the receiving vector
//    //- this is what takes care of maintaining the global dof order!
//    int* recv_shift = new int[mpi_size];
//    recv_shift[0] = 0;
//    for(int i=1;i< mpi_size; ++i)
//      recv_shift[i] = recv_shift[i-1] + n_masters_per_proc[i-1];
//    //now fill the current block's portion in the global rhs vector
//    double* temp_glob = &rhs_global.at(global_write_shift);
//    double* temp_loc = &master_values.at(block_write_shift);
//    MPI_Gatherv(temp_loc, block_n_masters_local, MPI_DOUBLE, //send
//                temp_glob, n_masters_per_proc, recv_shift, MPI_DOUBLE, //receive
//                root_rank, MPI_COMM_WORLD); //control

    int n_loc_masters_block = comms.at(index)->GetN_Master();
    comms.at(index)->GatherToRoot(
        &rhs_global.at(glob_dof_shift),        //receive
        &master_values.at(loc_master_shift), n_loc_masters_block,  //send
        root_rank);                                                //control

    // one block treated - count up the shifts
    loc_master_shift += n_loc_masters_block;
    loc_dof_shift += n_loc_dofs_block;
    int current_n_dofs_global;
    MPI_Allreduce(&n_loc_masters_block,&current_n_dofs_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    glob_dof_shift += current_n_dofs_global;
  }

  //2) let mumps do its jobs
  //2.1) give the data to mumps...
  // ..the (distributed) matrix
  id_.nz_loc  = matrix_.nz_loc;
  id_.irn_loc = &matrix_.irn_loc.at(0);
  id_.jcn_loc = &matrix_.jcn_loc.at(0);
  id_.a_loc   = &matrix_.a_loc.at(0);
  //...the (centralized) rhs
  if(i_am_root)
  {
    id_.rhs  = &rhs_global.at(0);
    id_.nrhs = 1;
    id_.lrhs = matrix_.n;
    id_.n    = matrix_.n;
  }

  //2.2) kick off the job
  kick_off_job(std::string("analyze"));
  kick_off_job(std::string("factorize"));
  kick_off_job(std::string("solve"));

  //3) distribute solution among processes (?)
  //three different shifts in this loop, all due to arrangement in blocks
  // - "block_dof_shift", which accounts for all dofs (master, slave, halo)
  // - "block_master_shift", which accounts for master dofs only
  // - "global_shift", which accounts for the global number of dofs
  //    per block (determined by adding up master dofs)
  loc_dof_shift = 0;
  loc_master_shift = 0;
  glob_dof_shift = 0;

  for(size_t index = 0; index < comms.size(); ++index) //loop over blocks
  {
    //receive all master values from root and store in master_values
    int n_loc_masters_block = comms.at(index)->GetN_Master();
    int n_glob_dofs_block;
    MPI_Allreduce(&n_loc_masters_block,&n_glob_dofs_block,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    comms.at(index)->ScatterFromRoot(
        &rhs_global.at(glob_dof_shift),                           //send
        &master_values.at(loc_master_shift), n_loc_masters_block, //receive
        root_rank);                                               //control

    //write master values into local solution vector
    const int* masters = comms.at(index)->GetMaster(); //TODO this should be a vector (in ParFECommunicator)!
    int n_local_dofs_block = comms.at(index)->GetNDof();
    int i_master = 0;
    for(int i_dof = 0; i_dof< n_local_dofs_block; ++i_dof)
    {//write solution values for all master dofs on this rank into solution vector (at correct place)
      if(masters[i_dof] == mpi_rank)
      {
        solution.at(loc_dof_shift + i_dof) = master_values.at(loc_master_shift + i_master);
        ++i_master;
      }
    }

    // Update all non-master values in solution vector. Fire!
    comms.at(index)->CommUpdate(&solution.at(loc_dof_shift));

    // count up the block shifts
    int n_loc_dofs_block = comms.at(index)->GetNDof();
    loc_dof_shift += n_loc_dofs_block;
    loc_master_shift += n_loc_masters_block;
    glob_dof_shift += n_glob_dofs_block;
  }
  return 0;
}

void MumpsWrapper::write_matrix_distributed(std::string filename) const
{
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  std::ofstream matrixfile;
  std::string file = filename + std::to_string(mpi_rank);
  matrixfile.open(file.c_str());

  //write the header line - coordinate format, real values, no symmetry used
  matrixfile << "%%MatrixMarket matrix coordinate real general \n";

  //write general matrix information
  matrixfile << matrix_.n << "\t" << matrix_.n << "\t" << matrix_.nz_loc << "\n";

  // loop and write info to file
  for (size_t index = 0; index < matrix_.nz_loc; ++index)
  {
    matrixfile << matrix_.irn_loc.at(index) << "\t" << matrix_.jcn_loc.at(index) << "\t"
        << matrix_.a_loc.at(index) << "\n";
  }
}


// Special member functions.
MumpsWrapper::MumpsWrapper(MumpsWrapper&&)
{
  //todo
  ;
}

MumpsWrapper& MumpsWrapper::operator=(MumpsWrapper&&)
{
  //todo
  ;
}

MumpsWrapper::~MumpsWrapper()
{
  id_.job=JOB_END;
  dmumps_c(&id_);
}

// Private functions.

void MumpsWrapper::check_input_solve(const BlockVector& rhs, const BlockVector& solution,
                                     const std::vector<TParFECommunicator3D*>& comms,
                                     int& n_masters_local_comms,
                                     int& n_dofs_global_comms)
{
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  //0) input checks
   n_masters_local_comms = 0;
   n_dofs_global_comms = 0;
   for (size_t i =0; i < comms.size(); ++i)
   {
     if(comms.at(i)->get_n_dim() != 1 )
       throw std::runtime_error("MumpsWrapper can only be used with "
           "ParFECommunicators of dimension 1");

     size_t local_n_dofs = comms.at(i)->GetNDof();
     size_t block_length_rhs = rhs.length(i);
     size_t block_length_sol = solution.length(i);

     if(local_n_dofs != block_length_rhs)
     {
       ErrThrow("Length of rhs block ", i, " does not fit n of local dofs in given communicator.");
     }
     if(local_n_dofs != block_length_sol)
     {
       ErrThrow("Length of sol block ", i, " does not fit n of local dofs in given communicator.");
     }
     // add up all masters
     n_masters_local_comms += comms.at(i)->GetN_Master();
   }
   MPI_Allreduce(&n_masters_local_comms, &n_dofs_global_comms, 1, MPI_INT,MPI_SUM, MPI_COMM_WORLD);
   if(n_dofs_global_comms != (int) matrix_.n)
   {
     ErrThrow("Total number of masters of given communicators not equal global stored matrix order.");
   }

}

void MumpsWrapper::kick_off_job(const std::string& job)
{
  int job_id = 0;
  //determine job
  if(job.compare("analyze") == 0)
    job_id = 1;
  else if(job.compare("factorize") == 0)
    job_id = 2;
  else if(job.compare("solve") == 0)
    job_id = 3;
  else
    ErrThrow("The string '",  job ,"' does not describe a MUMPS job!");

  // run the chosen job
  id_.job = job_id;
  dmumps_c(&id_);

  //inform about eventual errors
  if(id_.INFOG(1)<0)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool i_am_root = (rank == 0);
    if(i_am_root)
    {
      Output::print("MUMPS Analyze failed INFOG(1) = ",  id_.INFOG(1));
      Output::print("MUMPS Analyze failed INFOG(2) = ",  id_.INFOG(2));
    }
    MPI_Finalize();
    exit(-1);
  }

}

void MumpsWrapper::set_mumps_parameters()
{
  //input format parameters
  id_.ICNTL(5)  = 0; // matrices in "assembled" format
  id_.ICNTL(18) = 3; // structure AND entries are distributed among processes
  id_.ICNTL(20) = 0; // dense right hand side (stored in nrhs and lrhs)
  id_.ICNTL(21) = 0; //try centralized solution first!

  // parameters collected from old MumpsSolver.C
  id_.ICNTL(1) = 1;    // standard outstream for errors
  id_.ICNTL(2) = -1;   // no warnings output
  id_.ICNTL(3) = -1;   // no global info output
  id_.ICNTL(4) = 1;    // verbosity level
  //id_.ICNTL(6) = 1;  // matrix permutation control??
  //id_.ICNTL(7) = 5;  // matrix permutation control?
  id_.ICNTL(14)=400;   //estimated working space increase (%)
}

void MumpsWrapper::store_in_distributed_coordinate_form(
    const BlockFEMatrix& bmatrix,
    std::vector<const TParFECommunicator3D*> comms
)
{
  // no input checks - assume that this method is only called from
  // the constructor, which already performed checks

  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // 0) preparation
  // 0.0) determine global number of dofs by adding up the masters
  size_t n_dofs_global = 0; //sum over all local masters is number of global dofs
  {
    size_t n_masters_local = 0;
    for(auto comm : comms)
    {
      n_masters_local += comm->GetN_Master();
    }
    MPI_Allreduce(&n_masters_local, &n_dofs_global, 1,MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }

  // 0.1)
  //determine an array which holds index-shifts
  size_t nComms = comms.size(); //should equal the numbers of block rows/block columns/spaces
  std::vector<size_t> shifts(nComms,0);
  //for each space, find out how many dofs there are in total
  //(over all processors) - adding up those will give the shift
  for(size_t index =0; index< nComms - 1;++index)
  {
    int n_masters_local = comms.at(index)->GetN_Master();
    //sum up all nLocalMasters over all processes and store
    MPI_Allreduce(&n_masters_local, &shifts.at(index + 1), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //add up the shift with the previous one
    shifts.at(index + 1) += shifts.at(index);
  }

  // 0.2)
  // peel the master arrays from the Communicators -make copies!
  std::vector<std::vector<int>> dof_masters(nComms);
  for (size_t index = 0; index<nComms; ++index)
  {
    const int* master_array = comms.at(index)->GetMaster();
    size_t master_array_length = comms.at(index)->GetNDof(); //that correct??? - fuer N_Dim =1 wahrscheinlich schon
    dof_masters.at(index).resize(master_array_length);
    std::copy(master_array,
              master_array + master_array_length, //pointer arithmetic!
              dof_masters.at(index).begin());
  }

  // 0.3)
  // peel the local2global vectors from the communicators - make copies!
  std::vector<std::vector<int>> local2globals(nComms);
  for (size_t index = 0; index<nComms; ++index)
  {
    const int* l2g = comms.at(index)->Get_Local2Global();
    size_t array_length = comms.at(index)->GetNDof(); //that correct??? - fuer N_Dim =1 wahrscheinlich schon
    local2globals.at(index).resize(array_length);
    std::copy(l2g,
              l2g + array_length, //pointer arithmetic!
              local2globals.at(index).begin());
  }

  // 1) setting up the matrix
  // 1.1)
  //get a rough overview and reserve some space
  matrix_.n = n_dofs_global;
  matrix_.nz_loc = 0;
  {
    size_t expected_entries = bmatrix.get_n_total_entries();
    matrix_.irn_loc.reserve(expected_entries);
    matrix_.jcn_loc.reserve(expected_entries);
    matrix_.a_loc.reserve(expected_entries);
  }


  //loop through the block matrix
  //will always hold the transposed state of the last treated block
  bool transp;
  //will always hold the last treated block
  std::shared_ptr<const FEMatrix> block = bmatrix.get_block(0, 0, transp);

  for(size_t cell_row = 0; cell_row < nComms; ++cell_row) //loop over rows
  {
    const std::vector<int>& masters = dof_masters.at(cell_row); //fix masters
    size_t row_shift = shifts.at(cell_row);                     //fix shift
    const std::vector<int>& row_l2g  = local2globals.at(cell_row); //fix local-to-global mapping

    for(size_t cell_col = 0; cell_col < nComms; ++cell_col) //loop over columns
    {
      size_t col_shift = shifts.at(cell_col); //fix shift
      const std::vector<int>& col_l2g  = local2globals.at(cell_col); //fix local-to-global mapping

      //fetch all the info from the current block
      block = bmatrix.get_block(cell_row, cell_col, transp);
      const int* row_ptr = block->GetRowPtr();
      const int* k_col = block->GetKCol();
      const double* entries = block->GetEntries();

      if(!transp) //block is stored non-transposed
      {
        for (size_t row = 0; (int) row < block->GetN_Rows(); ++row)
        {//active rows

          size_t row_start = row_ptr[row];
          size_t row_end = row_ptr[row+1];

          for(size_t k = row_start; k<row_end ;++k)
          {
            size_t col = k_col[k];
            double entry = entries[k];

            // gather the conditions under which (AND) we will to put the entry
            // into the sparsity structure on this process
            bool is_nonzero = (entry != 0.0);
            bool is_in_master_row = masters.at(row) == mpi_rank;
            bool is_in_active_row = row < bmatrix.get_n_row_actives(cell_row);

            if (is_nonzero && is_in_master_row && is_in_active_row)
            {//put entry into the new mumps matrix
              matrix_.irn_loc.push_back(row_shift + row_l2g.at(row) + 1);
              matrix_.jcn_loc.push_back(col_shift + col_l2g.at(col) + 1);
              matrix_.a_loc.push_back(entry);
              ++matrix_.nz_loc;
            }
          }
        }
      }//end non-transposed case
      else if (transp) //block is stored transposed
      {
        if ( block->GetActiveBound() != block->GetN_Rows()) //check for security
        {
          ErrThrow("This block (",cell_row,",",cell_col,") has test space "
                   "non-actives, it should never have been stored in transposed state!");
        }
        for (int row = 0; row < block->GetN_Rows(); ++row)
        {
          size_t row_start = row_ptr[row];
          size_t row_end = row_ptr[row+1];

          for(size_t k = row_start; k<row_end ;++k)
          {
            size_t col = k_col[k];
            double entry = entries[k];

            // gather the conditions under which (AND) we will to put the entry
            // into the sparsity structure on this process
            bool is_nonzero = (entry != 0.0);
            bool is_in_master_row = masters.at(col) == mpi_rank;
            bool is_in_active_row = col < bmatrix.get_n_row_actives(cell_row);

            if (is_nonzero && is_in_master_row && is_in_active_row)
            {//fill entry in the sparsity structure
              matrix_.irn_loc.push_back(row_shift + row_l2g.at(col) + 1);
              matrix_.jcn_loc.push_back(col_shift + col_l2g.at(row) + 1);
              matrix_.a_loc.push_back(entry);
              ++matrix_.nz_loc;
            }
          }
        }
      }
      //on diagonal blocks, set the ones on diagonals of non-active rows
      // (same code both cases, but loop does only make iterations for
      // non-transp case, because transposed storage of blocks with testspace
      // non-actives is not allowed in BlockFEMatrix)
      if(cell_row == cell_col)
      {
        for(int row =   block->GetActiveBound(); row < block->GetN_Rows(); ++row)
        {
          if(masters.at(row)==mpi_rank)
          {//only add the entry if the current process is master of its row
            //put a 1 on the diagonal
            matrix_.irn_loc.push_back(row_shift + row_l2g.at(row) + 1);
            matrix_.jcn_loc.push_back(row_shift + row_l2g.at(row) + 1);
            matrix_.a_loc.push_back(1);
            ++matrix_.nz_loc;
          }
        }
      }//end treating dirichlet rows
    }//end loop over columns
  }//end loop over rows
}

#endif
