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
#include <Database.h>

#include <mpi.h>
#include <memory>

//two of the used mumps job codes
#define JOB_INIT -1
#define JOB_END -2

//macros doing the fortran shift for mumps control and info params
#define ICNTL(I) icntl[(I)-1]
#define INFOG(I) infog[(I)-1]

MumpsWrapper::MumpsWrapper(const BlockFEMatrix& bmatrix,
                           std::vector<double> pres0)
: analyzed_and_factorized(false)
{

  //check the input
  check_input_matrix(bmatrix);

  //copy the BlockFEMatrices communicators and store them
  comms_ = bmatrix.get_communicators();

  // initialize the mumps entity
  // set those "hard" parameters which must be set before initializing
  id_.par = 1; // root will take part in computation
  id_.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
  id_.sym = 0; //non-symmetric matrix
  id_.job=JOB_INIT;
  dmumps_c(&id_);

  // the "softer" parameters must be set after initializing
  set_mumps_parameters();

  // transform and store the matrix distributed in coordinate format
  store_in_distributed_coordinate_form(bmatrix, pres0);

}

MumpsWrapper::MumpsWrapper( const BlockMatrix& bmatrix,
                            std::vector<const TParFECommunicator3D*> comms,
                            std::vector<double> pres0,
                            std::vector<std::vector<int>> loc_to_seq)
: analyzed_and_factorized(false)
{

  //check the input
  check_input_matrix(bmatrix);

  //copy the BlockFEMatrices communicators and store them
  comms_ = comms;

  // initialize the mumps entity
  // set those "hard" parameters which must be set before initializing
  id_.par = 1; // root will take part in computation
  id_.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
  id_.sym = 0; //non-symmetric matrix
  id_.job=JOB_INIT;
  dmumps_c(&id_);

  // the "softer" parameters must be set after initializing
  set_mumps_parameters();

  // transform and store the matrix distributed in coordinate format
  store_in_distributed_coordinate_form(bmatrix, pres0, loc_to_seq);
}

void MumpsWrapper::solve(const BlockVector& rhs, BlockVector& solution)
{
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  int root_rank = 0;
  bool i_am_root = (mpi_rank == root_rank);


  //0) do input check and gather two important numbers
  int n_masters_local_comms;
  int n_dofs_global_comms;
  check_input_solve(rhs, solution, n_masters_local_comms, n_dofs_global_comms);

  //1) set up the global right hand side in root
  std::vector<double> master_values;
  master_values.reserve(n_masters_local_comms); //space for loc rhs and solution

  std::vector<double> rhs_global;
  if(i_am_root) //put up an array for the global right hand side
    rhs_global.resize(n_dofs_global_comms, 0);

  //three different shifts in this loop, all due to arrangement in blocks
  // - "loc_dof_shift", accounts for all local dofs (master, slave, halo)
  // - "loc_master_shift", which accounts for master dofs only
  // - "glob_dof_shift", which accounts for the global number of dofs per block
  //   (determined by adding up master dofs)
  size_t loc_master_shift = 0;
  size_t loc_dof_shift = 0;
  size_t glob_dof_shift = 0;

  for(size_t index = 0; index < comms_.size() ;++index) //loop over blocks
  {
    //fill the local right hand side with master dof rows only
    const int* masters = comms_.at(index)->GetMaster(); //TODO this should be a vector (in ParFECommunicator)!
    size_t n_loc_dofs_block = comms_.at(index)->GetNDof();
    for(size_t i = 0; i< n_loc_dofs_block; ++i)
    {//push rhs values for master dofs on this rank and block to master_values
      if(masters[i] == mpi_rank)
      {
       master_values.push_back(rhs.at(loc_dof_shift + i));
      }
    }

    //...and gather these local right hand sides globally
    int n_loc_masters_block = comms_.at(index)->GetN_Master();
    double* global_rhs_dummy = nullptr;
    if(i_am_root)
    {
      global_rhs_dummy = &rhs_global.at(glob_dof_shift);
    }

    //BUGFIX: make sure that master_values.at(loc_master_shift) does not throw,
    //even if the number of masters on this rank in this block was 0
    if(n_loc_masters_block == 0)
    {
      master_values.push_back(0);
    }

    gather_vector(
        global_rhs_dummy,                                          //receive
        &master_values.at(loc_master_shift), n_loc_masters_block,  //send
        root_rank);                                                //control

    if(n_loc_masters_block == 0)
    {
      master_values.pop_back();  //revert action due to BUGFIX above
    }

    // one block treated - count up the shifts
    loc_master_shift += n_loc_masters_block;
    loc_dof_shift += n_loc_dofs_block;
    int current_n_dofs_global;
    MPI_Allreduce(&n_loc_masters_block,&current_n_dofs_global,1,
                  MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    glob_dof_shift += current_n_dofs_global;
  }
  //2) let mumps do its jobs
  if(!analyzed_and_factorized)
  {
    id_.nz_loc  = matrix_.nz_loc;
    id_.irn_loc = &matrix_.irn_loc.at(0);
    id_.jcn_loc = &matrix_.jcn_loc.at(0);
    id_.a_loc   = &matrix_.a_loc.at(0);
  }

  if(i_am_root)
  {
    id_.rhs  = &rhs_global.at(0);
    id_.nrhs = 1;
    id_.lrhs = matrix_.n;
    id_.n    = matrix_.n;
  }

  if(!analyzed_and_factorized)
  {
    kick_off_job(std::string("analyze"));
    kick_off_job(std::string("factorize"));
    analyzed_and_factorized = true;
  }
  kick_off_job(std::string("solve"));

  //3) distribute solution among processes
  loc_dof_shift = 0;
  loc_master_shift = 0;
  glob_dof_shift = 0;

  for(size_t index = 0; index < comms_.size(); ++index) //loop over blocks
  {
    //receive all master values from root and store in master_values
    int n_loc_masters_block = comms_.at(index)->GetN_Master();
    int n_glob_dofs_block;
    MPI_Allreduce(&n_loc_masters_block,&n_glob_dofs_block,
                  1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    double* global_rhs_dummy = nullptr;
    if(i_am_root)
      global_rhs_dummy = &rhs_global.at(glob_dof_shift);

    //BUGFIX: make sure that master_values.at(loc_master_shift) does not throw,
    //even if the number of masters on this rank in this block was 0
    if(n_loc_masters_block == 0)
    {
      master_values.push_back(0);
    }

    scatter_vector(
        global_rhs_dummy,                                         //send
        &master_values.at(loc_master_shift), n_loc_masters_block, //receive
        root_rank);                                               //control

    //BUGFIX: make sure that master_values.at(loc_master_shift) does not throw,
    //even if the number of masters on this rank in this block was 0
    if(n_loc_masters_block == 0)
    {
      master_values.pop_back();
    }

    //write master values into local solution vector
    const int* masters = comms_.at(index)->GetMaster(); //TODO this should be a vector (in ParFECommunicator)!
    int n_local_dofs_block = comms_.at(index)->GetNDof();
    int i_master = 0;
    for(int i_dof = 0; i_dof< n_local_dofs_block; ++i_dof)
    {//write solution values for master dofs on this rank into solution vector
      if(masters[i_dof] == mpi_rank)
      {
        solution.at(loc_dof_shift + i_dof)
            = master_values.at(loc_master_shift + i_master);
        ++i_master;
      }
    }

    // Update all non-master values in solution vector. Big fire!
    comms_.at(index)->consistency_update(&solution.at(loc_dof_shift),3);

    // count up the block shifts
    int n_loc_dofs_block = comms_.at(index)->GetNDof();
    loc_dof_shift += n_loc_dofs_block;
    loc_master_shift += n_loc_masters_block;
    glob_dof_shift += n_glob_dofs_block;
  }
}

void MumpsWrapper::write_matrix_distributed(const std::string& filename) const
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


// Special member function.
MumpsWrapper::~MumpsWrapper()
{
  id_.job=JOB_END;
  dmumps_c(&id_);
}

// Private functions.
void MumpsWrapper::check_input_solve(
    const BlockVector& rhs, const BlockVector& solution,
    int& n_masters_local_comms, int& n_dofs_global_comms)
{
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  //0) input checks
   n_masters_local_comms = 0;
   n_dofs_global_comms = 0;
   for (size_t i =0; i < comms_.size(); ++i)
   {
     if(comms_.at(i)->get_n_dim() != 1 )
       throw std::runtime_error("MumpsWrapper can only be used with "
           "ParFECommunicators of dimension 1");

     size_t local_n_dofs = comms_.at(i)->GetNDof();
     size_t block_length_rhs = rhs.length(i);
     size_t block_length_sol = solution.length(i);

     if(local_n_dofs != block_length_rhs)
     {
       ErrThrow("Length of rhs block ", i, " does not"
           " fit n of local dofs in given communicator.");
     }
     if(local_n_dofs != block_length_sol)
     {
       ErrThrow("Length of sol block ", i, " does not fit n "
           "of local dofs in given communicator.");
     }
     // add up all masters
     n_masters_local_comms += comms_.at(i)->GetN_Master();
   }
   MPI_Allreduce(&n_masters_local_comms, &n_dofs_global_comms,
                 1, MPI_INT,MPI_SUM, MPI_COMM_WORLD);
   if(n_dofs_global_comms != (int) matrix_.n)
   {
     ErrThrow("Total number of masters of given communicators "
         "not equal global stored matrix order.");
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
      Output::print("MUMPS job ", job ," failed INFOG(1) = ",  id_.INFOG(1));
      Output::print("MUMPS job ", job ," failed INFOG(2) = ",  id_.INFOG(2));
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
  id_.ICNTL(21) = 0; // centralized solution TODO try to build in id_.ICNTL(21) = 1 - distributed solution

  // parameters collected from old MumpsSolver.C
  id_.ICNTL(1) = 6;    // standard outstream for errors
  id_.ICNTL(2) = 6;    // standard warnings output
  id_.ICNTL(3) = 0;    // global info output - surpressed ("max trans not allowed...")
  id_.ICNTL(4) = 1;    // verbosity level
  // NOTE: Clemens suggestion for this parameter to be 75 original was 20
  id_.ICNTL(14) = 75;   //estimated working space increase (%)

  //the following block is for choice of ordering tools in analysis phase
  // FIXME parellel ordering with parmetis is segfaulting!
  id_.ICNTL(28)=1;//request seq(1)/par(2) ordering in analysis phase
  //request mumps to independently choose a sequential ordering tool if id_.ICNTL(28)=1
  // in the current setup METIS and PORD are available plus some built-in orderings
  id_.ICNTL(7)=7;
  // if id_.ICNTL(28)=2 ordering is done in parallel, if a tool is available
  // FIXME in the current setup of the libraries, parmetis is segfaulting when
  // called - that's why id_.ICNTL(28) is set to 1 (sequential ordering) so far
  id_.ICNTL(29)=2;//request parmetis to do the ordering if id.ICNTL(28)=2

  // This parameter forces MUMPS to repeat the solution 2 times
  // using a simple iterative scheme. This (usually) makes the solution
  // more exact by several orders of magnitude.
  id_.ICNTL(10)=-2;
}

void MumpsWrapper::store_in_distributed_coordinate_form(
    const BlockMatrix& bmatrix,
    std::vector<double> pres0,
    std::vector<std::vector<int>> loc_to_seq
)
{
  //check what kind of block matrix we deal with
  bool is_block_fe_matrix = true;
  try
  {
    dynamic_cast<const BlockFEMatrix&>(bmatrix);
  }
  catch(const std::bad_cast& e)
  {//cast did not work
    is_block_fe_matrix = false;
  }

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
    for(auto comm : comms_)
    {
      n_masters_local += comm->GetN_Master();
    }
    MPI_Allreduce(&n_masters_local, &n_dofs_global,
                  1,MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }

  // 0.1)
  //determine an array which holds index-shifts
  size_t nComms = comms_.size();
  std::vector<size_t> shifts(nComms,0);
  //for each space, find out how many dofs there are in total
  //(over all processors) - adding up those will give the shift
  for(size_t index =0; index< nComms - 1;++index)
  {
    int n_masters_local = comms_.at(index)->GetN_Master();
    //sum up all nLocalMasters over all processes and store
    MPI_Allreduce(&n_masters_local, &shifts.at(index + 1),
                  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //add up the shift with the previous one
    shifts.at(index + 1) += shifts.at(index);
  }

  // 0.2)
  // peel the master arrays from the Communicators -make copies!
  std::vector<std::vector<int>> dof_masters(nComms);
  for (size_t index = 0; index<nComms; ++index)
  {
    const int* master_array = comms_.at(index)->GetMaster();
    size_t master_array_length = comms_.at(index)->GetNDof();
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
    const int* l2g = comms_.at(index)->Get_Local2Global();
    size_t array_length = comms_.at(index)->GetNDof();
    local2globals.at(index).resize(array_length);
    std::copy(l2g,
              l2g + array_length, //pointer arithmetic!
              local2globals.at(index).begin());
  }

  // 1) setting up the matrix
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
  std::shared_ptr<const TMatrix> block = bmatrix.get_block(0, 0, transp);

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
      bool is_bbt_type = ( block->get_sparse_type() == SparsityType::B_TIMES_BT );
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
            bool is_in_active_row = true; //true for standard block matrix
            if(is_block_fe_matrix)
              is_in_active_row  = (row < static_cast<const BlockFEMatrix&>(bmatrix).get_n_row_actives(cell_row) );

            // If this is BB^T case, the matrix is assumed to be
            // in additive storage and all entries are made count.
            if(is_bbt_type)
              is_in_master_row = true;

            if (is_nonzero && is_in_master_row && is_in_active_row)
            {//put entry into the new mumps matrix
              if(loc_to_seq.empty())
              {
                matrix_.irn_loc.push_back(row_shift + row_l2g.at(row) + 1);
                matrix_.jcn_loc.push_back(col_shift + col_l2g.at(col) + 1);
              }
              else
              {
                matrix_.irn_loc.push_back(row_shift + loc_to_seq.at(cell_row).at(row) + 1);
                matrix_.jcn_loc.push_back(col_shift + loc_to_seq.at(cell_col).at(col) + 1);
              }
              matrix_.a_loc.push_back(entry);
              ++matrix_.nz_loc;
            }
          }
        }
      }//end non-transposed case
      else if (transp) //block is stored transposed
      {
        if ( is_block_fe_matrix)
        {
          // static cast is okay, because we are sure to deal with a BlockFEMatrix
          if(std::static_pointer_cast<const FEMatrix>(block)->GetActiveBound() != block->GetN_Rows()) //check for security
          {
            ErrThrow("This block (",cell_row,",",cell_col,") has test space "
                     "non-actives, it should never have been stored in transposed state!");
          }
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
            bool is_in_active_row = true; //true for standard block matrix
            if(is_block_fe_matrix)
               is_in_active_row  = (col < static_cast<const BlockFEMatrix&>(bmatrix).get_n_row_actives(cell_row) );

            if (is_nonzero && is_in_master_row && is_in_active_row)
            {//fill entry in the sparsity structure
              if(loc_to_seq.empty())
              {
                matrix_.irn_loc.push_back(row_shift + row_l2g.at(col) + 1);
                matrix_.jcn_loc.push_back(col_shift + col_l2g.at(row) + 1);
              }
              else
              {
                matrix_.irn_loc.push_back(row_shift + loc_to_seq.at(cell_row).at(col) + 1);
                matrix_.jcn_loc.push_back(col_shift + loc_to_seq.at(cell_col).at(row) + 1);
              }
              matrix_.a_loc.push_back(entry);
              ++matrix_.nz_loc;
            }
          }
        }
      }
      if (is_block_fe_matrix)
      {
        //on diagonal blocks, set the ones on diagonals of non-active rows
        // (same code both cases, but loop does only make iterations for
        // non-transp case, because transposed storage of blocks with testspace
        // non-actives is not allowed in BlockFEMatrix)
        if(cell_row == cell_col)
        {
          for(int row = std::static_pointer_cast<const FEMatrix>(block)->GetActiveBound(); row < block->GetN_Rows(); ++row)
          {
            if(masters.at(row)==mpi_rank)
            {//only add the entry if the current process is master of its row
              //put a 1 on the diagonal
              if(loc_to_seq.empty())
              {
                matrix_.irn_loc.push_back(row_shift + row_l2g.at(row) + 1);
                matrix_.jcn_loc.push_back(row_shift + col_l2g.at(row) + 1);
              }
              else
              {
                matrix_.irn_loc.push_back(row_shift + loc_to_seq.at(cell_row).at(row) + 1);
                matrix_.jcn_loc.push_back(row_shift + loc_to_seq.at(cell_col).at(row) + 1);
              }
              matrix_.a_loc.push_back(1);
              ++matrix_.nz_loc;
            }
          }
        }//end treating dirichlet rows
      }
    }//end loop over columns
  }//end loop over rows

  if (is_block_fe_matrix)
  {
    // if the matrix comes from an enclosed flow problem,
    // an internal pressure row correction is necessary
    if (static_cast<const BlockFEMatrix&>(bmatrix).pressure_projection_enabled())
    {
      Output::root_info<5>("Pressure Projection", "MumpsWrapper applying pressure row correction.");
      pressure_row_correction(pres0);
    }
  }
}

void MumpsWrapper::pressure_row_correction(
    std::vector<double> pres0)
{
  if(!pres0.empty())
  { //
    const TFESpace3D* pres_space = this->comms_.at(3)->get_fe_space();
    double p0_x = pres0.at(0);
    double p0_y = pres0.at(1);
    double p0_z = pres0.at(2);
    int pressure_dof_to_correct = -1;
    for(int d = 0 ; d<pres_space->GetN_DegreesOfFreedom() ; ++d)
    {
      double x,y,z;
      double tol = 1e-6;
      pres_space->GetDOFPosition(d,x,y,z);
      if(    std::abs(p0_x - x) < tol
          && std::abs(p0_y - y) < tol
          && std::abs(p0_z - z) < tol)
      {
        pressure_dof_to_correct = d;
        break;
      }
    }

    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n_velo_dofs_local = comms_.at(0)->GetN_Master(); //master dofs
    int sendbuf[1]={n_velo_dofs_local};
    int recvbuf[1]={0};
    //gather total number of velo dofs in root 0
    MPI_Allreduce(sendbuf, recvbuf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    bool i_am_master = false;
    if(pressure_dof_to_correct != -1)
    {
        i_am_master = (comms_.at(3)->GetMaster()[pressure_dof_to_correct] == rank);
    }

    if(i_am_master)
    {//correct that entire row on "my" process
      Output::print("I AM YOUR MASTER!");
      int n_velo_dofs_global = 3*recvbuf[0]; //we're in 3D, thus multiply by 3
      //go through the matrix and remove all entries
      bool found_diagonal = false;
      for(size_t i = 0; i<matrix_.nz_loc; ++i)
      {
        if(matrix_.irn_loc.at(i) == n_velo_dofs_global + pressure_dof_to_correct + 1)
        {
          if(matrix_.jcn_loc.at(i) != n_velo_dofs_global + pressure_dof_to_correct + 1)
          {//off-diagonal entry - set to zero!
            matrix_.a_loc.at(i) = 0;
          }
          else
          {//diagonal entry - put to one!
            matrix_.a_loc.at(i) = 1;
            found_diagonal = true;
          }
        }
      }
      if(!found_diagonal)
      {
        matrix_.irn_loc.push_back(n_velo_dofs_global + pressure_dof_to_correct + 1);
        matrix_.jcn_loc.push_back(n_velo_dofs_global + pressure_dof_to_correct + 1);
        matrix_.a_loc.push_back(1);
        matrix_.nz_loc = matrix_.nz_loc + 1;
      }
    }

  }
  else
  {
    int size, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int n_velo_dofs_local = comms_.at(0)->GetN_Master(); //master dofs
    int sendbuf[1]={n_velo_dofs_local};
    int recvbuf[1]={0};
    //gather total number of velo dofs in root 0
    MPI_Reduce(sendbuf, recvbuf, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    // due to the used global ordering it is process 0 that holds the
    // globally first pressure row. that row should to be set to unit vector
    if(my_rank == 0)
    {
      int n_velo_dofs_global = 3*recvbuf[0]; //we're in 3D, thus multiply by 3
      int first_p_rw = n_velo_dofs_global;
      //go through the matrix and remove all entries
      for(size_t i = 0; i<matrix_.nz_loc; ++i)
      {
        if(matrix_.irn_loc.at(i)==first_p_rw)
        {
          if(matrix_.jcn_loc.at(i) != first_p_rw)
          {//off-diagonal entry - set to zero!
            matrix_.a_loc.at(i) = 0;
          }
          else
          {//diagonal entry - put to one!
            matrix_.a_loc.at(i) = 1;
          }
        }
      }
    }
  }

}

void MumpsWrapper::gather_vector(
    double* GlobalArray, double *LocalArray, int LocalSize, int root) const
{
  MPI_Comm Comm = MPI_COMM_WORLD;
  int size;
  MPI_Comm_size(Comm, &size);

  int *displ         = new int[size];
  int *N_ElementsAll = new int[size];
  //determine how many elements to receive per process - N_ElementsAll
  MPI_Allgather(&LocalSize, 1, MPI_INT, N_ElementsAll, 1, MPI_INT, Comm);

  displ[0] = 0;
  for(int i=1;i<size;i++)
    displ[i] = displ[i-1] + N_ElementsAll[i-1];

  MPI_Gatherv(LocalArray, LocalSize, MPI_DOUBLE, //send
              GlobalArray, N_ElementsAll, displ, MPI_DOUBLE, //receive
              root, Comm); //control

  delete [] displ;
  delete [] N_ElementsAll;
}

void MumpsWrapper::scatter_vector(
    double *GlobalArray, double *LocalArray, int LocalSize, int root) const
{
  MPI_Comm Comm = MPI_COMM_WORLD;
  int size;
  MPI_Comm_size(Comm, &size);

  int *displ         = new int[size];
  int *N_ElementsAll = new int[size];
  MPI_Allgather(&LocalSize, 1, MPI_INT, N_ElementsAll, 1, MPI_INT, Comm);

  displ[0] = 0;
  for(int i=1;i<size;i++)
    displ[i] = displ[i-1] + N_ElementsAll[i-1];

  MPI_Scatterv(GlobalArray, N_ElementsAll, displ, MPI_DOUBLE, //send
               LocalArray, LocalSize, MPI_DOUBLE,             //receive
               root, Comm);                                   //control

  delete [] displ;
  delete [] N_ElementsAll;
}

void MumpsWrapper::check_input_matrix(const BlockMatrix& bmatrix)
{
  // TODO This method should also react to the parameter
  // INTERNAL_FULL_MATRIX_STRUCTURE, which is important for algebraic flux
  // correction schemes of fem-fct and fem-tvd type. Because this is not
  // implemented yet, the MUMPS solver is currently not applicable for problems
  // that were stabilized with AFC.
  if(TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE)
    ErrThrow("MumpsWrapper cannot deal with INTERNAL_FULL_MATRIX_STRUCTURE "
        "(a parameter enabled, e.g., when FEM-FCT is used) yet. Use an iterative "
        "solver instead.")

  // Check that, if there is a B*B^T-matrix contained, it
  // is the only matrix and it is non-transposed.
  for(size_t row = 0 ; row<bmatrix.get_n_cell_rows(); ++row )
  {
    for(size_t col =0; col < bmatrix.get_n_cell_rows(); ++col )
    {
      bool transp;
      if(bmatrix.get_block(row,col,transp)->get_sparse_type()==SparsityType::B_TIMES_BT)
      {
        if(bmatrix.get_n_cell_rows() != 1 || bmatrix.get_n_cell_rows() != 1 || transp)
        {
          //currently I do not feel like thinking about this case
          ErrThrow("SparsityType::B_TIMES_BT matrix found in "
              "non-1x1 BlockMatrix (or in transposed state).");
        }
        else
        {
          Output::root_info("MumpsWrapper", "SparsityType::B_TIMES_BT matrix found.");
        }
      }
    }
  }


}

#endif
