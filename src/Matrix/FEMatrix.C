#include <FEMatrix.h>
#include <MooNMD_Io.h>
#include <algorithm>


FEMatrix::FEMatrix(std::shared_ptr<const TFESpace1D> space)
: FEMatrix(space, std::make_shared<TStructure>(space))
{
  
}

FEMatrix::FEMatrix(std::shared_ptr<const TFESpace2D> space)
: FEMatrix(space, std::make_shared<TStructure>(space))
{
  
}

#ifdef __3D__
FEMatrix::FEMatrix(std::shared_ptr<const TFESpace3D> space)
: FEMatrix(space, std::make_shared<TStructure>(space))
{
  
}
#endif // 3D

FEMatrix::FEMatrix(std::shared_ptr<const TFESpace2D> testspace,
                   std::shared_ptr<const TFESpace2D> ansatzspace, bool is_empty)
 : TMatrix(std::make_shared<TStructure>(testspace, ansatzspace, is_empty)),
   AnsatzSpace1D(nullptr), AnsatzSpace2D(ansatzspace), AnsatzSpace3D(nullptr),
   TestSpace1D(nullptr), TestSpace2D(testspace), TestSpace3D(nullptr)
{
  
}

#ifdef __3D__
FEMatrix::FEMatrix(std::shared_ptr<const TFESpace3D> testspace,
                   std::shared_ptr<const TFESpace3D> ansatzspace, bool is_empty)
: TMatrix(std::make_shared<TStructure>(testspace, ansatzspace, is_empty)),
  AnsatzSpace1D(nullptr), AnsatzSpace2D(nullptr), AnsatzSpace3D(ansatzspace),
  TestSpace1D(nullptr), TestSpace2D(nullptr), TestSpace3D(testspace)
{
 
}
#endif // 3D

FEMatrix::FEMatrix(std::shared_ptr<const TFESpace1D> space,
                   std::shared_ptr<TStructure> structure)
: TMatrix(structure),
  AnsatzSpace1D(space), AnsatzSpace2D(nullptr), AnsatzSpace3D(nullptr),
  TestSpace1D(space), TestSpace2D(nullptr), TestSpace3D(nullptr)
{
  if(!structure->isSquare())
  {
    ErrThrow("The structure must be square for this FEMatrix constructor");
  }
  if(space->GetN_DegreesOfFreedom() != structure->GetN_Rows())
  {
    ErrThrow("The given matrix and structure to not properly match");
  }
}

FEMatrix::FEMatrix(std::shared_ptr<const TFESpace2D> space,
                   std::shared_ptr<TStructure> structure)
: TMatrix(structure),
  AnsatzSpace1D(nullptr), AnsatzSpace2D(space), AnsatzSpace3D(nullptr),
  TestSpace1D(nullptr), TestSpace2D(space), TestSpace3D(nullptr)
{
  if(!structure->isSquare())
  {
    ErrThrow("The structure must be square for this FEMatrix constructor");
  }
  if(space->GetN_DegreesOfFreedom() != structure->GetN_Rows())
  {
    ErrThrow("The given matrix and structure to not properly match");
  }
}

#ifdef __3D__
FEMatrix::FEMatrix(std::shared_ptr<const TFESpace3D> space,
                   std::shared_ptr<TStructure> structure)
: TMatrix(structure),
  AnsatzSpace1D(nullptr), AnsatzSpace2D(nullptr), AnsatzSpace3D(space),
  TestSpace1D(nullptr), TestSpace2D(nullptr), TestSpace3D(space)
{
  if(!structure->isSquare())
  {
    ErrThrow("The structure must be square for this FEMatrix constructor");
  }
  if(space->GetN_DegreesOfFreedom() != structure->GetN_Rows())
  {
    ErrThrow("The given matrix and structure to not properly match");
  }
}
#endif // 3D

void FEMatrix::resetActive()
{
  // numer of entries in active rows
  int nActive = this->structure->getNActiveEntries();
  std::fill(this->entries.begin(), this->entries.begin()+nActive, 0.0);
}

void FEMatrix::resetNonActive()
{
  // numer of entries in active rows
  int nActive = this->structure->getNActiveEntries();
  std::fill(this->entries.begin() + nActive, this->entries.end(), 0.0);
}

void FEMatrix::scaleActive(double factor)
{
  if(factor == 1.0)
    return; // no scaling
  if(factor == 0.0)
    this->resetActive();
  
  // numer of entries in active rows
  int nActive = this->structure->getNActiveEntries();
  std::for_each(this->entries.begin(), this->entries.begin() + nActive,
                 [factor](double & a){ a = a*factor; } );
}

void FEMatrix::scale_non_active_diagonals(double factor)
{
  const int n_rows = this->GetN_Rows();
  const int active_bound = this->GetTestSpace()->GetActiveBound();
  const int * rowPtr = this->GetRowPtr();
  const int * colIndex = this->GetKCol();
  for(int i = active_bound; i < n_rows; ++i)
  {
    for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j)
    {
      if(colIndex[j] == i)
      {
        this->entries[j] *= factor;
      }
    }
  }
}


void FEMatrix::addActive(const FEMatrix& m, double factor)
{
  if(this->GetStructure() != m.GetStructure()) // compare objects
  {
    ErrThrow("FEMatrix::add : the two matrices do not match.");
  }
  
  // numer of entries in active rows
  int nActive = this->structure->getNActiveEntries();
  std::transform(this->entries.begin(), this->entries.begin() + nActive,
                 m.entries.begin(), this->entries.begin(), 
                 [factor](const double & a, const double & b)
                 { return a + factor * b; } );
}

void FEMatrix::multiplyActive(const double* x, double* y, double factor) const
{
  int nActive= this->GetActiveBound();
  const int * rowPtr = this->GetRowPtr();
  const int * colIndex = this->GetKCol();
  
  for(int i=0; i<nActive; ++i)
  {
    double val=0.;
    for(int j=rowPtr[i]; j<rowPtr[i+1]; ++j)
      val += this->entries[j]*x[colIndex[j]];
    y[i] += factor*val;
  }  
}

void FEMatrix::multiplyTransposedActive(const double *x, double *y, double factor) const
{
  //be careful! we have to rely on y's actives being as many as this ansatz spaces
  //FIXME this can be sped up of course, but for the moment do it like this
  //assume that y is as long as this has columns
  int n_actives = this->GetAnsatzSpace()->GetActiveBound();
  std::vector<double> y_non_actives(this->GetAnsatzSpace()->GetN_Dirichlet());
  for (int i= n_actives; i<this->GetN_Columns() ; ++i )
  {//store non-actives
    y_non_actives[i-n_actives]=y[i];
  }
  this->TMatrix::transpose_multiply(x,y,factor);
  for (int i= n_actives; i<this->GetN_Columns() ; ++i )
  {//put back non-actives
    y[i] = y_non_actives[i-n_actives];
  }

}


int FEMatrix::GetActiveBound() const
{
  return structure->GetActiveBound();
}

std::shared_ptr<const TFESpace1D> FEMatrix::GetTestSpace1D() const
{
  return TestSpace1D;
}

std::shared_ptr<const TFESpace1D> FEMatrix::GetAnsatzSpace1D() const
{
  return AnsatzSpace1D;
}

std::shared_ptr<const TFESpace2D> FEMatrix::GetTestSpace2D() const
{
  return TestSpace2D;
}

std::shared_ptr<const TFESpace2D> FEMatrix::GetAnsatzSpace2D() const
{
  return AnsatzSpace2D;
}

#ifdef __3D__
std::shared_ptr<const TFESpace3D> FEMatrix::GetTestSpace3D() const
{
  return TestSpace3D;
}

std::shared_ptr<const TFESpace3D> FEMatrix::GetAnsatzSpace3D() const
{
  return AnsatzSpace3D;
}
#endif // 3D

std::shared_ptr<const TFESpace> FEMatrix::GetTestSpace() const
{
  if(TestSpace1D)
  {
    return TestSpace1D;
  }
  else if(TestSpace2D)
  {
    return TestSpace2D;
  }
  else
  {
    return TestSpace3D;
  }
}

std::shared_ptr<const TFESpace> FEMatrix::GetAnsatzSpace() const
{
  if(AnsatzSpace1D)
  {
    return AnsatzSpace1D;
  }
  else if(AnsatzSpace2D)
  {
    return AnsatzSpace2D;
  }
  else
  {
    return AnsatzSpace3D;
  }
}

std::shared_ptr<const TFESpace1D> FEMatrix::GetFESpace1D() const
{
  if(this->structure->isSquare())
    return TestSpace1D;
  else
    ErrThrow("accessing FESpace for non-square matrix, but which one?");
}

std::shared_ptr<const TFESpace2D> FEMatrix::GetFESpace2D() const
{
  if(this->structure->isSquare())
    return TestSpace2D;
  else
    ErrThrow("accessing FESpace for non-square matrix, but which one?");
}

#ifdef __3D__
std::shared_ptr<const TFESpace3D> FEMatrix::GetFESpace3D() const
{
  if(this->structure->isSquare())
    return TestSpace3D;
  else
    ErrThrow("accessing FESpace for non-square matrix, but which one?");
}
#endif // 3D

#ifdef _MPI

#include "mpi.h"
#include <ParFECommunicator3D.h>

/* ************************************************************************ */
// The following are two functions which were used for debug reasons. They do
// not belong to the class FEMatrix.
// Yet they should not be lost entirely, so they were 'parked' here.

/**
 * This method can be used to check a matrix A in distributed storage.
 * A certain process ('sending_ps') will extract one master column after the other
 * from A, and ask all other processes for a consistency update of that vector.
 *
 * After that it will compare the updated column with the original one and inform
 * the programmer via Output::print on any differences.
 * It is then in the responsibility of the programmer to interpret the output.
 *
 * @param[in] A The matrix to check.
 * @param[in] sending_ps The process which checks its master columns.
 * @param[in] cons_level The consistency level to be updated to.
 */
void check_column_consistency(const FEMatrix& A, int sending_ps, int cons_level)
{

  const TParFECommunicator3D& comm = A.GetFESpace3D()->get_communicator();
  const int* masters = comm.GetMaster();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nDof = A.GetN_Rows();

  const char* markers = comm.get_dof_markers();

  //SENDER process executes the following code.
  if(rank == sending_ps)
  {
    Output::info("CHECK", "Checking master columns of process ", sending_ps, " on consistency level ", cons_level);
    int n_vecs_send = comm.GetN_Master();
    MPI_Bcast(&n_vecs_send, 1, MPI_INT, sending_ps, MPI_COMM_WORLD);
    for(int s = 0; s < nDof; ++s)
    {
      if(masters[s] == rank) //master col found, ping it
      {
        Output::suppressAll();
        comm.dof_ping(sending_ps, s);
        Output::setVerbosity(1);
        int dummy_sb = -1;
        int dof_remote = -1;
        MPI_Allreduce(&dummy_sb, &dof_remote, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if(dof_remote != -1)
        {//ping was received on at least one other ps, now set up the matrix row
          std::vector<double> col_master = A.get_matrix_column(s);
          std::vector<double> col_master_cpy = col_master;
          comm.consistency_update(col_master.data(), cons_level); //CONSIST UPDATE
          //now check the differences between updated and original column
          for(size_t i = 0 ; i < col_master.size(); ++i)
          {
            if(col_master[i] != col_master_cpy[i])
            {
              char type_col = markers[s];
              char type_row = markers[i];
              Output::print("Local entry (", i, "[", type_row ,"] , ", s, "[", type_col,"]) was changed by an update.");
            }

          }

        }
      }
    }
    Output::print("Test on process ", rank, " complete.");
  }
  //RECEIVER processes execute the following code.
  else
  {
    int n_vecs_recv;
    MPI_Bcast(&n_vecs_recv, 1, MPI_INT, sending_ps, MPI_COMM_WORLD);
    for(int r = 0; r< n_vecs_recv; ++r)
    {
      int dof_local = comm.dof_ping(sending_ps, -1);
      int dof_remote = -1;
      MPI_Allreduce(&dof_local, &dof_remote, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if(dof_remote != -1)
      {
        std::vector<double> col_slave = A.get_matrix_column(dof_local);
        comm.consistency_update(col_slave.data(), cons_level); //CONSIST UPDATE
      }
    }
  }
}

// Restores consistency of master columns.
void update_column_consistency(FEMatrix& B, int sending_ps, int cons_level)
{
  const TParFECommunicator3D& comm = B.GetFESpace3D()->get_communicator();
  const int* masters = comm.GetMaster();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nDof = B.GetN_Rows();

  const char* markers = comm.get_dof_markers();

  //SENDER process executes the following code.
  if(rank == sending_ps)
  {
    int n_vecs_send = comm.get_n_interface_master();
    MPI_Bcast(&n_vecs_send, 1, MPI_INT, sending_ps, MPI_COMM_WORLD);
    Output::info("MPI UPDATE", n_vecs_send, " interface master columns receive update.");
    for(int s = 0; s < nDof; ++s)
    {
      if(masters[s] == rank && markers[s] == 'm') //interface master col found, ping it
      {
        Output::suppressAll();
        comm.dof_ping(sending_ps, s);
        Output::setVerbosity(1);
        int dummy_sb = -1;
        int dof_remote = -1;
        MPI_Allreduce(&dummy_sb, &dof_remote, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if(dof_remote != -1)
        {//ping was received on at least one other ps, now set up the matrix row
          std::vector<double> col_master = B.get_matrix_column(s);
          comm.consistency_update(col_master.data(), cons_level); //CONSIST UPDATE
          B.set_matrix_column(s,col_master);

        }
      }
    }
  }
  //RECEIVER processes execute the following code.
  else
  {
    int n_vecs_recv;
    MPI_Bcast(&n_vecs_recv, 1, MPI_INT, sending_ps, MPI_COMM_WORLD);
    for(int r = 0; r< n_vecs_recv; ++r)
    {
      int dof_local = comm.dof_ping(sending_ps, -1);
      int dof_remote = -1;
      MPI_Allreduce(&dof_local, &dof_remote, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if(dof_remote != -1)
      {
        std::vector<double> col_slave = B.get_matrix_column(dof_local);
        comm.consistency_update(col_slave.data(), cons_level); //CONSIST UPDATE
      }
    }
  }
}
#endif //_MPI
