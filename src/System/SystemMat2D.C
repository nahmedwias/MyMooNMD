#include <SystemMat2D.h>
#include <string.h>

SystemMat2D::SystemMat2D(const SystemMat2D & mat)
 : fe_spaces(mat.fe_spaces), sq_matrices(mat.sq_matrices.size(), nullptr),
   rect_matrices(mat.rect_matrices.size(), nullptr), defect(mat.defect)
{
  // for each submatrix create a copy here
  //for(unsigned int i = 0; i < n_square_mat; ++i)
    ;//this->sq_matrices[i]->add_active(*A.get_square_matrix(i), factor);
  //for(unsigned int i = 0; i < n_rect_mat; ++i)
    ;//this->rect_matrices[i]->add_active(*A.get_rectangular_matrix(i), factor);
    ErrMsg("copy contructor for SystemMat2D not yet implemented, we need copy "
           << "constructors for TSquareMatrix2D and TMatrix2D first");
    throw("copy contructor for SystemMat2D not yet implemented");
}

SystemMat2D::~SystemMat2D()
{
  // to be implemented
}

void SystemMat2D::add_active(const SystemMat2D &A, double factor)
{
  unsigned int n_square_mat = this->sq_matrices.size();
  unsigned int n_rect_mat = this->rect_matrices.size();
  // check if matrix A can be added to this matrix at all
  if(this->fe_spaces.size() != A.fe_spaces.size())
  {
    ErrMsg("can't add matrices with different number of finite element spaces");
    throw("can't add matrices with different number of finite element spaces");
  }
  if(A.sq_matrices.size() != n_square_mat)
  {
    ErrMsg("can't add matrices with different number of square submatrices");
    throw("can't add matrices with different number of square submatrices");
  }
  if(A.rect_matrices.size() != n_rect_mat)
  {
    ErrMsg("can't add matrices with different number of rectangular matrices");
    throw("can't add matrices with different number of rectangular matrices");
  }
  
  if(factor == 0.0)
    // nothing needs to be done
    return;
  
  // add each submatrix
  for(unsigned int i = 0; i < n_square_mat; ++i)
    this->sq_matrices[i]->add_active(*A.get_square_matrix(i), factor);
  for(unsigned int i = 0; i < n_rect_mat; ++i)
    this->rect_matrices[i]->add_active(*A.get_rectangular_matrix(i), factor);
}

void SystemMat2D::add(const SystemMat2D &A, double factor)
{
  unsigned int n_square_mat = this->sq_matrices.size();
  unsigned int n_rect_mat = this->rect_matrices.size();
  // check if matrix A can be added to this matrix at all
  if(this->fe_spaces.size() != A.fe_spaces.size())
  {
    ErrMsg("can't add matrices with different number of finite element spaces");
    throw("can't add matrices with different number of finite element spaces");
  }
  if(A.sq_matrices.size() != n_square_mat)
  {
    ErrMsg("can't add matrices with different number of square submatrices");
    throw("can't add matrices with different number of square submatrices");
  }
  if(A.rect_matrices.size() != n_rect_mat)
  {
    ErrMsg("can't add matrices with different number of rectangular matrices");
    throw("can't add matrices with different number of rectangular matrices");
  }
  
  if(factor == 0.0)
    // nothing needs to be done
    return;
  
  // add each submatrix
  for(unsigned int i = 0; i < n_square_mat; ++i)
    this->sq_matrices[i]->add(*(TMatrix*)A.get_square_matrix(i), factor);
  for(unsigned int i = 0; i < n_rect_mat; ++i)
    this->rect_matrices[i]->add(*(TMatrix*)A.get_rectangular_matrix(i), factor);
}

void SystemMat2D::scale_active(double factor)
{
  //  scale each subblock
  for(TSquareMatrix2D* m : this->sq_matrices)
    m->scale_active(factor);
  for(TMatrix2D* m : this->rect_matrices)
    m->scale_active(factor);
}

void SystemMat2D::scale(double factor)
{
  //  scale each subblock
  for(TSquareMatrix2D* m : this->sq_matrices)
    m->scale(factor);
  for(TMatrix2D* m : this->rect_matrices)
    m->scale(factor);
}

void SystemMat2D::apply(const double *x, double *y, double factor) const
{
  unsigned int n_total_rows = this->n_total_rows();
  if(factor == 0.0)
  {
    // simply reset y, without actually multiplying the matrix with x
    memset(y, 0.0, n_total_rows*SizeOfDouble);
    return;
  }
  // this should be implemented in a derived class
  ErrMsg("no implementation of SystemMat2D::apply available");
  throw("no implementation of SystemMat2D::apply available");
}

void SystemMat2D::apply_scaled_add(const double *x, double *y, double factor)
  const
{
  if(factor == 0.0)
    // nothing needs to be done
    return;
  // this should be implemented in a derived class
  ErrMsg("no implementation of SystemMat2D::apply_scaled_add available");
  throw("no implementation of SystemMat2D::apply_scaled_add available");
}

unsigned int SystemMat2D::n_blocks() const
{
  return sq_matrices.size() + rect_matrices.size();
}

unsigned int SystemMat2D::n_rows() const
{
  // this should be implemented in a derived class
  // assuming a square block matrix, where every block is indeed stored
  return sqrt(sq_matrices.size() + rect_matrices.size());
}

unsigned int SystemMat2D::n_cols() const
{
  // this should be implemented in a derived class
  // assuming a square block matrix, where every block is indeed stored
  return sqrt(sq_matrices.size() + rect_matrices.size());
}

unsigned int SystemMat2D::n_total_rows() const
{
  // this should be implemented in a derived class
  if(this->n_blocks() == 1)
  {
    if(this->sq_matrices.size() == 1) // one square block
      return this->sq_matrices[0]->GetN_Rows();
    else if(this->rect_matrices.size() == 1) // one rectangular block
      return this->rect_matrices[0]->GetN_Rows();
  }
  ErrMsg("unable to determine the total number of rows in this SystemMat2D");
  throw("unable to determine the total number of rows in this SystemMat2D");
}

unsigned int SystemMat2D::n_total_cols() const
{
  // this should be implemented in a derived class
  if(this->n_blocks() == 1)
  {
    if(this->sq_matrices.size() == 1) // one square block
      return this->sq_matrices[0]->GetN_Columns();
    else if(this->rect_matrices.size() == 1) // one rectangular block
      return this->rect_matrices[0]->GetN_Columns();
  }
  ErrMsg("unable to determine the total number of columns in this SystemMat2D");
  throw("unable to determine the total number of columns in this SystemMat2D");
}
unsigned int SystemMat2D::n_total_entries() const
{
  // this should be implemented in a derived class
  unsigned int n_entries = 0;
  for(TSquareMatrix2D* m : this->sq_matrices)
    n_entries += m->GetN_Entries();
  for(TMatrix2D* m : this->rect_matrices)
    n_entries += m->GetN_Entries();
  return n_entries;
}
