#include <SystemMat2D.h>
#include <string.h>

SystemMat2D::~SystemMat2D()
{
  // to be implemented
}

void SystemMat2D::add_active(const SystemMat2D &A, double factor)
{
  
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
  
  if(factor != 1.0)
  {
    ErrMsg("SystemMat2D::add not yet implemented for a factor != 1.0");
    throw("SystemMat2D::add not yet implemented for a factor != 1.0");
  }
  
  // add each submatrix
  for(TSquareMatrix2D* m : this->sq_matrices)
    m->add(*(TMatrix*)m, factor);
  for(TMatrix2D* m : this->rect_matrices)
    m->add(*(TMatrix*)m, factor);
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
  ErrMsg("no implementation of SystemMat2D::apply available");
  throw("no implementation of SystemMat2D::apply available");
  /*
  unsigned int row_offset = 0;
  unsigned int n_rows = this->n_rows();
  unsigned int n_cols = this->n_cols();
  for(int i = 0; i < n_rows; i++)
  {
    int col_offset = 0;
    for(int j = 0; j < n_cols; j++)
    {
      TMatrix * current_block = _blocks[i * n_rows + j];
      current_block->multiply(x + col_offset, y + row_offset, 1.0);
      col_offset += current_block->GetN_Columns();
    }
    row_offset += _blocks[i * n_rows]->GetN_Rows();
  }
  */
}

void SystemMat2D::apply_scaled_add(const double *x, double *y, double factor)
  const
{
  if(factor == 0.0)
    // nothing needs to be done
    return;
  
  ErrMsg("no implementation of SystemMat2D::apply_scaled_add available");
  throw("no implementation of SystemMat2D::apply_scaled_add available");
}

unsigned int SystemMat2D::n_blocks() const
{
  return sq_matrices.size() + rect_matrices.size();
}

unsigned int SystemMat2D::n_rows() const
{
  // assuming a square block matrix, where every block is indeed stored
  return sqrt(sq_matrices.size() + rect_matrices.size());
}

unsigned int SystemMat2D::n_cols() const
{
  // assuming a square block matrix, where every block is indeed stored
  return sqrt(sq_matrices.size() + rect_matrices.size());
}

unsigned int SystemMat2D::n_total_rows() const
{
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
  unsigned int n_entries = 0;
  for(TSquareMatrix2D* m : this->sq_matrices)
    n_entries += m->GetN_Entries();
  for(TMatrix2D* m : this->rect_matrices)
    n_entries += m->GetN_Entries();
  return n_entries;
}
