#include <BlockMatrix3D.h>
#include <string.h>

BlockMatrix3D::~BlockMatrix3D()
{
  // to be implemented
}

void BlockMatrix3D::add_active(const BlockMatrix3D &A, double factor)
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
    ErrMsg("BlockMatrix3D::add not yet implemented for a factor != 1.0");
    throw("BlockMatrix3D::add not yet implemented for a factor != 1.0");
  }
  
  // add each submatrix
  for(unsigned int i = 0; i < n_square_mat; ++i)
    this->sq_matrices[i]->add_active(*A.get_square_matrix(i), factor);
  for(unsigned int i = 0; i < n_rect_mat; ++i)
    this->rect_matrices[i]->add_active(*A.get_rectangular_matrix(i), factor);
}

void BlockMatrix3D::add(const BlockMatrix3D &A, double factor)
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
    ErrMsg("BlockMatrix3D::add not yet implemented for a factor != 1.0");
    throw("BlockMatrix3D::add not yet implemented for a factor != 1.0");
  }
  
  // add each submatrix
  for(unsigned int i = 0; i < n_square_mat; ++i)
    this->sq_matrices[i]->add(*(TMatrix*)A.get_square_matrix(i), factor);
  for(unsigned int i = 0; i < n_rect_mat; ++i)
    this->rect_matrices[i]->add(*(TMatrix*)A.get_rectangular_matrix(i), factor);
}

void BlockMatrix3D::scale_active(double factor)
{
  //  scale each subblock
  for(TSquareMatrix3D* m : this->sq_matrices)
    m->scale_active(factor);
  for(TMatrix3D* m : this->rect_matrices)
    m->scale_active(factor);
}

void BlockMatrix3D::scale(double factor)
{
  //  scale each subblock
  for(TSquareMatrix3D* m : this->sq_matrices)
    m->scale(factor);
  for(TMatrix3D* m : this->rect_matrices)
    m->scale(factor);
}

void BlockMatrix3D::apply(const double *x, double *y, double factor) const
{
  unsigned int n_total_rows = this->n_total_rows();
  if(factor == 0.0)
  {
    // simply reset y, without actually multiplying the matrix with x
    memset(y, 0.0, n_total_rows*SizeOfDouble);
    return;
  }
  ErrMsg("no implementation of BlockMatrix3D::apply available");
  throw("no implementation of BlockMatrix3D::apply available");
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

void BlockMatrix3D::apply_scaled_add(const double *x, double *y, double factor)
  const
{
  if(factor == 0.0)
    // nothing needs to be done
    return;
  
  ErrMsg("no implementation of BlockMatrix3D::apply_scaled_add available");
  throw("no implementation of BlockMatrix3D::apply_scaled_add available");
}

unsigned int BlockMatrix3D::n_blocks() const
{
  return sq_matrices.size() + rect_matrices.size();
}

unsigned int BlockMatrix3D::n_rows() const
{
  // assuming a square block matrix, where every block is indeed stored
  return sqrt(sq_matrices.size() + rect_matrices.size());
}

unsigned int BlockMatrix3D::n_cols() const
{
  // assuming a square block matrix, where every block is indeed stored
  return sqrt(sq_matrices.size() + rect_matrices.size());
}

unsigned int BlockMatrix3D::n_total_rows() const
{
  if(this->n_blocks() == 1)
  {
    if(this->sq_matrices.size() == 1) // one square block
      return this->sq_matrices[0]->GetN_Rows();
    else if(this->rect_matrices.size() == 1) // one rectangular block
      return this->rect_matrices[0]->GetN_Rows();
  }
  ErrMsg("unable to determine the total number of rows in this BlockMatrix3D");
  throw("unable to determine the total number of rows in this BlockMatrix3D");
}

unsigned int BlockMatrix3D::n_total_cols() const
{
  if(this->n_blocks() == 1)
  {
    if(this->sq_matrices.size() == 1) // one square block
      return this->sq_matrices[0]->GetN_Columns();
    else if(this->rect_matrices.size() == 1) // one rectangular block
      return this->rect_matrices[0]->GetN_Columns();
  }
  ErrMsg("unable to determine the total number of columns in this BlockMatrix3D");
  throw("unable to determine the total number of columns in this BlockMatrix3D");
}
unsigned int BlockMatrix3D::n_total_entries() const
{
  unsigned int n_entries = 0;
  for(TSquareMatrix3D* m : this->sq_matrices)
    n_entries += m->GetN_Entries();
  for(TMatrix3D* m : this->rect_matrices)
    n_entries += m->GetN_Entries();
  return n_entries;
}
