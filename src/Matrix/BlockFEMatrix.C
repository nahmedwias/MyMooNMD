#include <BlockFEMatrix.h>
#include <MooNMD_Io.h>
#include <LinAlg.h>
#include <limits>

/** ************************************************************************* */
BlockFEMatrix::BlockFEMatrix(const Problem_type type, unsigned int space_dimension,
                bool mass_matrix)
: BlockMatrix(type, space_dimension,  mass_matrix),
  actives(block_pattern->n_blocks(), std::numeric_limits<unsigned int>::max())
{
  ;
}


/** ************************************************************************* */
void BlockFEMatrix::add_scaled_active(const BlockMatrix& A, double factor)
{
  unsigned int n_blocks = A.n_blocks();
  if(this->n_blocks() != n_blocks)
  {
    ErrThrow("BlockFEMatrix::add_scaled_active : the two BlockMatrix objects do ",
             "not have the same number of blocks.");
  }

  for(unsigned int b = 0; b < n_blocks; b++)
  {
    if(this->actives[b] <= (unsigned int) this->blocks[b]->GetN_Entries())
    {
      // note: this could be a method of TMatrix as well, if TMatrix knew its
      // actives
      Daxpy(this->actives[b], factor, A.block(b)->GetEntries(),
            this->blocks[b]->GetEntries());
    }
    else
      ErrThrow("adding actives where the number of actives has not been set");
  }
}

/* ************************************************************************* */
void BlockFEMatrix::scale_active(double factor)
{
  for(unsigned int b = 0; b < this->n_blocks(); ++b)
  {
    if(this->actives[b] <= (unsigned int) this->blocks[b]->GetN_Entries())
    {
      // note: this could be a method of TMatrix as well, if TMatrix knew its
      // actives
      Dscal(this->actives[b], factor, this->blocks[b]->GetEntries());
    }
    else
      ErrThrow("scaling actives where the number of actives has not been set");
  }
}

/* ************************************************************************* */
size_t BlockFEMatrix::get_n_row_actives( size_t index ) const
{
  // check if there is as many block rows
  if(index > block_pattern->n_rows())
  {
    ErrThrow("Index out of bounds.");
  }
  return ((FEMatrix&)block( index , 0)).GetTestSpace()->GetActiveBound();
}

/* ************************************************************************* */
size_t BlockFEMatrix::get_n_column_actives( size_t index ) const
{
  // check if there is as many block rows
  if(index > block_pattern->n_cols())
  {
    ErrThrow("Index out of bounds.");
  }
  return ((FEMatrix&)block( 0, index)).GetTestSpace()->GetActiveBound();
}
