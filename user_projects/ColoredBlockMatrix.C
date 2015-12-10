
#include <ColoredBlockMatrix.h>
#include <BlockVector.h>
#include <LinAlg.h>
#include <limits>

///* ************************************************************************* */
//ColoredBlockMatrix::ColoredBlockMatrix()
// : ColoredBlockMatrix(0, 0)
//{
//
//}

ColoredBlockMatrix::CellInfo::CellInfo()
: n_rows_{0}, n_columns_{0},
  //block_in_cell gets default initialised to null
  color_{std::numeric_limits<size_t>::max()}, // no color
  is_transposed_{false} //non-transposed

{
    //CB DEBUG
    Output::print("Debug: Default constructing CellInfo");
    //END DEBUG
}

///* ************************************************************************* */
ColoredBlockMatrix::CellInfo::CellInfo(size_t nRows, size_t nColumns)
: n_rows_{nRows}, n_columns_{nColumns}, //block_in_cell gets default initialised to null
  color_{std::numeric_limits<size_t>::max()}, // no colour
  is_transposed_{false} //non-transposed

{
    //CB DEBUG
    Output::print("Debug: Constructing CellInfo with n_rows_ ", n_rows_, " and n_columns_ ", n_columns_);
    //END DEBUG
}

/* ************************************************************************* */
ColoredBlockMatrix::ColoredBlockMatrix(std::vector<size_t> cellRowNumbers,
                         std::vector<size_t> cellColumnNumbers)
  : n_block_rows_{cellRowNumbers.size()},
    n_block_columns_{cellColumnNumbers.size()},
    cell_info_grid_(n_block_rows_, std::vector<CellInfo>(n_block_columns_)),
    n_colors_{0}
    // combinde_matrix gets default initialized to smart nullptr
{
  //traverse the cell info grid from left to right, top to bottom
  // and fill in newly constructed, correctly dimensioned zero matrices as blocks
 for(size_t i = 0; i < n_block_rows_ ; ++i )
 {
   //hold the number of rows each cell in this row will have
   size_t nRowsOfCell = cellRowNumbers[i];

   for(size_t j = 0; j < n_block_columns_ ; ++j )
   {
     //hold the number of columns each cell in this column will have
     size_t nColumnsOfCell = cellColumnNumbers[j];

     // construct the new cell info and the zero matrix it will hold
     CellInfo newInfo(nRowsOfCell, nColumnsOfCell);
     newInfo.block_in_cell_ = std::make_shared<TMatrix>(nRowsOfCell, nColumnsOfCell);

     // the next color is n_colors_ (because the colors 0 to n_colors_ -1 are already assigned)
     newInfo.color_ = n_colors_;

     //start as non-transposed
     newInfo.is_transposed_ = false;

     // put the new cell info to the correct place by copy assignment
     cell_info_grid_ [i][j] = newInfo;

     // an entirely new block has been added and colored,
     // so raise the number of colors by one
     ++n_colors_;
   }
 }
}

///* ************************************************************************* */
//ColoredBlockMatrix::ColoredBlockMatrix(unsigned int n_rows, unsigned int n_cols,
//                         std::vector<std::shared_ptr<TMatrix>> new_blocks)
// : ColoredBlockMatrix(n_rows, n_cols)
//{
//  if(new_blocks.size() < this->n_blocks())
//  {
//    ErrThrow("Creating a ColoredBlockMatrix with ", this->n_rows(), " rows and ",
//             this->n_cols(), " columns, but only ", new_blocks.size(),
//             " blocks given");
//  }
//  std::copy(new_blocks.begin(), new_blocks.end(), this->blocks.begin());
//
//  // check consistency
//  // for each row check if all blocks have the same number of rows
//  for(unsigned int row = 0; row < n_rows; ++row)
//  {
//    for(unsigned int col = 1; col < n_cols; ++col)
//    {
//      if(this->blocks[row*n_cols]->GetN_Rows()
//          != this->blocks[row*n_cols + col]->GetN_Rows())
//      {
//        ErrThrow("The blocks in row ", row,
//                 " do not have the same number of rows");
//      }
//    }
//  }
//  // for each column check if all blocks have the same number of columns
//  for(unsigned int col = 0; col < n_cols; ++col)
//  {
//    for(unsigned int row = 1; row < n_rows; ++row)
//    {
//      if(this->blocks[col]->GetN_Columns()
//          != this->blocks[row*n_cols + col]->GetN_Columns())
//      {
//        ErrThrow("The blocks in column ", col,
//                 " do not have the same number of columns");
//      }
//    }
//  }
//}
//
///* ************************************************************************* */
//ColoredBlockMatrix::ColoredBlockMatrix(const Problem_type type,
//                         unsigned int space_dimension, bool mass_matrix)
// : ColoredBlockMatrix(std::make_shared<const BlockPattern>(type, space_dimension,
//                                                    mass_matrix))
//{
//  // nothing more to do
//}
//
///* ************************************************************************* */
//ColoredBlockMatrix::ColoredBlockMatrix(std::shared_ptr<const BlockPattern> bp)
// : block_pattern(bp),
//   blocks(std::vector<std::shared_ptr<TMatrix>>(bp->n_blocks(), nullptr)),
//   combined_matrix(std::shared_ptr<TMatrix>())
//{
//  // matrices are not created here. You still have to do that
//}
//
///* ************************************************************************* */
//ColoredBlockMatrix::ColoredBlockMatrix(ColoredBlockMatrix& other)
// : block_pattern(other.block_pattern), blocks(other.blocks.size(), nullptr),
//   combined_matrix(other.combined_matrix)
//{
//  for(unsigned int b = 0; b < this->blocks.size(); ++b)
//  {
//    this->blocks[b] = other.blocks[b]; // set pointer
//  }
//}
//
///* ************************************************************************* */
//ColoredBlockMatrix::ColoredBlockMatrix(ColoredBlockMatrix&& other)
// : block_pattern(other.block_pattern), blocks(other.blocks.size(), nullptr),
//   combined_matrix(other.combined_matrix)
//{
//  for(unsigned int b = 0; b < this->blocks.size(); ++b)
//  {
//    this->blocks[b] = other.blocks[b]; // set pointer
//    // destructor on this sparse matrix is only called once from 'this', not
//    // from 'other'
//    other.blocks[b] = nullptr;
//  }
//}
//
///* ************************************************************************* */
//ColoredBlockMatrix::~ColoredBlockMatrix() noexcept
//{
//}
//
///* ************************************************************************* */
//void ColoredBlockMatrix::reset()
//{
//  for(unsigned int b = 0; b < this->n_blocks(); ++b)
//    this->blocks[b]->reset();
//}
//
///* ************************************************************************* */
//void ColoredBlockMatrix::add_scaled(const ColoredBlockMatrix& A, double factor)
//{
//  unsigned int n_blocks = A.n_blocks();
//  if(this->n_blocks() != n_blocks)
//  {
//    ErrThrow("ColoredBlockMatrix::add_scaled : the two ColoredBlockMatrix objects do ",
//             "not have the same number of blocks.");
//  }
//
//  for(unsigned int i = 0; i < n_blocks; i++)
//  {
//    this->block(i)->add_scaled(*A.block(i), factor);
//  }
//}
//
//
///** ************************************************************************* */
//void ColoredBlockMatrix::scale(double factor)
//{
//  for(unsigned int b = 0; b < this->n_blocks(); ++b)
//    this->blocks[b]->scale(factor);
//}
//
///* ************************************************************************* */
//void ColoredBlockMatrix::apply(const BlockVector & x, BlockVector & y) const
//{
//  unsigned int l = y.length();
//  if(l != this->n_total_cols() && l != 0)
//  {
//    ErrThrow("cannot multiply with matrix, dimension mismatch");
//  }
//  if(l == 0)
//  {
//    // BlockVector y is empty, set to to a suitable vector in the image of
//    // this ColoredBlockMatrix, true means y is in the image of this, rather
//    // than the pre-image
//    y.copy_structure(*this, true);
//    // all values of 'y' are set to 0
//  }
//  else
//    y.reset(); // set all values in 'y' to 0
//
//  this->apply_scaled_add(x, y, 1.0);
//}
//
///* ************************************************************************* */
//void ColoredBlockMatrix::apply_scaled_add(const BlockVector & x, BlockVector & y,
//                                         double a) const
//{
//  if(y.length() != this->n_total_rows())
//  {
//    ErrThrow("cannot multiply with matrix, dimension mismatch");
//  }
//  if(x.length() != this->n_total_cols())
//  {
//    ErrThrow("cannot multiply with matrix, dimension mismatch");
//  }
//
//  const double * xv = x.get_entries(); // array of values in x
//  double * yv = y.get_entries(); // array of values in y
//  unsigned int row_offset = 0;
//  // n_rows, n_cols are the number of block rows/columns
//  for(unsigned int i = 0, n_rows = this->n_rows(), n_cols = this->n_cols();
//      i < n_rows; ++i)
//  {
//    int col_offset = 0;
//    for(unsigned int j = 0; j < n_cols; j++)
//    {
//      auto current_block = this->block(i * n_rows + j);
//      current_block->multiply(xv + col_offset, yv + row_offset, a);
//      col_offset += current_block->GetN_Columns();
//    }
//    row_offset += this->block(i * n_rows)->GetN_Rows();
//  }
//}
//
///* ************************************************************************* */
//std::shared_ptr<TMatrix> ColoredBlockMatrix::get_combined_matrix()
//{
//  if(!this->combined_matrix)
//  {
//    // compute combined matrix
//    if(this->n_blocks() == 1)
//      this->combined_matrix = this->blocks.at(0);
//    else
//    {
//      // number of entries of the combined matrix
//      unsigned int n_comb_entries = this->n_total_entries();
//      unsigned int n_comb_rows = this->n_total_rows();
//      unsigned int n_comb_cols = this->n_total_cols();
//
//      // we will create a sparsity structure for the combined matrix. The
//      // following two vectors are needed for the constructor
//      int * column_of_entry = new int[n_comb_entries];
//      int * entries_in_rows = new int[n_comb_rows+1];
//      std::vector<double> comb_entries(n_comb_entries, 0.0);
//      entries_in_rows[0] = 0;
//
//      // filling the vectors:
//      unsigned int row_offset = 0;
//      // position of current entry in combined matrix
//      unsigned int pos = 0;
//      // loop over all block rows of this ColoredBlockMatrix
//      for(unsigned int block_row = 0, n_block_rows = this->n_rows();
//          block_row < n_block_rows; ++block_row)
//      {
//        // number of rows in this block_row
//        unsigned int n_rows = this->block(block_row, 0).GetN_Rows();
//        // loop over all rows in this block row
//        for(unsigned int row = 0; row < n_rows; ++row)
//        {
//          unsigned int column_offset = 0;
//          // loop over all block columns of this (block) row
//          for(unsigned int block_col = 0, n_block_col = this->n_cols();
//              block_col < n_block_col; ++block_col)
//          {
//            // current matrix block
//            const TMatrix& cm = this->block(block_row, block_col);
//            const int * row_ptr = cm.GetRowPtr();
//            const int * col_ptr = cm.GetKCol();
//            const double * entries = cm.GetEntries();
//            // loop over entire row in this block
//            for(int e = row_ptr[row]; e < row_ptr[row+1]; ++e)
//            {
//              comb_entries[pos] = entries[e];
//              column_of_entry[pos] = col_ptr[e] + column_offset;
//              ++pos;
//            }
//            column_offset += cm.GetN_Columns();
//          }
//          entries_in_rows[row_offset + row + 1] = pos;
//        }
//        row_offset += n_rows;
//      }
//
//      // create sparsity structure
//      std::shared_ptr<TStructure> sp(
//          new TStructure(n_comb_rows, n_comb_cols, n_comb_entries,
//                         column_of_entry, entries_in_rows));
//      // create Matrix
//      this->combined_matrix = std::make_shared<TMatrix>(sp);
//      this->combined_matrix->setEntries(comb_entries);
//    }
//  }
//  // else reuse the already computed combined matrix
//  return this->combined_matrix;
//}
//
///* ************************************************************************* */
//std::shared_ptr<const TMatrix> ColoredBlockMatrix::block(const unsigned int i) const
//{
//  if(i >= this->n_blocks())
//  {
//    ErrThrow("There are only ", this->n_blocks(),
//             " blocks in this ColoredBlockMatrix. Cannot access block ", i);
//  }
//  return this->blocks[i];
//}
//
///* ************************************************************************* */
//std::shared_ptr<TMatrix> ColoredBlockMatrix::block(const unsigned int i)
//{
//  if(i >= this->n_blocks())
//  {
//    ErrThrow("There are only ", this->n_blocks(),
//             " blocks in this ColoredBlockMatrix. Cannot access block ", i);
//  }
//  return this->blocks[i];
//}
//
///* ************************************************************************* */
//const TMatrix& ColoredBlockMatrix::block(const unsigned int r,
//                                  const unsigned int c) const
//{
//  if(r >= this->n_rows())
//  {
//    ErrThrow("There are only ", this->n_rows(),
//             " block rows in this ColoredBlockMatrix. Can not access a block in row ",
//             r);
//  }
//  if(c >= this->n_cols())
//  {
//    ErrThrow("There are only ", this->n_cols(),
//             " block columns in this ColoredBlockMatrix. Can not access a block in ",
//             "column ", c);
//  }
//  return *(this->blocks[r * this->n_cols() + c].get());
//}
//
///* ************************************************************************* */
//unsigned int ColoredBlockMatrix::n_total_rows() const
//{
//  unsigned int n_total_rows = 0;
//  for(unsigned int i = 0; i < this->n_rows(); i++)
//    n_total_rows += this->blocks[i * this->n_cols()]->GetN_Rows();
//
//  return n_total_rows;
//}
//
///* ************************************************************************* */
//unsigned int ColoredBlockMatrix::n_total_cols() const
//{
//  unsigned int n_total_cols = 0;
//  for(unsigned int i = 0; i < this->n_cols(); i++)
//    n_total_cols += this->blocks[i]->GetN_Columns();
//
//  return n_total_cols;
//}
//
///** ************************************************************************* */
//unsigned int ColoredBlockMatrix::n_total_entries() const
//{
//  unsigned int n_total_entries = 0;
//  unsigned int n_blocks = this->n_blocks();
//  for(unsigned int i = 0; i < n_blocks; i++)
//    n_total_entries += this->blocks[i]->GetN_Entries();
//
//  return n_total_entries;
//}
//
///* ************************************************************************* */
//double & ColoredBlockMatrix::operator()(unsigned int i, unsigned int j)
//{
//  unsigned int bI = this->block_of_index(i, j);
//  return this->blocks[bI]->operator()(i,j);
//}
//
///** ************************************************************************* */
//const double & ColoredBlockMatrix::operator()(unsigned int i, unsigned int j)
//  const
//{
//  int bI = this->block_of_index(i, j);
//  return this->blocks[bI]->operator()(i,j);
//}
//
///** ************************************************************************* */
//unsigned int ColoredBlockMatrix::block_of_index(unsigned int& i, unsigned int& j)
// const
//{
//  unsigned int n_block_rows = this->n_rows();
//  unsigned int n_block_columns = this->n_cols();
//  // find index of block where the (i,j)-th entry is located in
//  for(unsigned int block_row = 0; block_row < n_block_rows; block_row++)
//  {
//    // number of rows in this block row
//    unsigned int n_rows = blocks[block_row * n_block_columns]->GetN_Rows();
//    if( i < n_rows )
//    {
//      for(unsigned int block_column = 0; block_column < n_block_columns;
//          block_column++)
//      {
//        // this block
//        unsigned int index = block_row * n_block_columns + block_column;
//        unsigned int n_cols = blocks[index]->GetN_Columns();
//        if(j < n_cols)
//        {
//          return index;
//        }
//        j -= n_cols;
//      }
//    }
//    i -= n_rows;
//  }
//  // until here only in case of an error
//  ErrThrow("could not find the given index in the ColoredBlockMatrix");
//}
//
///** ************************************************************************* */
//void ColoredBlockMatrix::info(size_t verbose) const
//{
//  this->block_pattern->info(verbose);
//  if(verbose > 0 && verbose < 2)
//  {
//    for(unsigned int b = 0, n_b = this->n_blocks(); b < n_b; ++b)
//    {
//      const TMatrix& m = *this->blocks[b];
//      if(verbose < 2)
//      {
//        Output::print<1>(" block ", b, " has ", m.GetN_Rows(), " rows and ",
//                         m.GetN_Columns(), " columns");
//      }
//      else
//        m.info(verbose - 1);
//    }
//  }
//}
//
///** ************************************************************************* */

