#include <ColoredBlockFEMatrix.h>
#include <FEMatrix.h>
#include <BlockVector.h>

#include <limits>

/* ************************************************************************* */
// IMPLEMENTATION OF PUBLIC METHODS
/* ************************************************************************* */


ColoredBlockFEMatrix::ColoredBlockFEMatrix(
    std::vector< const TFESpace2D*  > spaces) :
    ColoredBlockMatrix(), //base class object is default (empty) constructed
    test_spaces_rowwise_(spaces),
    ansatz_spaces_columnwise_(spaces)
{
  // class invariant: testspaces are not allowed to hold hanging nodes,
  // as the only kind of non-active dofs this class can handle is Dirichlet dofs
  for(auto sp : spaces)
  {
    if (sp->GetN_SlaveDegrees() != 0 )
    {//there is  slave dofs, i.e. hanging nodes, quit
      ErrThrow("BlockFEMatrix cannot handle hanging nodes so far! "
          "Use spaces without hanging nodes.");
    }
  }

  //reset grid fittingly
  n_cell_rows_ = spaces.size();
  n_cell_columns_ = spaces.size();
  cell_grid_ = std::vector<std::vector<CellInfo>>(n_cell_rows_, std::vector<CellInfo>(n_cell_columns_));
  color_count_ = std::vector<size_t>(); //reset color counter vector

  //traverse the cell info grid from left to right, top to bottom
  // and fill in newly constructed, correctly dimensioned zero matrices as blocks
  for(size_t i = 0; i < n_cell_rows_ ; ++i )
  {
    //hold the testspace each cell in this row will have
    const TFESpace2D& testspace_of_row = *spaces[i];
    //hold the number of rows each cell in this row will have
    size_t n_rows_of_cell = testspace_of_row.GetN_DegreesOfFreedom();

    for(size_t j = 0; j < n_cell_columns_ ; ++j )
    {
      //hold the ansatzspace each cell in this column will have
      const TFESpace2D& ansatzspace_of_column = *spaces[j];
      //hold the number of columns each cell in this column will have
      size_t n_columns_of_cell = ansatzspace_of_column.GetN_DegreesOfFreedom();

      // construct the new cell info and the zero matrix it will hold
      CellInfo newInfo(n_rows_of_cell, n_columns_of_cell);

      //Create a new FEMatrix block.
      newInfo.block_ = std::make_shared<FEMatrix>(&testspace_of_row, &ansatzspace_of_column);

      // the next color is get_n_colors() (because the colors 0 to
      // get_n_colors() -1 are already assigned)
      newInfo.color_ = get_n_colors();

      //start as non-transposed
      newInfo.is_transposed_ = false;

      // put the new cell info to the correct place by copy assignment
      cell_grid_ [i][j] = newInfo;

      // one new block has been added and colored
      color_count_.push_back(1);

    }
  }


}

/* ************************************************************************* */

void ColoredBlockFEMatrix::add_scaled_actives(
    const FEMatrix& summand, double factor,
    const std::vector<std::vector<size_t>>& cell_positions,
    const std::vector<bool>& transposed_states)
{
  std::vector<std::tuple<size_t, size_t, bool>> input_tuples(
      check_and_tupelize_vector_input(cell_positions, transposed_states));


  add_scaled_actives(summand, factor, input_tuples);
}
/* ************************************************************************* */

void ColoredBlockFEMatrix::apply(const BlockVector & x, BlockVector & y) const
{
  //check if the vectors fit, if not so the program throws an error
  check_vector_fits_pre_image(x);
  check_vector_fits_image(y);

  //tests passed: reset all values in 'y' to 0 and delegate to apply_scaled_add
  apply_scaled_add(x, y, 1.0);
}
/* ************************************************************************* */

void ColoredBlockFEMatrix::apply_scaled_add(const BlockVector & x, BlockVector & y,
                      double a) const
{ //check if the vectors fit, if not so the program throws an error
  check_vector_fits_pre_image(x);
  check_vector_fits_image(y);

  const double * xv = x.get_entries(); // array of values in x
  double * yv = y.get_entries(); // array of values in y
  size_t row_offset = 0;
  // n_rows, n_cols are the number of cell rows/columns
  for(size_t i = 0; i < n_cell_rows_; ++i)
  {
    int col_offset = 0;
    for(size_t j = 0; j < n_cell_columns_; j++)
    {
      const FEMatrix& current_block = dynamic_cast<const FEMatrix&>(*cell_grid_[i][j].block_);
      bool transp_state = cell_grid_[i][j].is_transposed_;
      //non-transposed case
      if(transp_state == false)
      {
        current_block.multiplyActive(xv + col_offset, yv + row_offset, a);
      }
      else // in transposed state there is no active- non-active distinction,
           // since in cells with non-actives-containing test space no tranposed
           // matrices are stored!
      {
        current_block.transpose_multiply(xv + col_offset, yv + row_offset, a);
      }
      col_offset += cell_grid_[i][j].n_columns_;
    }
    row_offset += cell_grid_[i][0].n_rows_;
  }

  //...finally update the non-active rows.
  y.addScaledNonActive(x,a);

}
/* ************************************************************************* */

void ColoredBlockFEMatrix::check_pointer_types()
{
  for(size_t i = 0; i < n_cell_rows_;++i)
  {
    for(size_t j = 0; j < n_cell_columns_;++j)
    {
      auto bar = cell_grid_[i][j].block_;
      auto foo = std::dynamic_pointer_cast<FEMatrix>(bar);
      if(! foo)
      {
        Output::print("[",i," , ",j,"] That cast did not work.");
      }
      else
      {
        Output::print("[",i," , ",j,"] That cast is fine!");
      }
    }
  }
}
/* ************************************************************************* */

void ColoredBlockFEMatrix::check_vector_fits_image(const BlockVector& b) const
{
  //let the base class figure out block dimensions
  ColoredBlockMatrix::check_vector_fits_image(b);

  //check if non-actives fit
  for(size_t i = 0; i<b.n_blocks(); ++i)
  {//each vector block must have as many non-actives as the testspace of
   // the corresponding matrix row
    if((int)b.n_non_actives(i) != test_spaces_rowwise_.at(i)->GetN_Dirichlet())
    {
      ErrThrow("Number of non-actives in Block ", i, " of image BlockVector "
               "does not equal number of non-actives in row's testspace. ",
               b.n_non_actives(i), " != ", test_spaces_rowwise_.at(i)->GetN_Dirichlet());
    }
  }
}
/* ************************************************************************* */

void ColoredBlockFEMatrix::check_vector_fits_pre_image(const BlockVector& x) const
{
  //let the base class figure out block dimensions
  ColoredBlockMatrix::check_vector_fits_pre_image(x);

  //check if non-actives fit
  //each vector block must have as many non-actives as the testspace of
  // the corresponding matrix row (only due to symmetric spaces +
  // no hanging nodes condition,
  // otherwise the concept of non-actives in the factor vector does not make
  // any sense at all!!)
  for(size_t i = 0; i<x.n_blocks(); ++i)
  {
    if((int)x.n_non_actives(i) != test_spaces_rowwise_.at(i)->GetN_Dirichlet())
    {
      ErrThrow("Number of non-actives in Block ", i, " of pre-image BlockVector "
               "does not equal number of non-actives in row's testspace. ",
               x.n_non_actives(i), " != ", test_spaces_rowwise_.at(i)->GetN_Dirichlet());
    }
  }

}
/* ************************************************************************* */


std::shared_ptr<TMatrix> ColoredBlockFEMatrix::get_combined_matrix() const
{
  //let the base class put up the combined matrix as it can
  std::shared_ptr<TMatrix> combined_matrix = ColoredBlockMatrix::get_combined_matrix();

  // ...and revise all the dirichlet rows!

  //get pointers in the matrix
  const int* rowptr = combined_matrix->GetRowPtr();
  const int* kcolptr = combined_matrix->GetKCol();
  double* entries = combined_matrix->GetEntries();

  size_t row_offset = 0;

  //loop through all cell rows
  for(size_t i =0; i < this->n_cell_rows_ ;++i)
  {
    size_t n_actives = this->test_spaces_rowwise_.at(i)->GetN_ActiveDegrees();
    size_t n_non_actives = this->test_spaces_rowwise_.at(i)->GetN_Dirichlet();
    size_t n_local_rows =  n_actives + n_non_actives;

    //loop through all non-active rows
    for(size_t local_row = n_actives; local_row < n_local_rows; ++local_row)
    {
      size_t global_row = row_offset + local_row;
      size_t start = rowptr[global_row];
      size_t end = rowptr[global_row + 1];
      Output::print("b): ", global_row, " ", start , " ", end);

      //loop through all entries in the current global row
      for (size_t index = start; index < end; ++index)
      {
        // 1 on global diagonal, 0 elsewhere
        entries[index] = ( kcolptr[index] == (int) global_row ? 1 : 0);
      }

    }
    // go one block row further
    row_offset += n_local_rows;
  }

  return combined_matrix;
}


/* ************************************************************************* */

void ColoredBlockFEMatrix::replace_blocks(
    const TMatrix& new_block,
    const std::vector<std::vector<size_t>>& cell_positions,
    const std::vector<bool>& transposed_states)
{
  ErrThrow("Don't try to put a TMatrix into an FEMatrix!");
}
/* ************************************************************************* */

void ColoredBlockFEMatrix::replace_blocks(
    const FEMatrix& new_block,
    const std::vector<std::vector<size_t>>& cell_positions,
    const std::vector<bool>& transposed_states)
{
  // input checking lead to some code duping... TODO can be done by check_and_tupelize!
  if(cell_positions.size() != transposed_states.size())
  {
    ErrThrow("Number of grid positions must equal"
        "number of transposed states.");
  }

  for(size_t i = 0; i<cell_positions.size() ; ++i)
  {
    if(cell_positions[i].size() != 2)
    {
      ErrThrow("A grid position must have TWO values, habibi!");
    }

    // check if we do not try to store a matrix with testspace
    // non-actives as transposed
    size_t cell_row = cell_positions[i][0];
    size_t cell_column = cell_positions[i][1];
    if (get_test_space(cell_row, cell_column).GetN_ActiveDegrees() !=
        get_test_space(cell_row, cell_column).GetN_DegreesOfFreedom() )
    { //the testspace of the row has non-active dofs!
      if(transposed_states.at(i))
      {//this is not allowed!!!
        ErrThrow("I am not allowed to store an FEMatrix with "
            "test-space-non-actives in transposed state!")
      }
    }

    //check if the spaces fit
    const TFESpace2D* grid_test_space = &get_test_space(cell_row,cell_column);
    const TFESpace2D* grid_ansatz_space = &get_ansatz_space(cell_row,cell_column);
    const TFESpace2D* block_test_space = new_block.GetTestSpace2D();
    const TFESpace2D* block_ansatz_space = new_block.GetAnsatzSpace2D();
    if (!transposed_states.at(i))
    {//non-transposed state
      //check object identity by adress
      if (grid_test_space != block_test_space)
      {
        ErrThrow("Test spaces are not identical at (",cell_row,cell_column,")");
      }
      if (grid_ansatz_space != block_ansatz_space)
      {
        ErrThrow("Ansatz spaces are not identical at (",cell_row,cell_column,")");
      }
    }
    else
    {//transposed state
      //check object identity by adress
      if (grid_test_space != block_ansatz_space)
      {
        ErrThrow("Grid test space does not match transposed "
            "block's ansatz space at (",cell_row,cell_column,")");
      }
      if (grid_ansatz_space != block_test_space)
      {
        ErrThrow("Grid ansatz space does not match transposed "
            "block's test space at (",cell_row,cell_column,")");
      }
    }//space fitting check done

  }


  // if everything is alright with this class, do the block replacement in the base class
  ColoredBlockMatrix::replace_blocks(
      new_block, cell_positions, transposed_states);

}
/* ************************************************************************* */

void ColoredBlockFEMatrix::scale_blocks_actives(
    double factor,
    const std::vector<std::vector<size_t>>& cell_positions )
{
  std::vector<bool> transposed_states (cell_positions.size(), false);
  std::vector<std::tuple<size_t, size_t, bool>> input_tuples(
      check_and_tupelize_vector_input(cell_positions, transposed_states));

  scale_blocks_actives(factor, input_tuples);
}
/* ************************************************************************* */

/* ************************************************************************* */
// IMPLEMENTATION OF SPECIAL MEMBER FUNCTION(S)
/* ************************************************************************* */

ColoredBlockFEMatrix::ColoredBlockFEMatrix(const ColoredBlockFEMatrix& other)
: ColoredBlockMatrix::ColoredBlockMatrix(other),
  test_spaces_rowwise_(other.test_spaces_rowwise_),
  ansatz_spaces_columnwise_(other.ansatz_spaces_columnwise_)
{
  // each block instance has to be copied once and all the shared pointers of
  // the same color have to be set pointing to the new block
  std::vector<std::shared_ptr<TMatrix>> treated_colors{color_count_.size(), nullptr};
  for ( size_t i = 0; i < n_cell_rows_ ; ++i)
  {
    for ( size_t j = 0; j < n_cell_columns_; ++j)
    {
      size_t color = cell_grid_[i][j].color_;
      if(!treated_colors[color])
      { // this is our sign to make a copy
        treated_colors[color] = create_block_shared_pointer(*other.cell_grid_[i][j].block_);
        cell_grid_[i][j].block_ = treated_colors[color];
      }
      else
      { // a pointer is stored already
        cell_grid_[i][j].block_ = treated_colors[color];
      }
    }
  }

}

/* ************************************************************************* */
void swap(ColoredBlockFEMatrix& first, ColoredBlockFEMatrix& second)
{
  std::swap(first.n_cell_columns_, second.n_cell_columns_);
  std::swap(first.n_cell_rows_, second.n_cell_rows_);
  std::swap(first.cell_grid_, second.cell_grid_);
  std::swap(first.color_count_, second.color_count_);
  std::swap(first.ansatz_spaces_columnwise_, second.ansatz_spaces_columnwise_);
  std::swap(first.test_spaces_rowwise_, second.test_spaces_rowwise_);

}

/* ************************************************************************* */
ColoredBlockFEMatrix& ColoredBlockFEMatrix::operator=(ColoredBlockFEMatrix other)
{
  //do a swap with the copy constructed object "other"
  swap(*this, other);

  return *this;;
}

/* ************************************************************************* */
const TFESpace2D& ColoredBlockFEMatrix::get_test_space(size_t cell_row, size_t cell_column) const
{
  if(cell_column >= n_cell_columns_) //just to not let the cell_column go unnoticed
  {
    ErrThrow("That cell_column is out of bounds.")
  }
  return *test_spaces_rowwise_.at(cell_row);
}

/* ************************************************************************* */
const TFESpace2D& ColoredBlockFEMatrix::get_ansatz_space(size_t cell_row, size_t cell_column) const
{
  if(cell_row >= n_cell_rows_)
  {
    ErrThrow("That cell_row is out of bounds.")
  }
  return *ansatz_spaces_columnwise_.at(cell_column);
}

/* ************************************************************************* */

/* ************************************************************************* */
// IMPLEMENTATION OF PRIVATE METHODS
/* ************************************************************************* */

// Unfortunately add_scaled_actives and scale_blocks_actives
// are code dupes of the non-active base class methods...
void ColoredBlockFEMatrix::add_scaled_actives(
    const FEMatrix& summand, double scaling_factor,
    std::vector<grid_place_and_mode> row_column_transpose_tuples)
{
  // first of all check the input, modify if reparable or throw if not so.
  check_grid_fit(summand, row_column_transpose_tuples);

  // skip checking of FESpaces - this is no big deal here, since
  // FEMatrix::addActive must anyway figure out if the adding can be done

  // check if the replacement requires color splits and if so perform them.
  size_t colorToSplit = std::numeric_limits<size_t>::max();
  while (does_modification_require_color_split(colorToSplit, row_column_transpose_tuples))
  {
    mark_for_color_split(colorToSplit, row_column_transpose_tuples);
    split_color(colorToSplit);
  }

  // after the input is order nicely and the colors split,
  // check if all transposed states of the input match
  // the transposed states of the cells
  compare_transposed_mode(row_column_transpose_tuples);

  size_t searched_color = 0;


  // delegate the additions to the TMatrices
  for (auto it: row_column_transpose_tuples)
  {
    size_t cell_row = std::get<0>(it);
    size_t cell_column = std::get<1>(it);
    CellInfo& current_cell = cell_grid_[cell_row][cell_column];
    size_t cell_color = current_cell.color_;


    if (cell_color >= searched_color)
    { // we found an untreated color
      dynamic_cast<FEMatrix*>(current_cell.block_.get())->addActive(summand, scaling_factor);
      searched_color = cell_color + 1;
      continue;
    }


  }
}
/* ************************************************************************* */

std::shared_ptr<TMatrix> ColoredBlockFEMatrix::create_block_shared_pointer(const TMatrix& block)
{
  try
  { //try to cast the given TMatrix to an FEMatrix and make an FEMatrix copy of it
    const FEMatrix& fe_block_ref = dynamic_cast<const FEMatrix&>(block);
    std::shared_ptr<TMatrix> fe_block_ptr= std::make_shared<FEMatrix>(fe_block_ref);
    return fe_block_ptr;
  }
  catch (std::bad_cast e)
  {//cast did not work!
    ErrThrow("TMatrix given. Make sure to fill a ColoredBlockFEMatrix only with FEMatrices!");
  }
}
/* ************************************************************************* */

void ColoredBlockFEMatrix::scale_blocks_actives( double scaling_factor,
                   std::vector<grid_place_and_mode> row_column_transpose_tuples)
{
  // first of all check the input, modify if reparable or throw if not so.
  check_and_edit_input( row_column_transpose_tuples );

  //check for index-out-of-bound
  check_indices( row_column_transpose_tuples );

  // check if the replacement requires color splits and if so perform them.
  size_t colorToSplit = std::numeric_limits<size_t>::max();
  while (does_modification_require_color_split(colorToSplit, row_column_transpose_tuples))
  {
    mark_for_color_split(colorToSplit, row_column_transpose_tuples);
    split_color(colorToSplit);
  }

  size_t searched_color = 0;

  // delegate the scaling to the FEMatrices
  for (auto it: row_column_transpose_tuples)
  {
    size_t cell_row = std::get<0>(it);
    size_t cell_column = std::get<1>(it);
    CellInfo& current_cell = cell_grid_[cell_row][cell_column];
    size_t cell_color = current_cell.color_;

    if (cell_color >= searched_color)
    { // we found an untreated color
      dynamic_cast<FEMatrix*>(current_cell.block_.get())->scaleActive( scaling_factor );
      searched_color = cell_color + 1;
    }
  }
}
/* ************************************************************************* */
