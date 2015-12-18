#include <ColoredBlockFEMatrix.h>
#include <FEMatrix.h>

#include <limits>

/* ************************************************************************* */
// IMPLEMENTATION OF PUBLIC METHODS
/* ************************************************************************* */


ColoredBlockFEMatrix::ColoredBlockFEMatrix(
    std::vector<size_t> cell_row_numbers, std::vector<size_t> cell_column_numbers)
: ColoredBlockMatrix(cell_row_numbers, cell_column_numbers )
{
  // should there ever be a default FEMatrix,
  // fill with FEMatrix-zero-blocks.
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
: ColoredBlockMatrix::ColoredBlockMatrix(other)
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
        //treated_colors[color].reset(new TMatrix(*other.cell_grid_[i][j].block_));

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

  Output::print("Called derived class copy and store");
  try
  {
    const FEMatrix& fe_block_ref = dynamic_cast<const FEMatrix&>(block);
    std::shared_ptr<TMatrix> fe_block_ptr= std::make_shared<FEMatrix>(fe_block_ref);
    return fe_block_ptr;
  }
  catch (std::bad_cast e)
  {
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
