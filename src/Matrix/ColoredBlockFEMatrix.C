#include <ColoredBlockFEMatrix.h>
#include <FEMatrix.h>
#include <BlockVector.h>

#include <limits>

/* ************************************************************************* */
// IMPLEMENTATION OF PUBLIC METHODS
/* ************************************************************************* */

#ifdef __2D__
ColoredBlockFEMatrix::ColoredBlockFEMatrix(
    std::vector< const TFESpace2D*  > spaces) :
#elif __3D__
    ColoredBlockFEMatrix::ColoredBlockFEMatrix(
        std::vector< const TFESpace3D*  > spaces) :
#endif
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
#ifdef __2D__
    const TFESpace2D& testspace_of_row = *spaces[i];
#elif __3D__
    const TFESpace3D& testspace_of_row = *spaces[i];
#endif
    //hold the number of rows each cell in this row will have
    size_t n_rows_of_cell = testspace_of_row.GetN_DegreesOfFreedom();

    for(size_t j = 0; j < n_cell_columns_ ; ++j )
    {
      //hold the ansatzspace each cell in this column will have
#ifdef __2D__
      const TFESpace2D& ansatzspace_of_column = *spaces[j];
#elif __3D__
      const TFESpace3D& ansatzspace_of_column = *spaces[j];
#endif
      //hold the number of columns each cell in this column will have
      size_t n_columns_of_cell = ansatzspace_of_column.GetN_DegreesOfFreedom();

      // construct the new cell info and the zero matrix it will hold
      CellInfo newInfo(n_rows_of_cell, n_columns_of_cell);

      //Create a new FEMatrix block. "True" is a swatich to make an empty TStructure.
      newInfo.block_ = std::make_shared<FEMatrix>(&testspace_of_row, &ansatzspace_of_column, true);

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

ColoredBlockFEMatrix::ColoredBlockFEMatrix() :
ColoredBlockMatrix(), test_spaces_rowwise_(),
ansatz_spaces_columnwise_()
{
}

//named constructors
#ifdef __2D__
/* ************************************************************************* */
ColoredBlockFEMatrix ColoredBlockFEMatrix::CD2D( const TFESpace2D& space )
{
  ColoredBlockFEMatrix my_matrix({&space});
  //replace block by a block with non-empty TStructure
  my_matrix.replace_blocks(FEMatrix(&space), {{0,0}} , {false});
  return my_matrix;
}

ColoredBlockFEMatrix ColoredBlockFEMatrix::NSE2D_Type1( const TFESpace2D& velocity,
                                         const TFESpace2D& pressure)
{
  ColoredBlockFEMatrix my_matrix({&velocity, &velocity, &pressure});

  //create new blocks with correct structures filled with 0
  FEMatrix velo_velo(&velocity, &velocity); //A block
  FEMatrix velo_velo_zero(&velocity, &velocity, true); //velocity zero block

  FEMatrix pressure_velo_1(&pressure, &velocity);
  FEMatrix pressure_velo_2(pressure_velo_1); // copy constructed, shares TStructure!

  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo, {{0,0},{1,1}}, {false, false});
  my_matrix.replace_blocks(velo_velo_zero, {{1,0},{0,1}}, {false, false});

  // fill in the pressure_velo blocks B_1 and B_1^T
  my_matrix.replace_blocks(pressure_velo_1, {{2,0}, {0,2} }, {false, true});

  // fill in the pressure_velo blocks B_2 and B_2^T
  my_matrix.replace_blocks(pressure_velo_2, {{2,1}, {1,2} }, {false, true});

  //block (2,2) stays as default initialized

  return my_matrix;
}

ColoredBlockFEMatrix ColoredBlockFEMatrix::NSE2D_Type2(
    const TFESpace2D& velocity, const TFESpace2D& pressure)

{
  ColoredBlockFEMatrix my_matrix({&velocity, &velocity, &pressure});

  //create new blocks with correct structures filled with 0
  FEMatrix velo_velo(&velocity, &velocity); //A block
  FEMatrix velo_velo_zero(&velocity, &velocity, true); //velocity zero block

  FEMatrix pressure_velo_1(&pressure, &velocity);
  FEMatrix pressure_velo_2(pressure_velo_1); // copy constructed, shares TStructure!

  FEMatrix velo_pressure_1(&velocity, &pressure);
  FEMatrix velo_pressure_2(velo_pressure_1); // copy constructed, shares TStructure!

  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo, {{0,0},{1,1}}, {false, false});
  my_matrix.replace_blocks(velo_velo_zero, {{1,0},{0,1}}, {false, false});

  // fill in the pressure_velo blocks B_1 and B_1^T
  my_matrix.replace_blocks(pressure_velo_1, {{2,0}}, {false});
  my_matrix.replace_blocks(velo_pressure_1, {{0,2}}, {false});

  // fill in the pressure_velo blocks B_2 and B_2^T
  my_matrix.replace_blocks(pressure_velo_2, {{2,1}}, {false});
  my_matrix.replace_blocks(velo_pressure_2, {{1,2}}, {false});

  //block (2,2) stays as default initialized

  return my_matrix;
}

ColoredBlockFEMatrix ColoredBlockFEMatrix::NSE2D_Type3(
    const TFESpace2D& velocity, const TFESpace2D& pressure)
{
  ColoredBlockFEMatrix my_matrix({&velocity, &velocity, &pressure});

  //create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity); //A blocks
  FEMatrix velo_velo_0_1(velo_velo_0_0);
  FEMatrix velo_velo_1_0(velo_velo_0_0);
  FEMatrix velo_velo_1_1(velo_velo_0_0); //all copy constructed, share one TStructure

  FEMatrix pressure_velo_1(&pressure, &velocity);
  FEMatrix pressure_velo_2(pressure_velo_1); // copy constructed, shares TStructure!

  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}}, {false});
  my_matrix.replace_blocks(velo_velo_0_1, {{0,1}}, {false});
  my_matrix.replace_blocks(velo_velo_1_0, {{1,0}}, {false});
  my_matrix.replace_blocks(velo_velo_1_1, {{1,1}}, {false});


  // fill in the pressure_velo blocks B_1 and B_1^T
  my_matrix.replace_blocks(pressure_velo_1, {{2,0}, {0,2} }, {false, true});

  // fill in the pressure_velo blocks B_2 and B_2^T
  my_matrix.replace_blocks(pressure_velo_2, {{2,1}, {1,2} }, {false, true});

  //block (2,2) stays as default initialized

  return my_matrix;
}

ColoredBlockFEMatrix ColoredBlockFEMatrix::NSE2D_Type4(
    const TFESpace2D& velocity, const TFESpace2D& pressure)
{
  ColoredBlockFEMatrix my_matrix({&velocity, &velocity, &pressure});

  //create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity); //A blocks
  FEMatrix velo_velo_0_1(velo_velo_0_0);
  FEMatrix velo_velo_1_0(velo_velo_0_0);
  FEMatrix velo_velo_1_1(velo_velo_0_0); //all copy constructed, share one TStructure

  FEMatrix pressure_velo_1(&pressure, &velocity);
  FEMatrix pressure_velo_2(pressure_velo_1); // copy constructed, shares TStructure!

  FEMatrix velo_pressure_1(&velocity, &pressure);
  FEMatrix velo_pressure_2(velo_pressure_1); // copy constructed, shares TStructure!

  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}}, {false});
  my_matrix.replace_blocks(velo_velo_0_1, {{0,1}}, {false});
  my_matrix.replace_blocks(velo_velo_1_0, {{1,0}}, {false});
  my_matrix.replace_blocks(velo_velo_1_1, {{1,1}}, {false});

  // fill in the pressure_velo blocks B_1 and B_1^T
  my_matrix.replace_blocks(pressure_velo_1, {{2,0}}, {false});
  my_matrix.replace_blocks(velo_pressure_1, {{0,2}}, {false});

  // fill in the pressure_velo blocks B_2 and B_2^T
  my_matrix.replace_blocks(pressure_velo_2, {{2,1}}, {false});
  my_matrix.replace_blocks(velo_pressure_2, {{1,2}}, {false});

  //block (2,2) stays as default initialized

  return my_matrix;
}

ColoredBlockFEMatrix ColoredBlockFEMatrix::NSE2D_Type14(
    const TFESpace2D& velocity, const TFESpace2D& pressure)
{
  ColoredBlockFEMatrix my_matrix({&velocity, &velocity, &pressure});

  //create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity); //A blocks
  FEMatrix velo_velo_0_1(velo_velo_0_0);
  FEMatrix velo_velo_1_0(velo_velo_0_0);
  FEMatrix velo_velo_1_1(velo_velo_0_0); //all copy constructed, share one TStructure

  FEMatrix pressure_velo_1(&pressure, &velocity);
  FEMatrix pressure_velo_2(pressure_velo_1); // copy constructed, shares TStructure!

  FEMatrix velo_pressure_1(&velocity, &pressure);
  FEMatrix velo_pressure_2(velo_pressure_1); // copy constructed, shares TStructure!

  FEMatrix pressure_pressure(&pressure, &pressure);

  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}}, {false});
  my_matrix.replace_blocks(velo_velo_0_1, {{0,1}}, {false});
  my_matrix.replace_blocks(velo_velo_1_0, {{1,0}}, {false});
  my_matrix.replace_blocks(velo_velo_1_1, {{1,1}}, {false});

  // fill in the pressure_velo blocks B_1 and B_1^T
  my_matrix.replace_blocks(pressure_velo_1, {{2,0}}, {false});
  my_matrix.replace_blocks(velo_pressure_1, {{0,2}}, {false});

  // fill in the pressure_velo blocks B_2 and B_2^T
  my_matrix.replace_blocks(pressure_velo_2, {{2,1}}, {false});
  my_matrix.replace_blocks(velo_pressure_2, {{1,2}}, {false});

  // fill in the pressure-pressure block
  my_matrix.replace_blocks(pressure_pressure, {{2,2}}, {false});

  return my_matrix;
}

ColoredBlockFEMatrix ColoredBlockFEMatrix::Darcy2D( const TFESpace2D& velocity, const TFESpace2D& pressure)
{
  ColoredBlockFEMatrix my_matrix({&velocity, &pressure});

  //fill in the blocks with correct matrices constructed solely for them
  my_matrix.replace_blocks(FEMatrix(&velocity, &velocity), {{0,0}}, {false});
  my_matrix.replace_blocks(FEMatrix(&velocity, &pressure), {{0,1}}, {false});
  my_matrix.replace_blocks(FEMatrix(&pressure, &velocity), {{1,0}}, {false});
  my_matrix.replace_blocks(FEMatrix(&pressure, &pressure), {{1,1}}, {false});

  return my_matrix;
}

ColoredBlockFEMatrix ColoredBlockFEMatrix::Mass_NSE2D(const TFESpace2D& velocity)
{
  ColoredBlockFEMatrix my_matrix({&velocity, &velocity});
  
  my_matrix.replace_blocks(FEMatrix(&velocity, &velocity), {{0,0}, {1, 1}}, 
                           {false, false});

  return my_matrix;

}
#elif __3D__
//3D named constructors
ColoredBlockFEMatrix ColoredBlockFEMatrix::CD3D( const TFESpace3D& space )
{
  ColoredBlockFEMatrix my_matrix({&space});
  //replace block by a block with non-empty TStructure
  my_matrix.replace_blocks(FEMatrix(&space), {{0,0}} , {false});
  return my_matrix;
}
#endif

/* ************************************************************************* */

void ColoredBlockFEMatrix::add_matrix_actives(
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
  //reset all values in 'y' to 0 and delegate to apply_scaled_add
  y.reset();
  apply_scaled_add(x, y, 1.0);
}
/* ************************************************************************* */

void ColoredBlockFEMatrix::apply_scaled_add(const BlockVector & x,
                                            BlockVector & y, double a) const
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
      else
      {
        current_block.multiplyTransposedActive(xv + col_offset, yv + row_offset, a);
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
        //Output::print<2>("[",i," , ",j,"] That cast did not work.");
        //to use the method in tests, this error is important.
        // if you want to use it for debugging, outcomment it
        ErrThrow("[",i," , ",j,"] That cast did not work.");
      }
      else
      {
        Output::print<2>("[",i," , ",j,"] That cast is fine!");
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

std::vector<std::shared_ptr<const FEMatrix>> ColoredBlockFEMatrix::get_blocks() const
{
  std::vector<std::shared_ptr<const FEMatrix>> block_ptrs;

  for(size_t i =0 ; i < n_cell_columns_ ; ++i)
  {
    for(size_t j =0 ; j < n_cell_rows_ ; ++j)
    {
      //juggle the pointers around until we can store it the way we want...
      std::shared_ptr<const FEMatrix> shared //cast const and FEMatrix
      = std::dynamic_pointer_cast<const FEMatrix>(cell_grid_[i][j].block_);

      // make a weak pointer from it and push it back
      block_ptrs.push_back(shared);
    }
  }

  return block_ptrs;
}

/* ************************************************************************* */

std::vector<std::shared_ptr<FEMatrix>> ColoredBlockFEMatrix::get_blocks_TERRIBLY_UNSAFE()
{
  std::vector<std::shared_ptr<FEMatrix>> block_ptrs;

  for(size_t i =0 ; i < n_cell_columns_ ; ++i)
  {
    for(size_t j =0 ; j < n_cell_rows_ ; ++j)
    {
      //juggle the pointers around until we can store it the way we want...
      std::shared_ptr<FEMatrix> shared //cast FEMatrix
      = std::dynamic_pointer_cast<FEMatrix>(cell_grid_[i][j].block_);

      // make a weak pointer from it and push it back
      block_ptrs.push_back(shared);
    }
  }

  return block_ptrs;
}

/* ************************************************************************* */

std::vector<std::shared_ptr<FEMatrix>> ColoredBlockFEMatrix::get_blocks_uniquely(
    bool include_zeroes)
{
  //put up an all-in-input vector
  std::vector<std::vector<size_t>> cells;
  for (size_t i =0; i< n_cell_rows_ ; ++i)
  {
    for (size_t j =0; j< n_cell_rows_ ; ++j)
    {
      cells.push_back({i,j});
    }
  }
  //...and let the other implementation do the work
  return get_blocks_uniquely(cells, include_zeroes);
}

/* ************************************************************************* */

std::vector<std::shared_ptr<FEMatrix>> ColoredBlockFEMatrix::get_blocks_uniquely(
    std::vector<std::vector<size_t>> cells, bool include_zeroes)
{
  // stuff the input into a tuple which includes transposed states
  // - for our purpose it is not necessary, but thus we can make use
  // of the input checking methods, which always include tranposed state
  std::vector<bool> transp(cells.size(),false);

  //form to tuples
  std::vector<grid_place_and_mode> positions =
      check_and_tupelize_vector_input(cells, transp);

  //static input editing
  check_and_edit_input(positions);

  //check for index-out-of-bound issues
  check_indices(positions);

  //check if there is a color which is only partly affected by
  // the requested return
  size_t color = 0;
  if(does_modification_require_color_split(color, positions))
  {
    ErrThrow("The blocks you request do affect one color only partly."
        " This should be avoided, cause it can lead to unwanted side effects!")
  }

  //you reached here, everything fine? then let's prepare the output
  std::vector<std::shared_ptr<FEMatrix>> blocks;
  size_t next_highest_color = 0;
  //loop through the positions vector
  for (auto it : positions)
  {
    size_t cell_row = std::get<0>(it);
    size_t cell_column = std::get<1>(it);
    size_t color = cell_grid_[cell_row][cell_column].color_;
    if(color >= next_highest_color) //unstored color found
    {
      std::shared_ptr<FEMatrix>block =
          std::dynamic_pointer_cast<FEMatrix>(cell_grid_[cell_row][cell_column].block_);

      bool is_zero = (block->GetN_Entries() == 0); // block without entries has zero-map structure
      if(include_zeroes || !is_zero)
      {//if either zeroes are to be included or not zero - store the block!
        blocks.push_back(block);
      }
      next_highest_color = color + 1;
    }
  }

  return blocks;
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

size_t ColoredBlockFEMatrix::get_n_column_actives(size_t cell_column) const
{
  return ansatz_spaces_columnwise_.at(cell_column)->GetN_ActiveDegrees();
}
/* ************************************************************************* */

size_t ColoredBlockFEMatrix::get_n_row_actives(size_t cell_row) const
{
  return test_spaces_rowwise_.at(cell_row)->GetN_ActiveDegrees();
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
  // input checking lead to some code duping...
  if(cell_positions.size() != transposed_states.size())
  {
    ErrThrow("Number of grid positions must equal"
        "number of transposed states.");
  }

  //loop over all positions
  for(size_t i = 0; i<cell_positions.size() ; ++i)
  {
    if(cell_positions[i].size() != 2)
    {
      ErrThrow("A grid position must have TWO values, habibi!");
    }

    size_t cell_row = cell_positions[i].at(0); //hold indices
    size_t cell_column = cell_positions[i].at(1);

    //check if the spaces fit the cell grid
#ifdef __2D__
    const TFESpace2D* grid_test_space = &get_test_space(cell_row,cell_column);
    const TFESpace2D* grid_ansatz_space = &get_ansatz_space(cell_row,cell_column);
    const TFESpace2D* block_test_space = new_block.GetTestSpace2D();
    const TFESpace2D* block_ansatz_space = new_block.GetAnsatzSpace2D();
#elif __3D__
    const TFESpace3D* grid_test_space = &get_test_space(cell_row,cell_column);
    const TFESpace3D* grid_ansatz_space = &get_ansatz_space(cell_row,cell_column);
    const TFESpace3D* block_test_space = new_block.GetTestSpace3D();
    const TFESpace3D* block_ansatz_space = new_block.GetAnsatzSpace3D();
#endif
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
      // make sure we do not try to store a block with non-active rows in
      //transposed state
      int block_test_space_has_non_actives =
          block_test_space->GetN_DegreesOfFreedom() -
          block_test_space->GetN_ActiveDegrees();
      if (block_test_space_has_non_actives)
      { // the testspace of the matrix block has non-active dofs!
        // that's not allowed in transposed state
        Output::print(block_test_space_has_non_actives);
          ErrThrow("I am not allowed to store an FEMatrix with "
              "test-space-non-actives in transposed state. This would lead to"
              " 'non-active columns' and thus to loss of information.",
              "(Block (",cell_row, " , ", cell_column, ")");
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
  Output::print<5>("ColoredBlockFEMatrix copy constructor!");

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

ColoredBlockFEMatrix::ColoredBlockFEMatrix(ColoredBlockFEMatrix&& other)
: ColoredBlockMatrix::ColoredBlockMatrix(std::move(other)), //base class move
  test_spaces_rowwise_(std::move(other.test_spaces_rowwise_)),
  ansatz_spaces_columnwise_(std::move(other.ansatz_spaces_columnwise_))
{
  Output::print<5>("ColoredBlockFEMatrix move constructor!");

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
        treated_colors[color] = create_block_shared_pointer(*this->cell_grid_[i][j].block_);
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
  Output::print<5>("ColoredBlockFEMatrix swap!");

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
  Output::print<5>("ColoredBlockFEMatrix copy assignment!");

  //do a swap with the copy constructed object "other"
  swap(*this, other);

  return *this;
}

/* ************************************************************************* */
#ifdef __2D__
const TFESpace2D& ColoredBlockFEMatrix::get_test_space(size_t cell_row, size_t cell_column) const
#elif __3D__
const TFESpace3D& ColoredBlockFEMatrix::get_test_space(size_t cell_row, size_t cell_column) const
#endif
{
  if(cell_column >= n_cell_columns_) //just to not let the cell_column go unnoticed
  {
    ErrThrow("That cell_column is out of bounds.")
  }
  return *test_spaces_rowwise_.at(cell_row);
}

/* ************************************************************************* */
#ifdef __2D__
const TFESpace2D& ColoredBlockFEMatrix::get_ansatz_space(size_t cell_row, size_t cell_column) const
#elif __3D__
const TFESpace3D& ColoredBlockFEMatrix::get_ansatz_space(size_t cell_row, size_t cell_column) const
#endif
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
