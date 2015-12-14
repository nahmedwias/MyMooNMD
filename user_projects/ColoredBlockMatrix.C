
#include <ColoredBlockMatrix.h>
#include <BlockVector.h>
#include <LinAlg.h>

#include <limits>
#include <algorithm>
#include <list>


/* ************************************************************************* */
ColoredBlockMatrix::CellInfo::CellInfo()
: ColoredBlockMatrix::CellInfo::CellInfo(0, 0) //delegate construction
{
  //CB DEBUG
  //Output::print("Debug: Default constructing CellInfo");
  //END DEBUG
}

/* ************************************************************************* */
ColoredBlockMatrix::CellInfo::CellInfo(size_t nRows, size_t nColumns)
: n_rows_{nRows}, n_columns_{nColumns},
  //block_in_cell gets default initialised to null
  color_{std::numeric_limits<size_t>::max()}, // no colour
  is_transposed_{false}, //non-transposed
  re_color_{ReColoringFlag::KEEP}

  {
    //CB DEBUG
    //Output::print("Debug: Constructing CellInfo with n_rows_ ", n_rows_, " and n_columns_ ", n_columns_);
    //END DEBUG
  }

  /* ************************************************************************* */
  // IMPLEMENTATION OF PUBLIC METHODS
  /* ************************************************************************* */

  /* ************************************************************************* */
  ColoredBlockMatrix::ColoredBlockMatrix(std::vector<size_t> cellRowNumbers,
                                         std::vector<size_t> cellColumnNumbers)
  : n_block_rows_{cellRowNumbers.size()},
    n_block_columns_{cellColumnNumbers.size()},
    cell_info_grid_(n_block_rows_, std::vector<CellInfo>(n_block_columns_))
    // color_count_ gets default initialized as vector of size zero
    // combined_matrix gets default initialized to smart nullptr
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
          newInfo.block_ = std::make_shared<TMatrix>(nRowsOfCell, nColumnsOfCell);

          // the next color is get_n_colors() (because the colors 0 to
          // get_n_colors() -1 are already assigned)
          newInfo.color_ = get_n_colors();

          //start as non-transposed
          newInfo.is_transposed_ = false;

          // put the new cell info to the correct place by copy assignment
          cell_info_grid_ [i][j] = newInfo;

          // one new block has been added and colored
          color_count_.push_back(1);

        }
      }
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::add_matrix_to_blocks(const TMatrix& summand,
                                                  std::vector<grid_place_and_mode> row_column_transpose_tuples)
    {
      add_scaled_matrix_to_blocks(summand, 1.0, row_column_transpose_tuples);
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::add_scaled_matrix_to_blocks(
        const TMatrix& summand, double scaling_factor,
        std::vector<grid_place_and_mode> row_column_transpose_tuples)
    {
      //check and repair input (but give runtime warning!)
      // sort and remove duplicates (special function, pass reference to the vector,
      // especially check that no two elements with same row and column have different transpose)

      // check if dimensions fit to the cell where the matrix is supposed
      // to be added with the mode (special function, pass reference to the vector)

      //size_t color_to_split
      // while adding would lead to color class split( bool modification_requires_color_split(size_t& color_to_split) )
      // perform color class split for next class
      // CB DEBUG
      // checkColoring();
      // END DEBUG
      // endwhile (must finish because at one point all matrices have their own color)

      //perform actual addition
      // durchlauf durch row_column_transpose_tuples, betroffene farben heraussuchen
      // und in einem array zwischenspeichern: std::vector<bool> recieved_modification
      // wenn dort an der Stelle meiner Farbe "false" steht -> modifizieren!
      // wenn nicht: naechsten Kandidaten nehmen.
      ;
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::check_coloring() const
    {
      // check if the color_count_ array holds
      // the correct information
      check_color_count();

      // check the ascending order of the first appearing colors
      check_coloring_order();

      // check the equivalence of the relations "has the same color as"
      // and "holds a pointer to the same matrix as"
      check_equivalence_of_relations();

      // check if the recoloring flags are all in neutral state
      check_re_coloring_flags();

    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::print_and_check() const
    {
      ;
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::print_coloring_count() const
    {
      ;
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::print_coloring_pattern() const
    {
      ;
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::replace_blocks(
        const TMatrix& new_block,
        std::vector<grid_place_and_mode> row_column_transpose_tuples)
    {
      // first of all check the input, modify if reparable or throw if not so.
      check_grid_fit(new_block, row_column_transpose_tuples);

      // check if the replacement requires color splits and if so,
      // perform them.
      size_t colorToSplit = std::numeric_limits<size_t>::max();
      while (does_modification_require_color_split(
          colorToSplit, row_column_transpose_tuples
      ))
      {
        //set split markers
        mark_for_color_split(colorToSplit, row_column_transpose_tuples);
        // perform the actual color split
        split_color(colorToSplit);
      }

      // check if the replacement requires color merges and if so,
      // perform them.


      // now everything works out - do the actual replacement!
      // - don't forget to set the tranposed state correctly!
      ;
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::scale_blocks( double scaling_factor,
                                           std::vector<grid_place_and_mode> row_column_transpose_tuples)
    {
      ;
    }



    /* ************************************************************************* */
    // IMPLEMENTATION OF PRIVATE METHODS
    /* ************************************************************************* */

    /* ************************************************************************* */

    void ColoredBlockMatrix::check_grid_fit(
        const TMatrix& matrix,
        std::vector<grid_place_and_mode>& row_column_transpose_tuples
    ) const
    {
      // a) sort input and remove duplicates if need be
      check_and_edit_input(row_column_transpose_tuples);



      for (auto it: row_column_transpose_tuples)
      {
        // b) the index pair is not out of bounds
        size_t cell_row = std::get<0>(it);
        size_t cell_column = std::get<1>(it);
        if (cell_row > n_block_rows_ || cell_column > n_block_columns_)
        {
          ErrThrow("Cell index pair out of block matrix bounds!");
        }

        //  c) the matrix fits into the given cells
        const CellInfo& currentCell = cell_info_grid_[cell_row][cell_column];
        bool transposed = std::get<2>(it);

        if(!does_block_fit_cell(matrix, currentCell, transposed))
        {
          ErrThrow("The given matrix will not fit into cell "
              "[", std::get<0>(it) ," , ", std::get<1>(it) ,"]");
        }
      }
    }

    /* ************************************************************************* */
    bool ColoredBlockMatrix::does_block_fit_cell(
        const TMatrix& block,
        const CellInfo& cell,
        bool transposed
    )
    {
      size_t n_rows_cell = cell.n_rows_;
      size_t n_columns_cell = cell.n_columns_;

      size_t n_rows_block = transposed ? block.GetN_Columns() : block.GetN_Rows();
      size_t n_columns_block = transposed ? block.GetN_Rows() : block.GetN_Columns();

      // if both dimension fit return true, otherwise false
      return (n_rows_cell == n_rows_block && n_columns_cell ==  n_columns_block);
    }

    /* ************************************************************************* */
    bool ColoredBlockMatrix::does_color_match_block(const ColoredBlockMatrix::CellInfo& first,
                                                    const ColoredBlockMatrix::CellInfo& second)
    {
      return (first.color_ == second.color_) == (first.block_ == second.block_);
    }


    /* ************************************************************************* */
    void ColoredBlockMatrix::check_and_edit_input(
        std::vector<grid_place_and_mode>& row_column_transpose_tuples)
    {
      // check if the input is
      //  a) sorted w.r.t. two first indices

      // a lambda function for alphanumerical sort
      auto sort_alphnum = [] (grid_place_and_mode a, grid_place_and_mode b) -> bool {
        if (std::get<0>(a) < std::get<0>(b))
        {
          return true;
        }
        else if (std::get<0>(a) > std::get<0>(b))
        {
          return false;
        }
        else
        {
          return (std::get<1>(a) < std::get<1>(b));
        }
      };
      std::sort(row_column_transpose_tuples.begin(),
                row_column_transpose_tuples.end(),
                sort_alphnum);

      //  b) unique w.r.t. two first indices - throw if there is an ambiguity of the
      //     type "do for transposed and non-transposed"

      // lambda function for that purpose
      auto unique_wrt_position = [] (grid_place_and_mode a, grid_place_and_mode b) -> bool {
        // catch the case of one block being modified in transp and non-tranps mode
        if( std::get<0>(a) == std::get<0>(b)
            && std::get<1>(a) == std::get<1>(b)
            && std::get<2>(a) != std::get<2>(b))
        {
          ErrThrow("Block can not be"
              " modified in transposed AND non-transposed state!"
              " [", std::get<0>(a) ," , ", std::get<1>(a) ,"]" );
        }
        return (std::get<0>(a) == std::get<0>(b)
            && std::get<1>(a) == std::get<1>(b));
      };

      // move duplicates to the end of the vector and hold an iterator
      // to the first duplicate-position
      auto unique_end = std::unique(row_column_transpose_tuples.begin(),
                                    row_column_transpose_tuples.end(),
                                    unique_wrt_position);

      // give a warning if there was duplicate input
      if(unique_end != row_column_transpose_tuples.end())
      {
        Output::print("Warning! Duplicates in the input vector removed.");
      }

      //crop the vector to range of uniques
      row_column_transpose_tuples.erase(
          unique_end, row_column_transpose_tuples.end());

    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::check_color_count() const
    {
      //a list of the cell colors from left to right, top to bottom
      std::list<size_t> foundColors;
      for(size_t i = 0 ; i < n_block_rows_ ; ++i)
      {
        for(size_t j = 0; j < n_block_columns_; ++j)
        {
          foundColors.push_back(cell_info_grid_[i][j].color_);
        }
      }

      // check for every color if it is assigned and counted as often as
      // color_count_ states
      for(size_t color = 0; color < color_count_.size() ; ++color )
      {

        // count how often color appears in foundColors
        size_t n_finds = std::count(foundColors.begin(), foundColors.end(), color);

        if(n_finds == 0) //throw if something's wrong
        {
          ErrThrow("Here is an unassigned color: ", color);
        }
        if( n_finds != color_count_[color] )
        {
          ErrThrow("Number of found colors and stored color_count_"
              "do not match for color: ", color);
        }

      }
      //CB DEBUG
      Output::print("ColoredBlockMatrix passed test: check_color_count");
      //END DEBUG
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::check_equivalence_of_relations() const
    {
      //traverse the matrix to get the first element for the comparison
      for(size_t i_first = 0 ; i_first < n_block_rows_ ; ++i_first)
      {
        for(size_t j_first = 0; j_first < n_block_columns_; ++j_first)
        {
          if( is_last_index_pair(i_first, j_first) )
          {
            //last index pair reached
            //CB DEBUG
            Output::print("ColoredBlockMatrix passed test: check_equivalence_of_relations");
            //END DEBUG
            return;
          }

          //traverse the matrix to get the second element for the comparison

          // start in the same row as first element, unless that lies in the last column:
          // then go one row further down
          size_t i_second{i_first};
          size_t j_second{j_first};

          get_next_cell_grid_index(i_second, j_second);

          for( ; i_second < n_block_rows_ ; ++i_second)
          {
            // start one column right from first element, unless we are already in
            // a further donw row. then start at zero
            for( ; j_second < n_block_columns_; ++j_second)
            {
              const CellInfo first = cell_info_grid_[i_first][j_first];
              const CellInfo& second = cell_info_grid_[i_second][j_second];
              if (! does_color_match_block(first, second) )
              {
                ErrThrow(" The equaivalence of relations is broken in this BlockMatrix.");
              }
            }
          }
        }
      }
    }

    /* ************************************************************************* */

    void ColoredBlockMatrix::check_coloring_order() const
    {
      //we rely on condition 1 here
      size_t searched_for = 0;

      for(size_t i = 0 ; i < n_block_rows_ ; ++i)
      {
        for(size_t j = 0; j < n_block_columns_; ++j)
        {
          if (cell_info_grid_[i][j].color_ > searched_for)
          {
            ErrThrow("The color ordering of this ColoredBlockMatrix is incorrect.")
          }
          else if (cell_info_grid_[i][j].color_ == searched_for)
          {
            ++searched_for;
          }
        }
      }
      //CB DEBUG
      Output::print("ColoredBlockMatrix passed test: check_coloring_order");
      //END DEBUG
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::check_re_coloring_flags() const
    {
      for(size_t i = 0 ; i < n_block_rows_ ; ++i)
      {
        for(size_t j = 0; j < n_block_columns_; ++j)
        {
          if (cell_info_grid_[i][j].re_color_ != CellInfo::ReColoringFlag::KEEP)
          {
            ErrThrow("Oops! Somebody forgot to clean up a"
                "recoloring flag in cell [",i,",",j,"].");
          }
        }
      }
      //CB DEBUG
      Output::print("ColoredBlockMatrix passed test: check_re_coloring_flags");
      //END DEBUG
    }

    /* ************************************************************************* */
    bool ColoredBlockMatrix::does_modification_require_color_split(size_t& color_to_split,
                                                                   std::vector<grid_place_and_mode> row_column_transposed_tuples) const
    {
      // the part from here performs in O(m log m), when m is size of row_column_tuples

      // fill a list with the colors which will be affected
      // by the modification and sort it
      std::list<size_t> color_touches;
      for(auto it : row_column_transposed_tuples)
      {
        size_t color = cell_info_grid_[std::get<0>(it)][std::get<1>(it)].color_;
        color_touches.push_back( color );
      }
      color_touches.sort();

      // count out how many times each affected color gets affected
      // and compare to the number of its overall appearances
      size_t color = color_touches.front();
      size_t n_touches = 0;
      for (auto it : color_touches)
      {
        if(it != color) // finished counting of a certain color
        {
          if (n_touches - color_count_[color] != 0)
          { // we found a color class which is not entirely affected by
            // the modification - this class must undergo the splitting procedure

            // update the output and return true
            color_to_split = color;
            return true;
          }
          else
          {
            //reset number of touches and color
            color = it;
            n_touches = 0;
          }

        }
        ++n_touches;
      }

      // great, now color split needed - we finished in O(m log m)
      return false;
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::mark_for_color_split(size_t color_to_mark,
                                                  std::vector<grid_place_and_mode>
    row_column_transposed_tuples)
    {
      for (auto it : row_column_transposed_tuples)
      {
        if (cell_info_grid_[std::get<0>(it)][std::get<1>(it)].color_ == color_to_mark)
        {
          //the cell is of the color to split and belongs to the given subset
          // - so mark it for split
          cell_info_grid_[std::get<0>(it)][std::get<1>(it)].re_color_ =
              CellInfo::ReColoringFlag::SPLIT;
        }
      }
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::split_color(size_t color_to_split)
    {
      // a new color is going to appear - start with 0 count
      color_count_.push_back(0);

      size_t i = 0;
      size_t j = 0;

      find_first_appearance_of_color(color_to_split, i , j);

      //deep copy the matrix and make the shared_ptr temporarily responsible
      std::shared_ptr<TMatrix> matrix_copy(new TMatrix(*cell_info_grid_[i][j].block_.get()));

      // Whether the first found cell of the color to split
      // has flag SPLIT or flag KEEP determines, which part of the color
      // class gets the new numbering - the one marked with value
      // equal to first_split maintains its current color number
      CellInfo::ReColoringFlag first_split(cell_info_grid_[i][j].re_color_);

      if(is_last_index_pair(i,j))
      {
        ErrThrow("It should not happen that the first found element of a split "
            "color is the last cell grid entry. Something's very wrong!");
      }

      size_t new_color = color_to_split + 1;
      bool searching_first_cell_of_new_class{true};

      // traverse through the remainder of the matrix (after the first appearance)
      while (! is_last_index_pair(i,j) )
      {
        get_next_cell_grid_index(i,j); //iterate up

        CellInfo& current_cell = cell_info_grid_[i][j];


        if(current_cell.color_ >= new_color) // a cell of a higher color
        {//found a cell of higher color
          if (searching_first_cell_of_new_class)
          {//this means we have not yet found the first element
            // of the second part of the color class to split
            // and thus still have to determine its number
            //CB DEBUG
            if(current_cell.color_ > new_color)
            {
              ErrThrow("Error when trying to figure out the color"
                  " of the new color class!")
            }
            //END DEBUG
            // update new_color, because we have found another class
            // which lies between our two split parts
            ++new_color;
          }
          else
          { //we are not searching anymore, the new_color is determined
            //so recolor and update the number of colored matrices array
            color_count_[current_cell.color_] -= 1;
            ++current_cell.color_;
            color_count_[current_cell.color_] += 1;

          }
        } //end: cell of higher color class
        else if (current_cell.color_ == color_to_split)
        {//found another element of the split color

          if (current_cell.re_color_ == first_split)
          {//the element belongs to the first set and thus keeps its color
            //just reset the re color flags
            current_cell.re_color_ = CellInfo::ReColoringFlag::KEEP;
          }
          else
          {
            //end the search, by now we found the new color!
            // (now this is done every time, but nevermind)
            searching_first_cell_of_new_class = false;

            //The element will get assigned the new color
            color_count_[current_cell.color_] -= 1;
            current_cell.color_ = new_color;
            color_count_[current_cell.color_] += 1;

            current_cell.block_ = matrix_copy; //the new shared_ptr

            //..and resets the re color flag
            current_cell.re_color_ = CellInfo::ReColoringFlag::KEEP;
          }
        } //end: cell of split color class
      } //endwhile

    }


    /* ************************************************************************* */
    bool ColoredBlockMatrix::does_replacement_require_color_merge(
        std::vector<grid_place_and_mode> row_column_transposed_tuples)
    {
      ;
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::merge_colors()
    {
      ;
    }



    /* ************************************************************************* */
    // Some methods to navigate in the cell info grid
    /* ************************************************************************* */

    void ColoredBlockMatrix::find_first_appearance_of_color(
        size_t color_to_find, size_t& block_row , size_t& block_column) const
    {
      if( color_to_find >= get_n_colors() )
      {
        ErrThrow("That color does not exist in this ColoredBlockMatrix.");
      }

      for(size_t i = 0; i < n_block_rows_ ; ++i)
      {
        for(size_t j = 0; j < n_block_columns_ ; ++j)
        {
          if (cell_info_grid_[i][j].color_ == color_to_find)
          {
            // color was found
            block_row = i;
            block_column = j;
            return;
          }
        }
      }
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::get_next_cell_grid_index (
        size_t& block_row, size_t& block_column ) const
    {

      if( block_row > n_block_rows_
          || block_column > n_block_columns_)
      {
        throw std::runtime_error("Index is out of bounds!");
      }

      if( is_last_index_pair ( block_row, block_column ) )
      {
        throw std::logic_error("Next index pair after last required!");
      }

      size_t new_row = (block_column < n_block_columns_ - 1) ? block_row : block_row + 1;
      size_t new_column = (block_row == new_row) ? block_column + 1 : 0;

      block_row = new_row;
      block_column = new_column;

    }

    /* ************************************************************************* */
    bool ColoredBlockMatrix::is_last_index_pair(size_t block_row, size_t block_column) const
    {
      if( block_row > n_block_rows_
          || block_column > n_block_columns_)
      {
        throw std::runtime_error("Index is out of bounds!");
      }

      if( block_row == n_block_rows_ -1
          && block_column == n_block_columns_ - 1)
      {
        return true;
      }
      return false;
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

