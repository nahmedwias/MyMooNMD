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
}

/* ************************************************************************* */
ColoredBlockMatrix::CellInfo::CellInfo(size_t nRows, size_t nColumns)
: n_rows_{nRows}, n_columns_{nColumns},
  //block_in_cell gets default initialised to null
  color_{std::numeric_limits<size_t>::max()}, // no colour
  is_transposed_{false}, //non-transposed
  re_color_flag_{ReColoringFlag::KEEP}

  {
  }

  /* ************************************************************************* */
  // IMPLEMENTATION OF PUBLIC METHODS
  /* ************************************************************************* */

  /* ************************************************************************* */

  ColoredBlockMatrix::ColoredBlockMatrix()
  : ColoredBlockMatrix{{0},{0}}
  {
  }

  ColoredBlockMatrix::ColoredBlockMatrix(std::vector<size_t> cellRowNumbers,
                                         std::vector<size_t> cellColumnNumbers)
  : n_cell_rows_{cellRowNumbers.size()},
    n_cell_columns_{cellColumnNumbers.size()},
    cell_grid_(n_cell_rows_, std::vector<CellInfo>(n_cell_columns_))
    // color_count_ gets default initialized as vector of size zero
    // combined_matrix gets default initialized to smart nullptr
    {
      //traverse the cell info grid from left to right, top to bottom
      // and fill in newly constructed, correctly dimensioned zero matrices as blocks
      for(size_t i = 0; i < n_cell_rows_ ; ++i )
      {
        //hold the number of rows each cell in this row will have
        size_t nRowsOfCell = cellRowNumbers[i];

        for(size_t j = 0; j < n_cell_columns_ ; ++j )
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
          cell_grid_ [i][j] = newInfo;

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
    void ColoredBlockMatrix::add_unscaled_matrix(
        const TMatrix& summand,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states)
    {
      add_matrix(summand, 1.0, cell_positions, transposed_states);
    }



    /* ************************************************************************* */
    void ColoredBlockMatrix::add_matrix(
        const TMatrix& summand, double scaling_factor,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states)
    {
      std::vector<std::tuple<size_t, size_t, bool>> input_tuples(
          check_and_tupelize_vector_input(cell_positions, transposed_states));


      add_scaled_matrix_to_blocks(summand, scaling_factor, input_tuples);

    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::apply(const BlockVector & x, BlockVector & y) const
    {
      //check if the vectors fit, if not so the program throws an error
      check_vector_fits_pre_image(x);
      check_vector_fits_image(y);

      //tests passed: reset all values in 'y' to 0 and delegate to apply_scaled_add
      apply_scaled_add(x, y, 1.0);
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::apply_scaled_add(const BlockVector & x, BlockVector & y,
                                              double a) const
    {
      //check if the vectors fit, if not so the program throws an error
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
          std::shared_ptr<TMatrix> current_block = cell_grid_[i][j].block_;
          bool transp_state = cell_grid_[i][j].is_transposed_;
          //non-transposed case
          if(transp_state == false)
          {
            current_block->multiply(xv + col_offset, yv + row_offset, a);
          }
          else
          {
            current_block->transpose_multiply(xv + col_offset, yv + row_offset, a);
          }
          col_offset += cell_grid_[i][j].n_columns_;
        }
        row_offset += cell_grid_[i][0].n_rows_;
      }
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

    ///! Check whether a BlockVector b is fit to be the rhs b of the equation Ax=b.
    void ColoredBlockMatrix::check_vector_fits_image(const BlockVector& b) const
    {
      size_t n_vec_blocks = b.n_blocks();

      //check if number of blocks fits number of block rows
      if(n_vec_blocks != n_cell_rows_)
      {
        ErrThrow("Vector blocks number does not fit cell row number. ",
                 n_vec_blocks, "!=", n_cell_rows_);
      }

      //check if each vector block's length fits the number of rows in the corresponding row
      for(size_t i = 0; i<n_vec_blocks; ++i)
      {
        if (b.length(i) != get_n_rows_in_cell(i, 0))
        {
          ErrThrow("Length of Vector Block ", i, " is ", b.length(i),
                   "which does not fit n_rows_in_cell ", get_n_rows_in_cell(i, 0));
        }
        if(b.active(i)!=0)
        {
          //maybe put to virtual method: handle_discovery_of_vector_actives
          // give a warning if the vector has actives - the matrix has definitely not!
          Output::print<2>("Warning! The BlockVector has actives, but BlockMatrix does not."
              " Did you want to use a BlockFEMatrix instead?");
        }
      }
    }

    ///! Check whether a BlockVector x is fit to be to be factor x in the equation Ax=b.
    void ColoredBlockMatrix::check_vector_fits_pre_image(const BlockVector& x) const
    {
      size_t n_vec_blocks = x.n_blocks();

      //check if number of blocks fits number of block columns
      if(n_vec_blocks != n_cell_columns_)
      {
        ErrThrow("Vector blocks number does not fit cell column number. ",
                 n_vec_blocks, "!=", n_cell_columns_);
      }

      //check if each vector block's length fits the number of columns in the corresponding column
      for(size_t j = 0; j<n_vec_blocks; ++j)
      {
        if (x.length(j) != get_n_columns_in_cell(0,j))
        {
          ErrThrow("Length of Vector Block ", j, " is ", x.length(j),
                   "which does not fit n_columns_in_cell ", get_n_columns_in_cell(0,j));
        }
        if(x.active(j)!=0)
        {
          //maybe put to virtual method: handle_discovery_of_vector_actives
          // give a warning if the vector has actives - the matrix has definitely not!
          Output::print<2>("Warning! The BlockVector has actives, but BlockMatrix does not."
              " Did you want to use a BlockFEMatrix instead?");
        }
      }
    }

    /* ************************************************************************* */
    std::shared_ptr<TMatrix> ColoredBlockMatrix::get_combined_matrix() const
    {
      // fill an array with smart pointers to blocks,
      // such thate those stored as transposed really get transposed

      std::vector<std::vector<std::shared_ptr<TMatrix>>> temp_block_grid
      (n_cell_rows_, std::vector<std::shared_ptr<TMatrix>>(n_cell_rows_,nullptr));

      // store smart pointers to the already treated transposed colors
      std::vector<std::shared_ptr<TMatrix>> treated_transp_colors{color_count_.size(), nullptr};

      for ( size_t i = 0; i < n_cell_rows_ ; ++i)
      {
        for ( size_t j = 0; j < n_cell_columns_; ++j)
        {
          if (! cell_grid_[i][j].is_transposed_)
          { // non-transposed - just store the cells smart pointer
            temp_block_grid[i][j] = cell_grid_[i][j].block_;
          }
          else
          { //transposed case - check in the treated_transp_colors vector
            size_t color = cell_grid_[i][j].color_;

            if(!treated_transp_colors[color])
            { // this is our sign to make a transposed copy
              treated_transp_colors[color].reset(cell_grid_[i][j].block_->GetTransposed());
            }

            temp_block_grid[i][j] = treated_transp_colors[color];

          } // end transp -non-transp

        }// end cell_columns
      }// end cell_rows


      // number of entries of the combined matrix
      size_t n_comb_entries = get_n_total_entries();
      size_t n_comb_rows = get_n_total_rows();
      size_t n_comb_cols = get_n_total_columns();

      // we will create a sparsity structure for the combined matrix. The
      // following two vectors are needed for the constructor
      int * column_of_entry = new int[n_comb_entries];
      int * entries_in_rows = new int[n_comb_rows+1];
      std::vector<double> comb_entries(n_comb_entries, 0.0);
      entries_in_rows[0] = 0;

      // filling the vectors:
      size_t row_offset = 0;
      // position of current entry in combined matrix
      size_t pos = 0;
      // loop over all block rows of this ColoredBlockMatrix
      for(size_t block_row = 0; block_row < n_cell_rows_; ++block_row)
      {
        // number of matrix rows in this cell row
        size_t n_rows = cell_grid_[block_row][0].n_rows_;
        // loop over all rows in this cell row
        for(size_t row = 0; row < n_rows; ++row)
        {
          size_t column_offset = 0;
          // loop over all block columns of this (block) row
          for(size_t block_col = 0; block_col < n_cell_columns_; ++block_col)
          {
            // current matrix block
            const TMatrix& cm = *temp_block_grid[block_row][block_col];
            const int * row_ptr = cm.GetRowPtr();
            const int * col_ptr = cm.GetKCol();
            const double * entries = cm.GetEntries();
            // loop over entire row in this block
            for(int e = row_ptr[row]; e < row_ptr[row+1]; ++e)
            {
              comb_entries[pos] = entries[e];
              column_of_entry[pos] = col_ptr[e] + column_offset;
              ++pos;
            }
            column_offset += cm.GetN_Columns();
          }
          entries_in_rows[row_offset + row + 1] = pos;
        }
        row_offset += n_rows;
      }

      // create sparsity structure
      std::shared_ptr<TStructure> sp(
          new TStructure(n_comb_rows, n_comb_cols, n_comb_entries,
                         column_of_entry, entries_in_rows));

      //Structure copies these object
      delete[] column_of_entry;
      delete[] entries_in_rows;

      // create Matrix
      std::shared_ptr<TMatrix> combined (new TMatrix(sp));
      combined->setEntries(comb_entries);

      return combined;
    }

    /// @brief total number of columns(> n_block_columns)
    size_t ColoredBlockMatrix::get_n_total_columns() const
    { // sum the number of columns in the first block row
      size_t n_total_columns = 0;
      for (size_t j = 0; j < n_cell_columns_ ;++j)
      {
        n_total_columns += cell_grid_[0][j].n_columns_;
      }
      return n_total_columns;
    }

    /// @brief total number of entries
    size_t ColoredBlockMatrix::get_n_total_entries() const
    {
      size_t n_total_entries = 0;
      //sum the number of entries from all cells
      for(size_t i= 0; i < n_cell_rows_ ; ++i)
      {
        for (size_t j = 0; j < n_cell_columns_ ;++j)
        {
          n_total_entries += cell_grid_[i][j].block_->GetN_Entries();
        }
      }
      return n_total_entries;
    }

    /// @brief total number of rows (> n_block_rows)
    size_t  ColoredBlockMatrix::get_n_total_rows() const
    { // sum the number of rows in the first block column
      size_t n_total_rows = 0;
      for (size_t i = 0; i < n_cell_rows_ ;++i)
      {
        n_total_rows += cell_grid_[i][0].n_rows_;
      }
      return n_total_rows;
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::print_and_check(std::string matrix_name) const
    {
      // do both prints
      print_coloring_pattern(matrix_name);
      print_coloring_count(matrix_name);

      //...and perform the consistency check
      check_coloring();
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::print_coloring_count(std::string matrix_name) const
    {
      Output::print("----------------------");
      Output::print(" Color count: ", matrix_name);
      Output::print("----------------------");
      size_t color = 0;
      for (auto count : color_count_)
      {
        Output::print(color, "\t : \t", count );
        ++color;
      }
      Output::print("----------------------");
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::print_coloring_pattern(std::string matrix_name,
                                                    bool print_adress) const
    {
      Output::print("-------------------------");
      Output::print(" Coloring pattern: ", matrix_name);
      Output::print("-------------------------");
      for(size_t i = 0; i < n_cell_rows_; ++i)
      {
        std::stringstream out_row;
        out_row << "( ";
        for(size_t j = 0; j < n_cell_columns_; ++j)
        {
          out_row << "\t";
          if(cell_grid_[i][j].block_->GetN_Entries()==0)
          {//for a zero-map block the number is printed in brackets
            out_row << "(";
          }
          out_row << cell_grid_[i][j].color_;
          if(cell_grid_[i][j].is_transposed_)
          {
            out_row << "^T";
          }
          if(cell_grid_[i][j].block_->GetN_Entries()==0)
          {
            out_row << ")";
          }
          if(print_adress)
          {
            out_row << " ";
            out_row << cell_grid_[i][j].block_.get();
          }
        }
        out_row << "\t )";
        Output::print(out_row.str());
      }
      Output::print("-------------------------");
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::replace_blocks(
        const TMatrix& new_block,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states)
    {
      std::vector<std::tuple<size_t, size_t, bool>> input_tuples(
          check_and_tupelize_vector_input(cell_positions, transposed_states));

      replace_blocks(new_block, input_tuples );
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::scale(double factor)
    {
      std::vector<std::vector<size_t>> positions;
      for ( size_t i = 0; i < n_cell_rows_ ; ++i)
      {
        for ( size_t j = 0; j < n_cell_columns_; ++j)
        {
          positions.push_back({i,j});
        }
      }
      // scale it!
      scale_blocks(factor, positions);
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::scale_blocks(
        double scaling_factor,
        const std::vector<std::vector<size_t>>& cell_positions )
    {
      std::vector<bool> transposed_states (cell_positions.size(), false);
      std::vector<std::tuple<size_t, size_t, bool>> input_tuples(
          check_and_tupelize_vector_input(cell_positions, transposed_states));

      scale_blocks(scaling_factor, input_tuples);
    }

    /* ************************************************************************* */
    // IMPLEMENTATION OF SPECIAL MEMBER FUNCTIONS
    /* ************************************************************************* */
    ColoredBlockMatrix::ColoredBlockMatrix(const ColoredBlockMatrix& other)
    : n_cell_rows_(other.n_cell_rows_), n_cell_columns_(other.n_cell_columns_),
      cell_grid_(other.cell_grid_), color_count_(other.color_count_)

    {
      Output::print<3>("Base class copy constructor.");

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
            treated_colors[color].reset(new TMatrix(*other.cell_grid_[i][j].block_));

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

    void swap(ColoredBlockMatrix& first, ColoredBlockMatrix& second)
    {
      Output::print<3>("Base class swap.");

      std::swap(first.n_cell_columns_, second.n_cell_columns_);
      std::swap(first.n_cell_rows_, second.n_cell_rows_);
      std::swap(first.cell_grid_, second.cell_grid_);
      std::swap(first.color_count_, second.color_count_);
    }
    /* ************************************************************************* */

    ColoredBlockMatrix& ColoredBlockMatrix::operator=(ColoredBlockMatrix other)
    {
      Output::print<3>("Base class copy assignment.");

      //do a swap with the copy constructed object "other"
      swap(*this, other);

      return *this;
    }

    /* ************************************************************************* */
    // IMPLEMENTATION OF PRIVATE METHODS
    /* ************************************************************************* */

    /* ************************************************************************* */
    void ColoredBlockMatrix::add_scaled_matrix_to_blocks(
        const TMatrix& summand, double scaling_factor,
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

      //      // list of colors whose first instance is not in the same transposed
      //      // state as the summand - they have to be treated specially
      //      std::list<size_t> treat_special_colors;

      // delegate the additions to the TMatrices
      for (auto it: row_column_transpose_tuples)
      {
        size_t cell_row = std::get<0>(it);
        size_t cell_column = std::get<1>(it);
        CellInfo& current_cell = cell_grid_[cell_row][cell_column];
        size_t cell_color = current_cell.color_;


        if (cell_color >= searched_color)
        { // we found an untreated color
          //          bool addend_transp = current_cell.is_transposed_;
          //          bool summand_transp = std::get<2>(it);

          //          if(summand_transp == addend_transp) // transposed states match ?
          //          {
          current_cell.block_.get()->add_scaled(summand, scaling_factor);
          //          }
          //          else
          //          { //the color has to be kept in mind for later treatment
          //            treat_special_colors.push_back(cell_color);
          //            Output::print("Special treatment color found, transp not fitting"
          //                          "at ", cell_row ,",",cell_column);
          //          }

          searched_color = cell_color + 1;
          continue;
        }

        //        // check if the color is among those which need a special treatment
        //        auto hit = std::find(treat_special_colors.begin(),
        //                             treat_special_colors.end(), cell_color);
        //        if (hit != treat_special_colors.end())
        //        {
        //          // color needs special treatment
        //          bool addend_transp = current_cell.is_transposed_;
        //          bool summand_transp = std::get<2>(it);
        //
        //          if(summand_transp == addend_transp) // transposed states match ?
        //          {
        //            current_cell.block_.get()->add_scaled(summand, scaling_factor);
        //            //color does not need special treatment anymore
        //            treat_special_colors.erase(hit);
        //            Output::print("Special treatment color removed, found fitting transp"
        //                          "at ", cell_row ,",",cell_column);
        //          }
        //        }
      }

      //      //treat those colors which are still in the treat_special_colors list
      //      // - those where the transposed state of all appearances does not fit
      //      for (auto it : treat_special_colors)
      //      {
      //        //if this is nowhere the case, we have to add the transposed somewhere
      //        {
      //          size_t cell_row;
      //          size_t cell_column;
      //          find_first_appearance_of_color( it, cell_row, cell_column );
      //
      //          CellInfo& current_cell = cell_grid_[cell_row][cell_column];
      //
      //          //make a transposed version of the summand - this is expensive
      //          TMatrix* summand_transposed = summand.GetTransposed();
      //
      //          Output::print("Warning! Must transpose a TMatrix to perform addition. \n"
      //              "Consider different way to use ColoredBlockMatrix::add_scaled_matrix_to_blocks");
      //
      //          current_cell.block_.get()->add_scaled(*summand_transposed, scaling_factor);
      //
      //          delete summand_transposed;
      //        }
      //      }
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
    std::vector<std::tuple<size_t, size_t, bool>> ColoredBlockMatrix::
    check_and_tupelize_vector_input(
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states )
    {
      if (cell_positions.size() != transposed_states.size())
      {
        ErrThrow("Different input vector lenghts!");
      }

      std::vector<grid_place_and_mode> tuples;

      for (size_t i = 0; i < cell_positions.size(); ++i)
      {
        if (cell_positions[i].size() != 2)
        {
          ErrThrow("Cell_positions at ", i, " is not an index pair!");
        }
        size_t cell_row = cell_positions[i][0];
        size_t cell_column = cell_positions[i][1];
        bool transp = transposed_states[i];

        tuples.push_back(std::make_tuple(cell_row, cell_column, transp));
      }

      return tuples;

    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::check_color_count() const
    {
      //a list of the cell colors from left to right, top to bottom
      std::list<size_t> foundColors;
      for(size_t i = 0 ; i < n_cell_rows_ ; ++i)
      {
        for(size_t j = 0; j < n_cell_columns_; ++j)
        {
          foundColors.push_back(cell_grid_[i][j].color_);
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
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::check_coloring_order() const
    {
      //we rely on condition 1 here
      size_t searched_for = 0;

      for(size_t i = 0 ; i < n_cell_rows_ ; ++i)
      {
        for(size_t j = 0; j < n_cell_columns_; ++j)
        {
          if (cell_grid_[i][j].color_ > searched_for)
          {
            ErrThrow("The color ordering of this ColoredBlockMatrix is incorrect.")
          }
          else if (cell_grid_[i][j].color_ == searched_for)
          {
            ++searched_for;
          }
        }
      }
    }



    /* ************************************************************************* */
    void ColoredBlockMatrix::check_equivalence_of_relations() const
    {
      //traverse the matrix to get the first element for the comparison
      for(size_t i_first = 0 ; i_first < n_cell_rows_ ; ++i_first)
      {
        for(size_t j_first = 0; j_first < n_cell_columns_; ++j_first)
        {
          if( is_last_index_pair(i_first, j_first) )
          {
            //last index pair reached
            return;
          }

          //traverse the matrix to get the second element for the comparison

          // start in the same row as first element, unless that lies in the last column:
          // then go one row further down
          size_t i_second{i_first};
          size_t j_second{j_first};

          get_next_cell_grid_index(i_second, j_second);

          for( ; i_second < n_cell_rows_ ; ++i_second)
          {
            // start one column right from first element, unless we are already in
            // a further donw row. then start at zero
            for( ; j_second < n_cell_columns_; ++j_second)
            {
              const CellInfo first = cell_grid_[i_first][j_first];
              const CellInfo& second = cell_grid_[i_second][j_second];
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

    void ColoredBlockMatrix::check_grid_fit(
        const TMatrix& matrix,
        std::vector<grid_place_and_mode>& row_column_transpose_tuples
    ) const
    {
      // sort input and remove duplicates if need be
      check_and_edit_input(row_column_transpose_tuples);

      // check for index out of bounds
      check_indices(row_column_transpose_tuples);

      //  check if the matrix fits into the given cells
      for (auto it: row_column_transpose_tuples)
      {
        size_t cell_row = std::get<0>(it);
        size_t cell_column = std::get<1>(it);

        const CellInfo& currentCell = cell_grid_[cell_row][cell_column];
        bool transposed = std::get<2>(it);

        if(!does_block_fit_cell(matrix, currentCell, transposed))
        {
          ErrThrow("The given matrix will not fit into cell "
              "[", std::get<0>(it) ," , ", std::get<1>(it) ,"]");
        }
      }
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::check_indices(
        std::vector<grid_place_and_mode>& row_column_transpose_tuples
    ) const
    {
      for (auto it: row_column_transpose_tuples)
      {
        size_t cell_row = std::get<0>(it);
        size_t cell_column = std::get<1>(it);
        if (cell_row >= n_cell_rows_ || cell_column >= n_cell_columns_)
        {
          ErrThrow("Cell index pair out of block matrix bounds! "
              "[", std::get<0>(it) ," , ", std::get<1>(it) ,"]");
        }
      }
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::check_re_coloring_flags() const
    {
      for(size_t i = 0 ; i < n_cell_rows_ ; ++i)
      {
        for(size_t j = 0; j < n_cell_columns_; ++j)
        {
          if (cell_grid_[i][j].re_color_flag_ != CellInfo::ReColoringFlag::KEEP)
          {
            ErrThrow("Oops! Somebody forgot to clean up a"
                "recoloring flag in cell [",i,",",j,"].");
          }
        }
      }
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::compare_transposed_mode(
        std::vector<grid_place_and_mode>& row_column_transpose_tuples) const
    {
      for ( auto it : row_column_transpose_tuples)
      {
        size_t i = std::get<0>(it);
        size_t j = std::get<1>(it);
        bool you_are_transposed = std::get<2>(it);
        bool i_am_transposed=cell_grid_[i][j].is_transposed_;
        if(you_are_transposed != i_am_transposed)
        {
          ErrThrow("Transposed state of input does not fit transposed state of cell."
              "Rearrange input!");
        }
      }
    }

    /* ************************************************************************* */
    std::shared_ptr<TMatrix> ColoredBlockMatrix::create_block_shared_pointer(const TMatrix& block)
    {

      Output::print("Called base class copy and store");
      return std::make_shared<TMatrix>(block);
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
    bool ColoredBlockMatrix::does_modification_require_color_split(
        size_t& color_to_split,
        std::vector<grid_place_and_mode> row_column_transposed_tuples) const
    {
      // the part from here performs in O(m log m), when m is size of row_column_tuples

      // fill a list with the colors which will be affected
      // by the modification and sort it
      std::list<size_t> color_touches;
      for(auto it : row_column_transposed_tuples)
      {
        size_t color = cell_grid_[std::get<0>(it)][std::get<1>(it)].color_;
        color_touches.push_back( color );
      }
      color_touches.sort();

      //push back one more element as "stop" element
      color_touches.push_back(std::numeric_limits<size_t>::max());

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


      // great, no color split needed - we finished in O(m log m)
      return false;
    }

    /* ************************************************************************* */
    bool ColoredBlockMatrix::does_replace_require_color_merge(
        size_t& color_a, size_t& color_b,
        std::vector<grid_place_and_mode> row_column_transposed_tuples) const
    {
      // grab the first color as "fixed_color"
      size_t first_cell_row = std::get<0>(*row_column_transposed_tuples.begin());
      size_t first_cell_column = std::get<1>(*row_column_transposed_tuples.begin());
      size_t fixed_color = cell_grid_[first_cell_row][first_cell_column].color_;

      // find out if any other color is affected by the replacement
      for(auto it : row_column_transposed_tuples)
      {
        size_t cell_row = std::get<0>(it);
        size_t cell_column = std::get<1>(it);
        size_t cell_color = cell_grid_[cell_row][cell_column].color_;

        if (fixed_color != cell_color)
        { // found two different participating colors!
          color_a = fixed_color;
          color_b = cell_color;
          return true;
        }

      }
      // all participating cells belong to the same color class already
      color_a = fixed_color;
      color_b = fixed_color;
      return false;

    }


    /* ************************************************************************* */
    void ColoredBlockMatrix::find_first_appearance_of_color(
        size_t color_to_find, size_t& block_row , size_t& block_column) const
    {
      if( color_to_find >= get_n_colors() )
      {
        ErrThrow("That color does not exist in this ColoredBlockMatrix.");
      }

      for(size_t i = 0; i < n_cell_rows_ ; ++i)
      {
        for(size_t j = 0; j < n_cell_columns_ ; ++j)
        {
          if (cell_grid_[i][j].color_ == color_to_find)
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
    void ColoredBlockMatrix::find_first_appearance_of_color_and_mode(
        size_t color_to_find, bool find_transposed,
        size_t& block_row , size_t& block_column ) const
    {
      if( color_to_find >= get_n_colors() )
      {
        ErrThrow("That color does not exist in this ColoredBlockMatrix.");
      }

      for(size_t i = 0; i < n_cell_rows_ ; ++i)
      {
        for(size_t j = 0; j < n_cell_columns_ ; ++j)
        {
          if (cell_grid_[i][j].color_ == color_to_find
              && cell_grid_[i][j].is_transposed_ == find_transposed)
          {
            // color was found
            block_row = i;
            block_column = j;
            return;
          }
        }
      }
      throw std::runtime_error("Could not find a cell of that color and mode.");
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::get_next_cell_grid_index (
        size_t& block_row, size_t& block_column ) const
    {

      if( block_row > n_cell_rows_
          || block_column > n_cell_columns_)
      {
        throw std::runtime_error("Index is out of bounds!");
      }

      if( is_last_index_pair ( block_row, block_column ) )
      {
        throw std::logic_error("Next index pair after last required!");
      }

      size_t new_row = (block_column < n_cell_columns_ - 1) ? block_row : block_row + 1;
      size_t new_column = (block_row == new_row) ? block_column + 1 : 0;

      block_row = new_row;
      block_column = new_column;

    }

    /* ************************************************************************* */
    bool ColoredBlockMatrix::is_last_index_pair(size_t block_row, size_t block_column) const
    {
      if( block_row > n_cell_rows_
          || block_column > n_cell_columns_)
      {
        throw std::runtime_error("Index is out of bounds!");
      }

      if( block_row == n_cell_rows_ -1
          && block_column == n_cell_columns_ - 1)
      {
        return true;
      }
      return false;
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::mark_for_color_split(size_t color_to_mark,
                                                  std::vector<grid_place_and_mode>
    row_column_transposed_tuples)
    {
      for (auto it : row_column_transposed_tuples)
      {
        if (cell_grid_[std::get<0>(it)][std::get<1>(it)].color_ == color_to_mark)
        {
          //the cell is of the color to split and belongs to the given subset
          // - so mark it for split
          cell_grid_[std::get<0>(it)][std::get<1>(it)].re_color_flag_ =
              CellInfo::ReColoringFlag::SPLIT;
        }
      }
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::merge_colors(size_t color_a, size_t color_b)
    {
      //hold the color of the merged class
      size_t new_color = std::min(color_a, color_b);

      // hold the color which will get replaced
      size_t replaced_color = std::max(color_a, color_b);

      // a shared pointer to the matrix which the merged color class will share
      std::shared_ptr<TMatrix> shared_matrix;
      bool searching_shared_matrix = true;

      //traverse the entire cell grid
      for (size_t i = 0; i < n_cell_rows_; ++i)
      {
        for (size_t j = 0; j < n_cell_columns_; ++j)
        {
          CellInfo& currentCell = cell_grid_[i][j];


          if (currentCell.color_ == color_a || currentCell.color_ == color_b)
          {
            // found  matrix to be merged into new class
            if(searching_shared_matrix)
            {
              //found the first appeareance of the lower class, hold its block
              shared_matrix = currentCell.block_;
              searching_shared_matrix = false;
            } // end if searching shared matrix

            // replace matrix and color and update color count (it's even done when unnecessary)
            color_count_[currentCell.color_] -= 1;
            currentCell.color_ = new_color;
            currentCell.block_ = shared_matrix;
            color_count_[currentCell.color_] += 1;

          } //end: found cell of new merged class

          else if (currentCell.color_ > replaced_color)
          {
            // we found a class whose color has to be reduced by one
            color_count_[currentCell.color_] -= 1;
            currentCell.color_ -= 1;
            color_count_[currentCell.color_] += 1;
          }
        }
      } //end traversing entire matrix

      // remove the last color class - should be empty (0)!
      if(color_count_.back() != 0)
      {
        ErrThrow("Something's terribly wrong here - last color class is not empty!")
      }

      color_count_.pop_back();
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::replace_blocks(
        const TMatrix& new_block,
        std::vector<grid_place_and_mode> row_column_transpose_tuples)
    {
      // first of all check the input, modify if reparable or throw if not so.
      check_grid_fit(new_block, row_column_transpose_tuples);

      // check if the replacement requires color splits and if so perform them.
      size_t colorToSplit = std::numeric_limits<size_t>::max();
      while (does_modification_require_color_split(colorToSplit, row_column_transpose_tuples))
      {
        mark_for_color_split(colorToSplit, row_column_transpose_tuples);
        split_color(colorToSplit);
      }

      // check if the replacement requires color merges and if so perform them.
      size_t colorA = std::numeric_limits<size_t>::max();
      size_t colorB = std::numeric_limits<size_t>::max();
      while(does_replace_require_color_merge(colorA, colorB, row_column_transpose_tuples))
      {
        merge_colors(colorA, colorB);
      }

      // now everything works out - do the actual replacement!

      //wrap shared ptr around (correctly typed) copy of new_block
      std::shared_ptr<TMatrix> new_block_shared = create_block_shared_pointer(new_block);

      for (auto it : row_column_transpose_tuples)
      {
        size_t cell_row = std::get<0>(it);
        size_t cell_column = std::get<1>(it);
        bool transposed = std::get<2>(it);
        //copy assign new shared pointer
        cell_grid_[cell_row][cell_column].block_ = new_block_shared;
        // set the transposed state correctly
        cell_grid_[cell_row][cell_column].is_transposed_ = transposed;
      }
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::scale_blocks(
        double scaling_factor,
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

      // delegate the scaling to the TMatrices
      for (auto it: row_column_transpose_tuples)
      {
        size_t cell_row = std::get<0>(it);
        size_t cell_column = std::get<1>(it);
        CellInfo& current_cell = cell_grid_[cell_row][cell_column];
        size_t cell_color = current_cell.color_;

        if (cell_color >= searched_color)
        { // we found an untreated color
          current_cell.block_.get()->scale( scaling_factor );
          searched_color = cell_color + 1;
        }
      }
    }

    /* ************************************************************************* */
    void ColoredBlockMatrix::split_color(size_t color_to_split)
    {
      Output::print<2>("A color of this BlockMatrix had to be split."
          " Is that what you intended?");

      // a new color is going to appear - start with 0 count
      color_count_.push_back(0);

      size_t i = 0;
      size_t j = 0;

      find_first_appearance_of_color(color_to_split, i , j);

      //deep copy the matrix and make the shared_ptr temporarily responsible
      std::shared_ptr<TMatrix> matrix_copy = create_block_shared_pointer(*cell_grid_[i][j].block_.get());

      // Whether the first found cell of the color to split
      // has flag SPLIT or flag KEEP determines, which part of the color
      // class gets the new numbering - the one marked with value
      // equal to first_split maintains its current color number
      CellInfo::ReColoringFlag first_split(cell_grid_[i][j].re_color_flag_);

      //reset recoloring flag in this first found cell
      cell_grid_[i][j].re_color_flag_ = CellInfo::ReColoringFlag::KEEP;

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

        CellInfo& current_cell = cell_grid_[i][j];


        if(current_cell.color_ >= new_color) // a cell of a higher color
        {//found a cell of higher color
          if (searching_first_cell_of_new_class)
          {//this means we have not yet found the first element
            // of the second part of the color class to split
            // and thus still have to determine its number

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

          if (current_cell.re_color_flag_ == first_split)
          {//the element belongs to the first set and thus keeps its color
            //just reset the re color flags
            current_cell.re_color_flag_ = CellInfo::ReColoringFlag::KEEP;
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
            current_cell.re_color_flag_ = CellInfo::ReColoringFlag::KEEP;
          }
        } //end: cell of split color class
      } //endwhile

    }
