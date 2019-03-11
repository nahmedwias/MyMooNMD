#include "BlockFEMatrixPr.h"
#include "../include/General/Database.h"
#include "../include/Matrix/FEMatrix.h"
#include "../include/General/MooNMD_Io.h"


#include <limits>
#include <algorithm>

BlockFEMatrixPr::BlockFEMatrixPr(std::vector< const TFESpace2D* > spaces_rows, 
                    std::vector< const TFESpace2D* > spaces_cols)
: BlockFEMatrix()
{
  this->test_spaces_rowwise_ = spaces_rows;
  this->ansatz_spaces_columnwise_ = spaces_cols;
  
  Output::print<5>("BlockFEMatrix constructor");
  // class invariant: testspaces are not allowed to hold hanging nodes,
  // as the only kind of non-active dofs this class can handle is Dirichlet dofs
  for(auto sp : spaces_rows)
  {
    if (sp->GetN_SlaveDegrees() != 0 )
    {//there is  slave dofs, i.e. hanging nodes, quit
      ErrThrow("BlockFEMatrix cannot handle hanging nodes so far! "
          "Use spaces without hanging nodes.");
    }
  }
  
  for(auto sp : spaces_cols)
  {
    if (sp->GetN_SlaveDegrees() != 0 )
    {//there is  slave dofs, i.e. hanging nodes, quit
      ErrThrow("BlockFEMatrix cannot handle hanging nodes so far! "
          "Use spaces without hanging nodes.");
    }
  }
  

  //reset grid fittingly
  n_cell_rows_ = spaces_rows.size();
  n_cell_columns_ = spaces_cols.size();
  cell_grid_ = std::vector<std::vector<CellInfo>>(n_cell_rows_, std::vector<CellInfo>(n_cell_columns_));
  color_count_ = std::vector<size_t>(); //reset color counter vector

  //traverse the cell info grid from left to right, top to bottom
  // and fill in newly constructed, correctly dimensioned zero matrices as blocks
  for(size_t i = 0; i < n_cell_rows_ ; ++i )
  {
    //hold the testspace each cell in this row will have
#ifdef __2D__
    const TFESpace2D& testspace_of_row = *spaces_rows[i];
#elif __3D__
    const TFESpace3D& testspace_of_row = *spaces_rows[i];
#endif
    //hold the number of rows each cell in this row will have
    size_t n_rows_of_cell = testspace_of_row.GetN_DegreesOfFreedom();

    for(size_t j = 0; j < n_cell_columns_ ; ++j )
    {
      //hold the ansatzspace each cell in this column will have
#ifdef __2D__
      const TFESpace2D& ansatzspace_of_column = *spaces_cols[j];
#elif __3D__
      const TFESpace3D& ansatzspace_of_column = *spaces_cols[j];
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


BlockFEMatrixPr BlockFEMatrixPr::Projection_NSE2D(const TFESpace2D& velocity, const TFESpace2D& projection, const TFESpace2D& pressure)
{
  BlockFEMatrixPr my_matrix({&velocity}, {&projection, &projection});
  
  if( (velocity.GetCollection()) != projection.GetCollection() )
  {
    Output::print<1>(velocity.GetCollection()->GetN_Cells(), " ", 
                     projection.GetCollection()->GetN_Cells());
    ErrThrow("grid is different for the test and ansatz sapces");
  }
  
  // first and second column block
  my_matrix.replace_blocks(FEMatrix(&velocity, &projection), {{0,0}}, {false});
  my_matrix.replace_blocks(FEMatrix(&velocity, &projection), {{0,1}}, {false});
  
  return my_matrix;
}
