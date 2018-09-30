#include <BlockFEMatrix.h>
#include <FEMatrix.h>
#include <BlockVector.h>
#include <Database.h>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#endif

#include <limits>
#include <algorithm>

#ifdef __2D__
bool determine_need_for_pressure_row_correction(std::vector<const TFESpace2D*> spaces);
#else //__3D__
bool determine_need_for_pressure_row_correction(std::vector<const TFESpace3D*> spaces);
#endif

/* ************************************************************************* */
// IMPLEMENTATION OF PUBLIC METHODS
/* ************************************************************************* */

#ifdef __2D__
BlockFEMatrix::BlockFEMatrix(
    std::vector< const TFESpace2D*  > spaces) :
#else
    BlockFEMatrix::BlockFEMatrix(
        std::vector< const TFESpace3D*  > spaces) :
#endif
    BlockMatrix(), //base class object is default (empty) constructed
    test_spaces_rowwise_(spaces),
    ansatz_spaces_columnwise_(spaces)
{
  Output::print<5>("BlockFEMatrix constructor");
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
#else
    const TFESpace3D& testspace_of_row = *spaces[i];
#endif
    //hold the number of rows each cell in this row will have
    size_t n_rows_of_cell = testspace_of_row.GetN_DegreesOfFreedom();

    for(size_t j = 0; j < n_cell_columns_ ; ++j )
    {
      //hold the ansatzspace each cell in this column will have
#ifdef __2D__
      const TFESpace2D& ansatzspace_of_column = *spaces[j];
#else
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

  // TODO figure out whether pressure correction is needed
  use_pressure_projection_ = determine_need_for_pressure_row_correction(ansatz_spaces_columnwise_);

}

BlockFEMatrix::BlockFEMatrix() :
BlockMatrix(), test_spaces_rowwise_(),
ansatz_spaces_columnwise_()
{
}

//named constructors
#ifdef __2D__
/* ************************************************************************* */
BlockFEMatrix BlockFEMatrix::CD2D( const TFESpace2D& space )
{
  BlockFEMatrix my_matrix({&space});
  //replace block by a block with non-empty TStructure
  my_matrix.replace_blocks(FEMatrix(&space), {{0,0}} , {false});
  return my_matrix;
}

BlockFEMatrix BlockFEMatrix::NSE2D_Type1( const TFESpace2D& velocity,
                                         const TFESpace2D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &pressure});

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

BlockFEMatrix BlockFEMatrix::NSE2D_Type2(
    const TFESpace2D& velocity, const TFESpace2D& pressure)

{
  BlockFEMatrix my_matrix({&velocity, &velocity, &pressure});

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

BlockFEMatrix BlockFEMatrix::NSE2D_Type3(
    const TFESpace2D& velocity, const TFESpace2D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &pressure});

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

BlockFEMatrix BlockFEMatrix::NSE2D_Type4(
    const TFESpace2D& velocity, const TFESpace2D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &pressure});

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

BlockFEMatrix BlockFEMatrix::NSE2D_Type14(
    const TFESpace2D& velocity, const TFESpace2D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &pressure});

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

BlockFEMatrix BlockFEMatrix::Darcy2D( const TFESpace2D& velocity, const TFESpace2D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &pressure});

  //fill in the blocks with correct matrices constructed solely for them
  my_matrix.replace_blocks(FEMatrix(&velocity, &velocity), {{0,0}}, {false});
  my_matrix.replace_blocks(FEMatrix(&velocity, &pressure), {{0,1}}, {false});
  my_matrix.replace_blocks(FEMatrix(&pressure, &velocity), {{1,0}}, {false});
  my_matrix.replace_blocks(FEMatrix(&pressure, &pressure), {{1,1}}, {false});

  return my_matrix;
}

BlockFEMatrix BlockFEMatrix::Mass_NSE2D(const TFESpace2D& velocity)
{
  BlockFEMatrix my_matrix({&velocity, &velocity});
  
  my_matrix.replace_blocks(FEMatrix(&velocity, &velocity), {{0,0}, {1, 1}}, 
                           {false, false});

  return my_matrix;

}

BlockFEMatrix BlockFEMatrix::Mass_Matrix_NSE2D(const TFESpace2D& velocity, 
                              const TFESpace2D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &pressure});
  
  //create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity); //A blocks
  FEMatrix velo_velo_0_1(velo_velo_0_0);
  FEMatrix velo_velo_1_0(velo_velo_0_0);
  FEMatrix velo_velo_1_1(velo_velo_0_0); //
  
  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}}, {false});
  my_matrix.replace_blocks(velo_velo_0_1, {{0,1}}, {false});
  my_matrix.replace_blocks(velo_velo_1_0, {{1,0}}, {false});
  my_matrix.replace_blocks(velo_velo_1_1, {{1,1}}, {false});
  
  return my_matrix;
}

BlockFEMatrix BlockFEMatrix::Mass_NSE2D_Type1( const TFESpace2D& velocity,
                                         const TFESpace2D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &pressure});

  //create new blocks with correct structures filled with 0
  FEMatrix velo_velo(&velocity, &velocity); //A block
  FEMatrix velo_velo_zero(&velocity, &velocity, true); //velocity zero block

  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo, {{0,0},{1,1}}, {false, false});
  my_matrix.replace_blocks(velo_velo_zero, {{1,0},{0,1}}, {false, false});

  return my_matrix;
}

BlockFEMatrix BlockFEMatrix::Mass_NSE2D_Type2(
    const TFESpace2D& velocity, const TFESpace2D& pressure)

{
  BlockFEMatrix my_matrix({&velocity, &velocity, &pressure});

  //create new blocks with correct structures filled with 0
  FEMatrix velo_velo(&velocity, &velocity); //A block
  FEMatrix velo_velo_zero(&velocity, &velocity, true); //velocity zero block

  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo, {{0,0},{1,1}}, {false, false});
  my_matrix.replace_blocks(velo_velo_zero, {{1,0},{0,1}}, {false, false});

  return my_matrix;
}


BlockFEMatrix BlockFEMatrix::Mass_NSE2D_Type3(
    const TFESpace2D& velocity, const TFESpace2D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &pressure});

  //create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity); //A blocks
  FEMatrix velo_velo_0_1(velo_velo_0_0);
  FEMatrix velo_velo_1_0(velo_velo_0_0);
  FEMatrix velo_velo_1_1(velo_velo_0_0); //all copy constructed, share one TStructure

  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}}, {false});
  my_matrix.replace_blocks(velo_velo_0_1, {{0,1}}, {false});
  my_matrix.replace_blocks(velo_velo_1_0, {{1,0}}, {false});
  my_matrix.replace_blocks(velo_velo_1_1, {{1,1}}, {false});

  // B's and C stays as default initialized

  return my_matrix;
}

BlockFEMatrix BlockFEMatrix::Mass_NSE2D_Type4(
    const TFESpace2D& velocity, const TFESpace2D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &pressure});

  //create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity); //A blocks
  FEMatrix velo_velo_0_1(velo_velo_0_0);
  FEMatrix velo_velo_1_0(velo_velo_0_0);
  FEMatrix velo_velo_1_1(velo_velo_0_0); //all copy constructed, share one TStructure

  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}}, {false});
  my_matrix.replace_blocks(velo_velo_0_1, {{0,1}}, {false});
  my_matrix.replace_blocks(velo_velo_1_0, {{1,0}}, {false});
  my_matrix.replace_blocks(velo_velo_1_1, {{1,1}}, {false});

  //B's and C stays as default initialized

  return my_matrix;
}
#elif __3D__
//3D named constructors
BlockFEMatrix BlockFEMatrix::CD3D( const TFESpace3D& space )
{
  BlockFEMatrix my_matrix({&space});
  //replace block by a block with non-empty TStructure
  my_matrix.replace_blocks(FEMatrix(&space), {{0,0}} , {false});
  return my_matrix;
}

BlockFEMatrix BlockFEMatrix::Darcy3D(const TFESpace3D& velocity,
                                     const TFESpace3D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &pressure});

  //fill in the blocks with correct matrices constructed solely for them
  my_matrix.replace_blocks(FEMatrix(&velocity, &velocity), {{0,0}}, {false});
  my_matrix.replace_blocks(FEMatrix(&velocity, &pressure), {{0,1}}, {false});
  my_matrix.replace_blocks(FEMatrix(&pressure, &velocity), {{1,0}}, {false});
  my_matrix.replace_blocks(FEMatrix(&pressure, &pressure), {{1,1}}, {false});

  return my_matrix;
}

BlockFEMatrix BlockFEMatrix::NSE3D_Type1( const TFESpace3D& velocity, const TFESpace3D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &velocity, &pressure});

  // create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity);
  
  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}, {1,1}, {2,2}},
                                          {false, false, false});

  // zero velocity blocks
  FEMatrix velo_velo_zero(&velocity, &velocity, true);
  my_matrix.replace_blocks(velo_velo_zero, {{0,1}, {0,2}, {1,0},
                                            {1,2}, {2,0}, {2,1}},
                                            {false, false, false,
                                             false, false, false});


  FEMatrix pressure_velo_1(&pressure, &velocity);
  FEMatrix pressure_velo_2(pressure_velo_1);
  FEMatrix pressure_velo_3(pressure_velo_1); // copy constructed, shares TStructure!
 
  // fill in the pressure_velo blocks B_1 and B_1^T
  my_matrix.replace_blocks(pressure_velo_1, {{3,0}, {0,3} }, {false, true});
 
  // fill in the pressure_velo blocks B_2 and B_2^T
  my_matrix.replace_blocks(pressure_velo_2, {{3,1}, {1,3} }, {false, true});

  // fill in the pressure_velo blocks B_3 and B_3^T
  my_matrix.replace_blocks(pressure_velo_3, {{3,2}, {2,3} }, {false, true});  

  return my_matrix;
}

BlockFEMatrix BlockFEMatrix::NSE3D_Type2( const TFESpace3D& velocity, const TFESpace3D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &velocity, &pressure});
  // create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity);
  
  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}, {1,1}, {2,2}},
                                          {false, false, false});

  // zero velocity blocks
  FEMatrix velo_velo_zero(&velocity, &velocity, true);
  my_matrix.replace_blocks(velo_velo_zero, {{0,1}, {0,2}, {1,0},
                                            {1,2}, {2,0}, {2,1}},
                                            {false, false, false,
                                             false, false, false});
  
  FEMatrix pressure_velo_1(&pressure, &velocity);
  FEMatrix pressure_velo_2(pressure_velo_1); 
  FEMatrix pressure_velo_3(pressure_velo_1); // copy constructed, shares TStructure!

  FEMatrix velo_pressure_1(&velocity, &pressure);
  FEMatrix velo_pressure_2(velo_pressure_1); 
  FEMatrix velo_pressure_3(velo_pressure_1); // copy constructed, shares TStructure!


  // fill in the pressure_velo blocks B_1 and B_1^T
  my_matrix.replace_blocks(pressure_velo_1, {{3,0}}, {false});
  my_matrix.replace_blocks(velo_pressure_1, {{0,3}}, {false});  

  // fill in the pressure_velo blocks B_2 and B_2^T
  my_matrix.replace_blocks(pressure_velo_2, {{3,1}}, {false});
  my_matrix.replace_blocks(velo_pressure_2, {{1,3}}, {false});
  
  // fill in the pressure_velo blocks B_3 and B_3^T
  my_matrix.replace_blocks(pressure_velo_3, {{3,2}}, {false});
  my_matrix.replace_blocks(velo_pressure_3, {{2,3}}, {false});  

  return my_matrix;
}

BlockFEMatrix BlockFEMatrix::NSE3D_Type3( const TFESpace3D& velocity, const TFESpace3D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &velocity, &pressure});

  // create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity); //A block
  FEMatrix velo_velo_0_1(velo_velo_0_0);
  FEMatrix velo_velo_0_2(velo_velo_0_0);
  FEMatrix velo_velo_1_0(velo_velo_0_0);
  FEMatrix velo_velo_1_1(velo_velo_0_0);
  FEMatrix velo_velo_1_2(velo_velo_0_0);
  FEMatrix velo_velo_2_0(velo_velo_0_0);
  FEMatrix velo_velo_2_1(velo_velo_0_0);
  FEMatrix velo_velo_2_2(velo_velo_0_0);   //all copy constructed, share one TStructure
  
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}}, {false});
  my_matrix.replace_blocks(velo_velo_0_1, {{0,1}}, {false});
  my_matrix.replace_blocks(velo_velo_0_2, {{0,2}}, {false});
  my_matrix.replace_blocks(velo_velo_1_0, {{1,0}}, {false});
  my_matrix.replace_blocks(velo_velo_1_1, {{1,1}}, {false});
  my_matrix.replace_blocks(velo_velo_1_2, {{1,2}}, {false});
  my_matrix.replace_blocks(velo_velo_2_0, {{2,0}}, {false});
  my_matrix.replace_blocks(velo_velo_2_1, {{2,1}}, {false});
  my_matrix.replace_blocks(velo_velo_2_2, {{2,2}}, {false});

  FEMatrix pressure_velo_1(&pressure, &velocity);
  FEMatrix pressure_velo_2(pressure_velo_1); 
  FEMatrix pressure_velo_3(pressure_velo_1); // copy constructed, shares TStructure!
  
  // fill in the pressure_velo blocks B_1 and B_1^T
  my_matrix.replace_blocks(pressure_velo_1, {{3,0}, {0,3} }, {false, true});
 
  // fill in the pressure_velo blocks B_2 and B_2^T
  my_matrix.replace_blocks(pressure_velo_2, {{3,1}, {1,3} }, {false, true});

  // fill in the pressure_velo blocks B_3 and B_3^T
  my_matrix.replace_blocks(pressure_velo_3, {{3,2}, {2,3} }, {false, true});

  return my_matrix;
}

BlockFEMatrix BlockFEMatrix::NSE3D_Type4( const TFESpace3D& velocity, const TFESpace3D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &velocity, &pressure});

  // create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity); //A block
  FEMatrix velo_velo_0_1(velo_velo_0_0);
  FEMatrix velo_velo_0_2(velo_velo_0_0);
  FEMatrix velo_velo_1_0(velo_velo_0_0);
  FEMatrix velo_velo_1_1(velo_velo_0_0);
  FEMatrix velo_velo_1_2(velo_velo_0_0);
  FEMatrix velo_velo_2_0(velo_velo_0_0);
  FEMatrix velo_velo_2_1(velo_velo_0_0);
  FEMatrix velo_velo_2_2(velo_velo_0_0);   //all copy constructed, share one TStructure
  
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}}, {false});
  my_matrix.replace_blocks(velo_velo_0_1, {{0,1}}, {false});
  my_matrix.replace_blocks(velo_velo_0_2, {{0,2}}, {false});
  my_matrix.replace_blocks(velo_velo_1_0, {{1,0}}, {false});
  my_matrix.replace_blocks(velo_velo_1_1, {{1,1}}, {false});
  my_matrix.replace_blocks(velo_velo_1_2, {{1,2}}, {false});
  my_matrix.replace_blocks(velo_velo_2_0, {{2,0}}, {false});
  my_matrix.replace_blocks(velo_velo_2_1, {{2,1}}, {false});
  my_matrix.replace_blocks(velo_velo_2_2, {{2,2}}, {false});
  
  FEMatrix pressure_velo_1(&pressure, &velocity);
  FEMatrix pressure_velo_2(pressure_velo_1); 
  FEMatrix pressure_velo_3(pressure_velo_1); // copy constructed, shares TStructure!
  
  FEMatrix velo_pressure_1(&velocity, &pressure);
  FEMatrix velo_pressure_2(velo_pressure_1); 
  FEMatrix velo_pressure_3(velo_pressure_1); // copy constructed, shares TStructure!
  
  // fill in the pressure_velo blocks B_1 and B_1^T
  my_matrix.replace_blocks(pressure_velo_1, {{3,0}}, {false});
  my_matrix.replace_blocks(velo_pressure_1, {{0,3}}, {false});  
  
  // fill in the pressure_velo blocks B_2 and B_2^T
  my_matrix.replace_blocks(pressure_velo_2, {{3,1}}, {false});
  my_matrix.replace_blocks(velo_pressure_2, {{1,3}}, {false});
  
  // fill in the pressure_velo blocks B_3 and B_3^T
  my_matrix.replace_blocks(pressure_velo_3, {{3,2}}, {false});
  my_matrix.replace_blocks(velo_pressure_3, {{2,3}}, {false});  
  
  return my_matrix;

}

BlockFEMatrix BlockFEMatrix::NSE3D_Type14( const TFESpace3D& velocity, const TFESpace3D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &velocity, &pressure});

  // create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity); //A block
  FEMatrix velo_velo_0_1(velo_velo_0_0);
  FEMatrix velo_velo_0_2(velo_velo_0_0);
  FEMatrix velo_velo_1_0(velo_velo_0_0);
  FEMatrix velo_velo_1_1(velo_velo_0_0);
  FEMatrix velo_velo_1_2(velo_velo_0_0);
  FEMatrix velo_velo_2_0(velo_velo_0_0);
  FEMatrix velo_velo_2_1(velo_velo_0_0);
  FEMatrix velo_velo_2_2(velo_velo_0_0);   //all copy constructed, share one TStructure
  
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}}, {false});
  my_matrix.replace_blocks(velo_velo_0_1, {{0,1}}, {false});
  my_matrix.replace_blocks(velo_velo_0_2, {{0,2}}, {false});
  my_matrix.replace_blocks(velo_velo_1_0, {{1,0}}, {false});
  my_matrix.replace_blocks(velo_velo_1_1, {{1,1}}, {false});
  my_matrix.replace_blocks(velo_velo_1_2, {{1,2}}, {false});
  my_matrix.replace_blocks(velo_velo_2_0, {{2,0}}, {false});
  my_matrix.replace_blocks(velo_velo_2_1, {{2,1}}, {false});
  my_matrix.replace_blocks(velo_velo_2_2, {{2,2}}, {false});
  
  FEMatrix pressure_velo_1(&pressure, &velocity);
  FEMatrix pressure_velo_2(pressure_velo_1); 
  FEMatrix pressure_velo_3(pressure_velo_1); // copy constructed, shares TStructure!
  
  FEMatrix velo_pressure_1(&velocity, &pressure);
  FEMatrix velo_pressure_2(velo_pressure_1); 
  FEMatrix velo_pressure_3(velo_pressure_1); // copy constructed, shares TStructure!
  
  // fill in the pressure_velo blocks B_1 and B_1^T
  my_matrix.replace_blocks(pressure_velo_1, {{3,0}}, {false});
  my_matrix.replace_blocks(velo_pressure_1, {{0,3}}, {false});  
  
  // fill in the pressure_velo blocks B_2 and B_2^T
  my_matrix.replace_blocks(pressure_velo_2, {{3,1}}, {false});
  my_matrix.replace_blocks(velo_pressure_2, {{1,3}}, {false});
  
  // fill in the pressure_velo blocks B_3 and B_3^T
  my_matrix.replace_blocks(pressure_velo_3, {{3,2}}, {false});
  my_matrix.replace_blocks(velo_pressure_3, {{2,3}}, {false});  
  
  FEMatrix pressure_pressure(&pressure, &pressure);
  // fill in the pressure-pressure block
  my_matrix.replace_blocks(pressure_pressure, {{3,3}}, {false});
  
  return my_matrix;

}

BlockFEMatrix BlockFEMatrix::Mass_NSE3D(const TFESpace3D& velocity)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &velocity});

  my_matrix.replace_blocks(FEMatrix(&velocity, &velocity),
                           {{0,0}, {1,1}, {2,2}},
                           {false, false, false});

  return my_matrix;

}

BlockFEMatrix BlockFEMatrix::Mass_NSE3D_Type1( const TFESpace3D& velocity, const TFESpace3D& pressure)
{ 
  BlockFEMatrix my_matrix({&velocity, &velocity, &velocity, &pressure});
   
   // create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity);
  // fill in velocity-velocity blocks
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}, {1,1}, {2,2}},
                                          {false, false, false});
  // zero velocity blocks
  FEMatrix velo_velo_zero(&velocity, &velocity, true);
  // zero blocks
  my_matrix.replace_blocks(velo_velo_zero, {{0,1}, {0,2}, {1,0},
                                            {1,2}, {2,0}, {2,1}},
                                            {false, false, false,
                                             false, false, false});
  
  return my_matrix;
}

BlockFEMatrix BlockFEMatrix::Mass_NSE3D_Type2( const TFESpace3D& velocity, const TFESpace3D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &velocity, &pressure});
  // create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity);
  
  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}, {1,1}, {2,2}},
                                          {false, false, false});

    // zero velocity blocks
  FEMatrix velo_velo_zero(&velocity, &velocity, true);
  my_matrix.replace_blocks(velo_velo_zero, {{0,1}, {0,2}, {1,0},
                                            {1,2}, {2,0}, {2,1}},
                                            {false, false, false,
                                             false, false, false});

  return my_matrix;
}

BlockFEMatrix BlockFEMatrix::Mass_NSE3D_Type3( const TFESpace3D& velocity, const TFESpace3D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &velocity, &pressure});

  // create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity); //A block
  FEMatrix velo_velo_0_1(velo_velo_0_0);
  FEMatrix velo_velo_0_2(velo_velo_0_0);
  FEMatrix velo_velo_1_0(velo_velo_0_0);
  FEMatrix velo_velo_1_1(velo_velo_0_0);
  FEMatrix velo_velo_1_2(velo_velo_0_0);
  FEMatrix velo_velo_2_0(velo_velo_0_0);
  FEMatrix velo_velo_2_1(velo_velo_0_0);
  FEMatrix velo_velo_2_2(velo_velo_0_0);   //all copy constructed, share one TStructure
  
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}}, {false});
  my_matrix.replace_blocks(velo_velo_0_1, {{0,1}}, {false});
  my_matrix.replace_blocks(velo_velo_0_2, {{0,2}}, {false});
  my_matrix.replace_blocks(velo_velo_1_0, {{1,0}}, {false});
  my_matrix.replace_blocks(velo_velo_1_1, {{1,1}}, {false});
  my_matrix.replace_blocks(velo_velo_1_2, {{1,2}}, {false});
  my_matrix.replace_blocks(velo_velo_2_0, {{2,0}}, {false});
  my_matrix.replace_blocks(velo_velo_2_1, {{2,1}}, {false});
  my_matrix.replace_blocks(velo_velo_2_2, {{2,2}}, {false});
  
  return my_matrix;
}

BlockFEMatrix BlockFEMatrix::Mass_NSE3D_Type4( const TFESpace3D& velocity, const TFESpace3D& pressure)
{
  BlockFEMatrix my_matrix({&velocity, &velocity, &velocity, &pressure});

  // create new blocks with correct structures filled with 0
  FEMatrix velo_velo_0_0(&velocity, &velocity); //A block
  FEMatrix velo_velo_0_1(velo_velo_0_0);
  FEMatrix velo_velo_0_2(velo_velo_0_0);
  FEMatrix velo_velo_1_0(velo_velo_0_0);
  FEMatrix velo_velo_1_1(velo_velo_0_0);
  FEMatrix velo_velo_1_2(velo_velo_0_0);
  FEMatrix velo_velo_2_0(velo_velo_0_0);
  FEMatrix velo_velo_2_1(velo_velo_0_0);
  FEMatrix velo_velo_2_2(velo_velo_0_0);   //all copy constructed, share one TStructure
  
  my_matrix.replace_blocks(velo_velo_0_0, {{0,0}}, {false});
  my_matrix.replace_blocks(velo_velo_0_1, {{0,1}}, {false});
  my_matrix.replace_blocks(velo_velo_0_2, {{0,2}}, {false});
  my_matrix.replace_blocks(velo_velo_1_0, {{1,0}}, {false});
  my_matrix.replace_blocks(velo_velo_1_1, {{1,1}}, {false});
  my_matrix.replace_blocks(velo_velo_1_2, {{1,2}}, {false});
  my_matrix.replace_blocks(velo_velo_2_0, {{2,0}}, {false});
  my_matrix.replace_blocks(velo_velo_2_1, {{2,1}}, {false});
  my_matrix.replace_blocks(velo_velo_2_2, {{2,2}}, {false});
  
  return my_matrix;
}
#endif

/* ************************************************************************* */

void BlockFEMatrix::add_matrix_actives(
    const FEMatrix& summand, double factor,
    const std::vector<std::vector<size_t>>& cell_positions,
    const std::vector<bool>& transposed_states)
{
  std::vector<std::tuple<size_t, size_t, bool>> input_tuples(
      check_and_tupelize_vector_input(cell_positions, transposed_states));


  add_scaled_actives(summand, factor, input_tuples);
}
/* ************************************************************************* */

void BlockFEMatrix::apply(const BlockVector & x, BlockVector & y) const
{
  //reset all values in 'y' to 0 and delegate to apply_scaled_add
  y.reset();
  apply_scaled_add(x, y, 1.0);
}
/* ************************************************************************* */

void BlockFEMatrix::apply_scaled_add(const BlockVector & x,
                                            BlockVector & y, double a) const
{
  // first do the multiplication for active rows only
  apply_scaled_add_actives(x,y,a);

  //and then update the non-active rows.
  if (!TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE)
    y.addScaledNonActive(x,a);

  /// FIXME It is an awful hack to refer to INTERNAL_FULL_MATRIX_STRUCTURE.
  /// Sometimes TStructures are assembled as if there were no non-actives
  /// (especially when using algebraic flux correction). In that case, the method
  /// apply_scaled_add_actives will think, that all rows of a stored FEMatrix
  /// are active, if then additionally addScaledNonActive is called this will
  /// produce wrong results. Two big design problems get apparent here
  ///
  /// 1) Referring to global scope for such a localized issue is a no-no.
  /// 2) The handling of non-actives is still a mess and must be fixed
  ///    in the whole programm, in the process of restructuring the Assembling.
  ///    Especially there may not be multiple classes which care for the
  ///    non-active issue (as now FESpace, TStructure AND BlockFEMatrix).
  ///    I'd like to entirely remove the information from TStructure and FEMatrix,
  ///    and place it in FESpace, while leaving the handling to the BlockFEMatrix.
  ///

}

/* ************************************************************************* */

void BlockFEMatrix::apply_scaled_add_actives(const BlockVector & x, BlockVector & y,
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
}

/* ************************************************************************* */
void BlockFEMatrix::apply_scaled_submatrix(const BlockVector & x, BlockVector & y,
                                        size_t sub_row, size_t sub_col,
                                        double a) const
{ 
  //check if the vectors fit, if not so the program throws an error
  //FIXME commented due to Mass_NSE2D matrix does not fit here  
  //FIXME Also the "=" sign in the check statement have been removed 
  //due to the Mass_NSE2D matrices multiplication: 
  // This reduces the coding in the main class as well easy for time stepping schemes
  
  // check_vector_fits_pre_image(x); 
  // check_vector_fits_image(y);
  
  if(sub_row > this->n_cell_rows_
      || sub_col > this->n_cell_columns_)
  {
    ErrThrow("Submatrix is not correctly specified!", sub_row, " ", this->n_cell_rows_,
             " ", sub_col, "  " , this->n_cell_columns_ );
  }

  const double * xv = x.get_entries(); // array of values in x
  double * yv = y.get_entries(); // array of values in y
  size_t row_offset = 0;
  // n_rows, n_cols are the number of cell rows/columns
  for(size_t i = 0; i < sub_row; ++i)
  {
    int col_offset = 0;
    for(size_t j = 0; j < sub_col; j++)
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
}

/* ************************************************************************* */
void BlockFEMatrix::handle_discovery_of_vector_actives(const int nActive, 
                                                     const int spaceNumber) const
{
  if(nActive != this->get_row_space(spaceNumber).GetN_ActiveDegrees())
  {
    ErrThrow("Number of actives in this vector's block ", spaceNumber,
             " does not fit the number of non-actives in the corresponding test space",
             " of this matrix!");
  }
}
/* ************************************************************************* */

void BlockFEMatrix::check_pointer_types()
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

void BlockFEMatrix::check_vector_fits_image(const BlockVector& b) const
{
  //let the base class figure out block dimensions
  BlockMatrix::check_vector_fits_image(b);

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

void BlockFEMatrix::check_vector_fits_pre_image(const BlockVector& x) const
{
  //let the base class figure out block dimensions
  BlockMatrix::check_vector_fits_pre_image(x);

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

std::shared_ptr<const FEMatrix> BlockFEMatrix::get_block(
    size_t cell_row, size_t cell_col,  bool& is_transposed) const
{
  //find out the transposed state
  is_transposed = cell_grid_.at(cell_row).at(cell_col).is_transposed_;

  //cast const and FEMatrix (range check is done via "at")
  std::shared_ptr<const FEMatrix> shared
  = std::dynamic_pointer_cast<const FEMatrix>(cell_grid_.at(cell_row).at(cell_col).block_);
  return shared;
}


/* ************************************************************************* */

std::vector<std::shared_ptr<const FEMatrix>> BlockFEMatrix::get_blocks() const
{
  std::vector<std::shared_ptr<const FEMatrix>> block_ptrs;

  for(size_t i =0 ; i < n_cell_columns_ ; ++i)
  {
    for(size_t j =0 ; j < n_cell_rows_ ; ++j)
    {
      //juggle the pointers around until we can store it the way we want...
      std::shared_ptr<const FEMatrix> shared //cast const and FEMatrix
      = std::dynamic_pointer_cast<const FEMatrix>(cell_grid_[i][j].block_);

      // push it back
      block_ptrs.push_back(shared);
    }
  }

  return block_ptrs;
}

/* ************************************************************************* */

std::vector<std::shared_ptr<FEMatrix>> BlockFEMatrix::get_blocks_TERRIBLY_UNSAFE()
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

std::vector<std::shared_ptr<FEMatrix>> BlockFEMatrix::get_blocks_uniquely(
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

std::vector<std::shared_ptr<FEMatrix>> BlockFEMatrix::get_blocks_uniquely(
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

std::shared_ptr<TMatrix> BlockFEMatrix::get_combined_matrix() const
{
  //delegate the business far down
  return get_combined_submatrix({0,0},{n_cell_rows_-1, n_cell_columns_-1});
}
/* ************************************************************************* */

std::shared_ptr<TMatrix> BlockFEMatrix::get_combined_submatrix(
    std::pair<size_t,size_t> upper_left,
    std::pair<size_t,size_t> lower_right) const
{
  //let base class do as much work as possible
  std::shared_ptr<TMatrix> sub_cmat =
      this->BlockMatrix::get_combined_submatrix(upper_left, lower_right);

  size_t r_first = upper_left.first;
  size_t r_last  = lower_right.first;
  size_t c_first = upper_left.second;
  size_t c_last  = lower_right.second;


  //check for empty diagonal blocks with non-actives TODO Is there a need to fix this?
  for(size_t diag = 0; diag < n_cell_rows_; ++diag)
  {
    bool diag_in =     (r_first <= diag) && (r_last >= diag)
                    && (c_first <= diag) && (c_last >= diag);

    if(!diag_in)
      continue; // this diag block is not requested anyway
    if(this->cell_grid_[diag][diag].block_->GetN_Entries()==0)
    {//zero-block on diagonal
      if(get_n_row_actives(diag) != get_n_rows_in_cell(diag, diag))
      {// This is trouble, because the baseclass will not place any entries at all
        // into the combined matrix where this empty block stands.
        // This means, that there will be no entries on the diagonal for the
        // BlockFEMatrix to put ones to - will result in 0 rows!
        Output::print("Warning! Trying to get combined matrix of a BlockFEMatrix "
            "with a zero block on a diagonal with test-space non-actives.");
      }
    }
  }

  // A: CORRECTIONS DUE TO DIRICHLET ROWS
  //get pointers in the matrix
  const int* rowptr = sub_cmat->GetRowPtr();
  int* kcolptr = sub_cmat->GetKCol();
  double* entries = sub_cmat->GetEntries();

  size_t row_offset = 0;

  //loop through all relevant cell rows
  for(size_t i = r_first; i <= r_last ;++i)
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

  //B: CORRECTIONS DUE TO PRESSURE PROJECTION
  if(use_pressure_projection_)
  {
#ifdef __2D__
    size_t dim = 2;
#endif
#ifdef __3D__
    size_t dim = 3;
#endif
    if( use_pressure_projection_ && n_cell_rows_ == dim + 1
        && n_cell_columns_ == dim + 1 && r_first <= dim && r_last >= dim)
    {
#ifdef _MPI
      int my_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

      // let 0 send a ping to all other processes, to determine which pressure
      // row to eliminate
      int verb = Output::getVerbosity();
      Output::suppressAll();
      int p_dof = get_communicators().back()->dof_ping(0,0);
      Output::setVerbosity(verb);

      if(my_rank==0)
      {
#endif
        // the relevant block row is contained
        // number of velocity dofs
        int n_rows = this->get_blocks().at(0)->GetN_Rows();

        // find the first row of the third block column
        size_t skipped_block_rows = r_first;
        int begin = rowptr[(dim - skipped_block_rows )*n_rows];
        int end   = rowptr[(dim - skipped_block_rows)*n_rows+1]; //now we have the row

        size_t skipped_block_cols = c_first;
        int future_diagonal = end -1;
        int diagonal_relative = (dim-skipped_block_cols)*n_rows;
        for(int j = begin; j<end;++j)
        {
          entries[j]=0;
          if(kcolptr[j] >= (int)dim*diagonal_relative && future_diagonal == end-1)
            future_diagonal = j;
        }
        if(diagonal_relative < sub_cmat->GetN_Columns() && begin != end)
        {
          entries[future_diagonal] = 1;
          kcolptr[future_diagonal] = diagonal_relative;
        }
        else if(begin == end)
        {
          // the matrix row is empty, therefore we have to change the matrix
          // structure to include a value here.
          // This typically happens if you try to get only the C-block in a 
          // Navier-Stokes problem, which is (usually) an all-zero block with 
          // no entries at all.
          if(sub_cmat->GetN_Entries() == 0)
          {
            // there are no entries at all: create a matrix of the same size
            // with a single entry
            // new_rows array is one everywhere except at the beginning, this 
            // means there is exactly one entry in the first row.
            std::vector<int> new_rows(sub_cmat->GetN_Rows()+1, 1);
            new_rows[0] = 0;
            // that entry has index 0, so the entry is on the diagonal.
            std::vector<int> new_cols(1, 0); 
            // new structure
            auto s = std::make_shared<TStructure>(sub_cmat->GetN_Rows(),
                                                  sub_cmat->GetN_Columns(),
                                                  1, &new_cols[0], 
                                                  &new_rows[0]);
            sub_cmat.reset(new TMatrix(s));
            sub_cmat->GetEntries()[0] = 1; // set diagonal entry
          }
          else
          {
            // this should not happen.
            ErrThrow("Unable to return a pressure-pressure matrix C, which has "
                     "has no entries in the first row, but in other rows.");
          }
        }
#ifdef _MPI
      }
      else if (p_dof != -1) //this rank knows the affected dof and has to set its row to 0
      {//determine the row to be swept
        int n_rows_before = 0;
        for(size_t br = r_first; br < n_cell_rows_-1; ++br)
          n_rows_before += this->get_row_space(br).GetN_DegreesOfFreedom();
        int row_glob = n_rows_before + p_dof;
        auto b_row = sub_cmat->get_row_array().at(row_glob);
        auto e_row = sub_cmat->GetRowPtr()[row_glob+1];
        auto entries = sub_cmat->GetEntries();
        for(int i = b_row; i < e_row ; ++i)
        {
          entries[i] = 0; //nuke the entries.
        }
      }
#endif
    }
  }


  // Remove all zero entries from the structure and the entries array
  // TODO doing this here is very slow and it should be changed,
  // but removing zeroes is important for the interface with direct solvers.
  sub_cmat->remove_zeros(0);

  return sub_cmat;

}

/* ************************************************************************* */

#ifdef _MPI
/// Return a list of the FE communicators belonging to the FESpaces of
/// the rows/columns.
std::vector<const TParFECommunicator3D*> BlockFEMatrix::get_communicators() const
{
  std::vector<const TParFECommunicator3D*> comms;
  for(auto sp : test_spaces_rowwise_ )
  {
    comms.push_back(&sp->get_communicator());
  }
  return comms;
}
#endif

/* ************************************************************************* */

BlockFEMatrix BlockFEMatrix::get_sub_blockfematrix(size_t first, size_t last) const
{
  //check input
  if(first >= n_cell_rows_ || last >= n_cell_rows_ )
  {
    ErrThrow("Out of bounds in get_sub_blockfematrix: ",
             first, ", ", last, ", ", n_cell_rows_);
  }
  // step 1: construct a new block fematrix with correct fe spaces
#ifdef __2D__
    std::vector< const TFESpace2D*  > spaces;
#elif __3D__
    std::vector< const TFESpace3D*  > spaces;
#endif
  for( size_t sp = first; sp < last + 1; ++sp)
  {
    spaces.push_back(this->test_spaces_rowwise_.at(sp));
  }
  BlockFEMatrix sub_matrix(spaces); //construct empty blockfematrix!

  // step 2: fill in the blocks - maintaining the correct coloring
  // TODO this is copy-paste identical to corresp. part in
  // BlockMatrix::get_subblockmatrix - put in private method!
  std::vector< int > known_colors;
  std::vector< std::shared_ptr<TMatrix> > known_mats; //actually FEMatrices...

  for( size_t r = first; r < last + 1; ++r)
  {
    for( size_t c = first; c < last + 1; ++c)
    {

      int old_color = this->cell_grid_[r][c].color_;
      auto known = std::find(known_colors.begin(), known_colors.end(), old_color);

      if(known == known_colors.end())
      {//case: this color appears for the first time
        known_colors.push_back(old_color);

        //this is actually a shared ptr to FEMatrix...
        std::shared_ptr<TMatrix> new_mat =
            create_block_shared_pointer(*cell_grid_[r][c].block_.get());

        known_mats.push_back(new_mat);
        //reset known
        known = known_colors.end()-1;
      }

      size_t new_color = std::distance(known_colors.begin(), known); //position in vector
      size_t transp = this->cell_grid_[r][c].is_transposed_;
      size_t new_r = r-first;       // force element and color
      size_t new_c = c-first;       // into the new sub matrix' cell grid
      size_t old_color_sub =  sub_matrix.cell_grid_[new_r][new_c].color_;
      sub_matrix.cell_grid_[new_r][new_c].color_ = new_color;
      sub_matrix.cell_grid_[new_r][new_c].is_transposed_ = transp;
      sub_matrix.cell_grid_[new_r][new_c].block_ = known_mats.at(new_color);
      //color count
      ++sub_matrix.color_count_.at(new_color);
      --sub_matrix.color_count_.at(old_color_sub);
    }
  }

  //tidy up the color count vector
  sub_matrix.color_count_.erase(
      std::remove( sub_matrix.color_count_.begin(),
                   sub_matrix.color_count_.end(), 0),
                   sub_matrix.color_count_.end()
  );

  return sub_matrix;
}

/* ************************************************************************* */
size_t BlockFEMatrix::get_n_column_actives(size_t cell_column) const
{
  return ansatz_spaces_columnwise_.at(cell_column)->GetN_ActiveDegrees();
}
/* ************************************************************************* */

size_t BlockFEMatrix::get_n_row_actives(size_t cell_row) const
{
  return test_spaces_rowwise_.at(cell_row)->GetN_ActiveDegrees();
}

/* ************************************************************************* */
void BlockFEMatrix::print_matrix_info(std::string name) const
{
  //Gather some information to be printed.
  int n_spaces_row = this->get_n_cell_rows();
  int n_spaces_col = this->get_n_cell_columns();

  int my_rank = 0;
  int n_dof = this->get_n_total_rows();

#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  n_dof = 0;
  for(auto c : this->get_communicators())
  {
    n_dof += c->get_n_global_dof();
  }
#endif

  if(my_rank == 0)
  {
    Output::info("BlockFEMatrix", name, " (distributed)");
    Output::dash(n_spaces_row, " x ", n_spaces_col, " block structure");
    Output::dash(n_dof, " total d.o.f (across all processors)");
  }

}

/* ************************************************************************* */

void BlockFEMatrix::replace_blocks(
    const TMatrix& new_block,
    const std::vector<std::vector<size_t>>& cell_positions,
    const std::vector<bool>& transposed_states)
{
  ErrThrow("Don't try to put a TMatrix into an FEMatrix!");
}
/* ************************************************************************* */

void BlockFEMatrix::replace_blocks(
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
            "block's ansatz space at (",cell_row,",",cell_column,")");
      }
      if (grid_ansatz_space != block_test_space)
      {
        ErrThrow("Grid ansatz space does not match transposed "
            "block's test space at (",cell_row,",",cell_column,")");
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
  BlockMatrix::replace_blocks(
      new_block, cell_positions, transposed_states);

}
/* ************************************************************************* */

void BlockFEMatrix::scale_blocks_actives(
    double factor,
    const std::vector<std::vector<size_t>>& cell_positions )
{
  std::vector<bool> transposed_states (cell_positions.size(), false);
  std::vector<std::tuple<size_t, size_t, bool>> input_tuples(
      check_and_tupelize_vector_input(cell_positions, transposed_states));

  scale_blocks_actives(factor, input_tuples);
}

/* ************************************************************************* */
void BlockFEMatrix::enable_pressure_projection() const
{
  this->use_pressure_projection_ = true;
}

/* ************************************************************************* */
void BlockFEMatrix::disable_pressure_projection() const
{
  this->use_pressure_projection_ = false;
}

/* ************************************************************************* */
bool BlockFEMatrix::pressure_projection_enabled() const
{
  return this->use_pressure_projection_;
}

/* ************************************************************************* */

/* ************************************************************************* */
// IMPLEMENTATION OF SPECIAL MEMBER FUNCTION(S)
/* ************************************************************************* */

BlockFEMatrix::BlockFEMatrix(const BlockFEMatrix& other)
: BlockMatrix::BlockMatrix(other),
  test_spaces_rowwise_(other.test_spaces_rowwise_),
  ansatz_spaces_columnwise_(other.ansatz_spaces_columnwise_)
{
  Output::print<5>("BlockFEMatrix copy constructor!");

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

  this->use_pressure_projection_ = other.use_pressure_projection_;
}

BlockFEMatrix::BlockFEMatrix(BlockFEMatrix&& other)
: BlockMatrix::BlockMatrix(std::move(other)), //base class move
  test_spaces_rowwise_(std::move(other.test_spaces_rowwise_)),
  ansatz_spaces_columnwise_(std::move(other.ansatz_spaces_columnwise_))
{
  Output::print<5>("BlockFEMatrix move constructor!");

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
  this->use_pressure_projection_ = other.use_pressure_projection_;
}

/* ************************************************************************* */
void swap(BlockFEMatrix& first, BlockFEMatrix& second)
{
  Output::print<5>("BlockFEMatrix swap!");

  std::swap(first.n_cell_columns_, second.n_cell_columns_);
  std::swap(first.n_cell_rows_, second.n_cell_rows_);
  std::swap(first.cell_grid_, second.cell_grid_);
  std::swap(first.color_count_, second.color_count_);
  std::swap(first.ansatz_spaces_columnwise_, second.ansatz_spaces_columnwise_);
  std::swap(first.test_spaces_rowwise_, second.test_spaces_rowwise_);
  std::swap(first.use_pressure_projection_, second.use_pressure_projection_);
}

/* ************************************************************************* */
BlockFEMatrix& BlockFEMatrix::operator=(BlockFEMatrix other)
{
  Output::print<5>("BlockFEMatrix copy assignment!");

  //do a swap with the copy constructed object "other"
  swap(*this, other);

  return *this;
}

/* ************************************************************************* */
#ifdef __2D__
const TFESpace2D& BlockFEMatrix::get_test_space(size_t cell_row, size_t cell_column) const
#elif __3D__
const TFESpace3D& BlockFEMatrix::get_test_space(size_t cell_row, size_t cell_column) const
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
const TFESpace2D& BlockFEMatrix::get_ansatz_space(size_t cell_row, size_t cell_column) const
#elif __3D__
const TFESpace3D& BlockFEMatrix::get_ansatz_space(size_t cell_row, size_t cell_column) const
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
void BlockFEMatrix::add_scaled_actives(
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

std::shared_ptr<TMatrix> BlockFEMatrix::create_block_shared_pointer(const TMatrix& block) const
{
  try
  { //try to cast the given TMatrix to an FEMatrix and make an FEMatrix copy of it
    const FEMatrix& fe_block_ref = dynamic_cast<const FEMatrix&>(block);
    std::shared_ptr<TMatrix> fe_block_ptr= std::make_shared<FEMatrix>(fe_block_ref);
    return fe_block_ptr;
  }
  catch (std::bad_cast e)
  {//cast did not work!
    ErrThrow("TMatrix given. Make sure to fill a BlockFEMatrix only with FEMatrices!");
  }
}
/* ************************************************************************* */

void BlockFEMatrix::scale_blocks_actives( double scaling_factor,
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
void BlockFEMatrix::add_blockfe_matrix(const BlockFEMatrix& Matrix, double factor)
{
  if(this->test_spaces_rowwise_ != Matrix.test_spaces_rowwise_)
  {
    ErrThrow("test spaces are not same");
  }
  if(this->ansatz_spaces_columnwise_ != Matrix.ansatz_spaces_columnwise_)
  {
    ErrThrow("ansatz spaces are not same");
  }
  
  size_t colors = Matrix.color_count_.size();
  
  for(size_t c = 0; c<colors; ++c)
  {
    FEMatrix* matrix_to_add;
    std::vector<std::vector<size_t>> blocks_to_add_to;
    std::vector<bool> is_transposed;
    for(size_t br=0; br<Matrix.n_cell_rows_; ++br)
    {
      for(size_t bc=0; bc<Matrix.n_cell_columns_; ++ bc)
      {
        if(Matrix.cell_grid_[br][bc].color_ == c)
        {
          matrix_to_add = dynamic_cast<FEMatrix*>(Matrix.cell_grid_[br][bc].block_.get());          
          blocks_to_add_to.push_back({static_cast<unsigned long>(br),static_cast<unsigned long>(bc)});
          is_transposed.push_back(Matrix.cell_grid_[br][bc].is_transposed_);
        }
      }
    }
    if(matrix_to_add->GetN_Entries() == 0)
      continue;
    this->add_matrix_actives(*matrix_to_add, factor, blocks_to_add_to,is_transposed);
  }
}
/* ************************************************************************* */

// Implementation of a helper method, which is used to figure out,
// whether this is an enclosed flow problem, and hence pressure correction is
// needed.
#ifdef __2D__
bool determine_need_for_pressure_row_correction(std::vector<const TFESpace2D*> spaces)
#else //__3D__
bool determine_need_for_pressure_row_correction(std::vector<const TFESpace3D*> spaces)
#endif
{
  int my_rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  // if both theses conditions are fulfilled, (enclosed flow problem)
  // the matrix requires pressure row correction
  bool is_saddle_point_problem = false;
  bool is_enclosed_flow = false;

  // CHECK IF SADDLE POINT PROBLEM
  //determine which space is the last in the sequence of equals
  #ifdef __2D__
  auto p = std::adjacent_find(spaces.begin(), spaces.end(), std::not_equal_to<const TFESpace2D*>());
  auto first_space = spaces.front();
  #else //__3D__
  auto p = std::adjacent_find(spaces.begin(), spaces.end(), std::not_equal_to<const TFESpace3D*>());
  auto first_space = spaces.front();
  #endif

  auto penultimate_space = spaces.end() - 2;
  if(p == penultimate_space)
  {//the last space in the sequence of equals is the penultimate space
    is_saddle_point_problem = true;
  }

  if(is_saddle_point_problem)
  {
    // CHECK IF ENCLOSED FLOW
    int n_velo_neumann = first_space->get_n_neumann_dof();
    int n_velo_robin = first_space->get_n_robin_dof();
    /// ... add further types of bdry conditions to check?

    if(n_velo_neumann == 0 && n_velo_robin == 0)
      is_enclosed_flow = true;

#ifdef _MPI
    bool sbuf = is_enclosed_flow;
    bool rbuf = false;
    MPI_Allreduce(&sbuf, &rbuf, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
    //rbuf will contain "true", if any of the processes held "true" before
    is_enclosed_flow = rbuf;
#endif

  }

  bool needs_prc = is_saddle_point_problem && is_enclosed_flow;
  
  if(TDatabase::ParamDB->RE_NR == 180 || TDatabase::ParamDB->RE_NR == 395.
    || TDatabase::ParamDB->RE_NR == 280000 )
  {
    needs_prc = false;
  }


  if(needs_prc && my_rank == 0)
      Output::info("Pressure Projection","BlockFEMatrix identified as enclosed-flow saddle point matrix. "
                   "Pressure projection enabled.");

  //TODO: Fix me please
  if(TDatabase::ParamDB->RE_NR == 180. || TDatabase::ParamDB->RE_NR == 395.)
    return false;
  else
    return needs_prc;
  

}
