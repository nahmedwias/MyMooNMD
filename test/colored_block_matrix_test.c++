/*
 * Unit testing of the colored block matrix.
 *
 *  Created on: Dec 10, 2015
 *      Author: bartsch
 */

#include <BlockMatrix.h>

#include <vector>
#include <tuple>
#include <MooNMD_Io.h>

int main(int argc, char* argv[])
{
  {
  //check correct splitting and merging of colors when replacing blocks
  ColoredBlockMatrix myMatrix({5 , 5 , 5}, {5, 5 , 5});

  TMatrix placeHolderA(5,5);
  auto t_1 = std::make_tuple(0,0,false);
  auto t_2 = std::make_tuple(0,1,false);
  auto t_3 = std::make_tuple(1,0,false);
  auto t_4 = std::make_tuple(1,1,false);
  std::vector<std::tuple<size_t, size_t, bool>> vecA = {t_1, t_2, t_3, t_4};

  TMatrix placeHolderB(5,5);
  std::vector<std::tuple<size_t, size_t, bool>> vecB = {t_1, t_4};


  myMatrix.replace_blocks(placeHolderA,
                          {{0,0}, {0,1}, {1,0}, {1,1} },
                          {false,false,false,false});
  myMatrix.check_coloring();

  myMatrix.replace_blocks(placeHolderB,
                          { {0,0}, {1,1} },
                          { false, false } );
  myMatrix.check_coloring();

  myMatrix.replace_blocks(placeHolderA,
                          {{0,0}, {0,1}, {1,0}, {1,1} },
                          {false,false,false,false});
  myMatrix.check_coloring();

  }

  {
  //check if it works with tranposition
  ColoredBlockMatrix myBlockMatrix({3,7},{3,7});
  myBlockMatrix.check_coloring();

  TMatrix placeHolderA(7,3);

  myBlockMatrix.replace_blocks(placeHolderA,
                               {{0,1}, {1,0}},
                               {true, false}  );
  myBlockMatrix.print_and_check("insert of transp");
  myBlockMatrix.print_coloring_pattern("insert of transp", true);

  }

  //check blockwise addition
  {
    /**
     * A = [ 1,  0,  0,  0  ;
     *       0,  1,  0,  0  ;
     *       0,  0,  1,  0  ;
     *       2,  0,  0,  1 ]
     */

    //Matrix A
    int * RowA = new int[5];
    int * ColA = new int[5];
    ColA[0] = 0; ColA[1] = 1; ColA[2] =  2; ColA[3] = 0; ColA[4] = 3;
    RowA[0] = 0; RowA[1] = 1; RowA[2] = 2; RowA[3] = 3; RowA[4] = 5;
    std::shared_ptr<TStructure> structureA(new TStructure(4, 5, ColA, RowA));
    TMatrix matA(structureA);
    matA.GetEntries()[0] = 1;
    matA.GetEntries()[1] = 1;
    matA.GetEntries()[2] = 1;
    matA.GetEntries()[3] = 1;
    matA.GetEntries()[4] = 1;

    ColoredBlockMatrix myMatrix({4 , 2 , 4}, {4 , 2, 4});

    myMatrix.replace_blocks(
        matA,
        {{0,0}, {0,2}, {2,0}, {2,2}},
        {false, false, true, false} );

    myMatrix.print_and_check("after replace");
    myMatrix.print_coloring_pattern("myMatrix", true);
    //myMatrix.get_combined_matrix()->PrintFull();

    //check copy constructor
    ColoredBlockMatrix yourMatrix(myMatrix);
    yourMatrix.print_and_check("yourMatrix");
    yourMatrix.print_coloring_pattern("yourMatrix", true);

    TMatrix summand(matA);
    summand.GetEntries()[0] = 1;
    summand.GetEntries()[1] = 1;
    summand.GetEntries()[2] = 1;
    summand.GetEntries()[3] = 1;
    summand.GetEntries()[4] = 1;

    myMatrix.add_unscaled_matrix(summand,
                                  {{2,0}, {2,2}},
                                  {true, false}  );

    myMatrix.print_and_check("addition to replace");
    myMatrix.get_combined_matrix()->PrintFull();

    myMatrix.add_unscaled_matrix(
        summand ,
        { {2,0}, {2,2} },
        {true, false} );

    myMatrix.print_and_check("addition transpose");
    myMatrix.get_combined_matrix()->PrintFull();

    std::vector<std::vector<size_t>> positions = {{0,0},{2,0}};
    std::vector<bool> transp_states = {false,true};

    myMatrix.add_matrix(
        summand, 1.0, positions, transp_states);


    myMatrix.print_and_check("addition 2");
    myMatrix.get_combined_matrix()->PrintFull();

    //check copy assignment
    ColoredBlockMatrix yourAssignedMatrix({0},{0});
    yourAssignedMatrix = myMatrix;
    yourAssignedMatrix.print_and_check("yourAssignedMatrix");
    yourAssignedMatrix.print_coloring_pattern("yourAssignedMatrix", true);


    delete[] RowA;
    delete[] ColA;
  }

  Output::print("Test program finished.");

}


