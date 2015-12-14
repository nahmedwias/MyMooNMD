/*
 * Unit testing of the colored block matrix.
 *
 *  Created on: Dec 10, 2015
 *      Author: bartsch
 */

#include <ColoredBlockMatrix.h>
#include <vector>
#include <tuple>
#include <MooNMD_Io.h>

int main(int argc, char* argv[])
{
  {
  //check correct splitting and mergingf of colors when replacing blocks
  std::vector<size_t> cell_row_numbers = {5 , 5 , 2};
  std::vector<size_t> cell_column_numbers = {5, 5 , 2, 3};

  ColoredBlockMatrix myMatrix(cell_row_numbers, cell_column_numbers);

  TMatrix placeHolderA(5,5);
  auto t_1 = std::make_tuple(0,0,false);
  auto t_2 = std::make_tuple(0,1,false);
  auto t_3 = std::make_tuple(1,0,false);
  auto t_4 = std::make_tuple(1,1,false);
  std::vector<std::tuple<size_t, size_t, bool>> vecA = {t_1, t_2, t_3, t_4};

  TMatrix placeHolderB(5,5);
  std::vector<std::tuple<size_t, size_t, bool>> vecB = {t_1, t_4};

  myMatrix.replace_blocks(placeHolderA, vecA);
  myMatrix.check_coloring();

  myMatrix.replace_blocks(placeHolderB, vecB);
  myMatrix.check_coloring();

  myMatrix.replace_blocks(placeHolderA, vecA);
  myMatrix.check_coloring();

  // check dealing with multiple and unsorted input
  std::vector<std::tuple<size_t, size_t, bool>> vecC = {t_3, t_1, t_1, t_2, t_1,t_2, t_1};
  // check if method catches pathological input: index out of bound, transpose AND non-transposed
  //auto t_false_1 = std::make_tuple(3,1,false);
  //std::vector<std::tuple<size_t, size_t, bool>> vecC = {t_1, t_1, t_2, t_3, t_1,t_2, t_1, t_false_1};
  //auto t_false_2 = std::make_tuple(3,1,true);
  //std::vector<std::tuple<size_t, size_t, bool>> vecC = {t_1, t_1, t_2, t_3, t_1,t_2, t_1, t_false_1, t_false_2};
  myMatrix.replace_blocks(placeHolderB, vecC);
  myMatrix.print_coloring_pattern("After multiple input", true);
  myMatrix.print_and_check("After multiple input");
  }

  {
  //check if it works with tranposition
  ColoredBlockMatrix myBlockMatrix({3,7},{3,7});
  myBlockMatrix.check_coloring();
  auto t_2 = std::make_tuple(0,1,true);
  auto t_3 = std::make_tuple(1,0,false);
  std::vector<std::tuple<size_t, size_t, bool>> vec = {t_2,t_3};

  TMatrix placeHolderA(7,3);

  myBlockMatrix.replace_blocks(placeHolderA, vec);
  myBlockMatrix.print_and_check("insert of transp");

  }
  Output::print("Test program finished.");

}


