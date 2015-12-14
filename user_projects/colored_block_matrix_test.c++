/*
 * Unit testing of the colored block matrix.
 *
 *  Created on: Dec 10, 2015
 *      Author: bartsch
 */

#include <ColoredBlockMatrix.h>
#include <vector>
#include <MooNMD_Io.h>

int main(int argc, char* argv[])
{
  // construct a 3 times 4 colored block matrix filled with zero blocks
  // of the given dimension

  std::vector<size_t> cell_row_numbers = {5 , 5 , 2};
  std::vector<size_t> cell_column_numbers = {5, 5 , 2, 3};

  ColoredBlockMatrix myMatrix(cell_row_numbers, cell_column_numbers);

  myMatrix.check_coloring();

  Output::print("Test program finished.");
}


