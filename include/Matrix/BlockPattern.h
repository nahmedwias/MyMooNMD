#ifndef __BLOCKPATTERN__
#define __BLOCKPATTERN__

#include <vector>

/** ************************************************************************ 
*
* @class      BlockPattern
* @brief      represent a pattern of (possibly local) matrices 
* 
*             This class describes the block structure of a BlockMatrix.
*             One main purpose is to identify which blocks are formed by which
*             Finite Element spaces.
*
*             
* 
* @todo extend this to handle NSTYPE 1,2,...
* 
* @author     Ulrich Wilbrandt
* @date       08.09.2015
*
* @ruleof0
*
****************************************************************************/

/** @brief the supported problem types when creating a BlockPattern */
enum class Problem_type 
{ ConvDiffReac, NavierStokes, Darcy, dummy };
std::ostream& operator<<(std::ostream& out, const Problem_type value);
std::string to_string(const Problem_type value);

class BlockPattern
{
 private:
  /** @brief the type of problem 
   * 
   * This is not strictly needed to use this class, i.e. it can be set to 
   * Problem_type::dummy.
   */
  Problem_type type;
  
  /** @brief the space dimension
   * 
   * This is not strictly needed to use this class. Then it is set to 0.
   */
  unsigned int space_dim;
   
  /** @brief the number of Finite Element spaces used, 
   * 
   * There is at least one and usually no more than two spaces.
   * 
   * In fact you can also create a BlockPattern without any Finite Element 
   * spaces, then n_spaces is typically the number of different block 
   * dimensions.
   */
  unsigned int n_spaces;
  
  /// @brief the number of block rows
  unsigned int n_block_rows;
  /// @brief the number of block columns
  unsigned int n_block_columns;
  
  /// @brief index of the row space for each block (smaller than n_spaces)
  std::vector<unsigned int> row_spaces;
  /// @brief index of the column space for each block (smaller than n_spaces)
  std::vector<unsigned int> column_spaces;
  
  /** @brief store which different sparsity structures exist
   * 
   * The length of this vector is equal to the number of different structures, 
   * which is usually the number of Finite Element spaces squared. In other 
   * words the number of possible combinations which occur as test/ansatz spaces
   * in a BlockMatrix using this BlockPattern.
   * 
   * If this->patterns[i] = j, then j >= i and the i-th pattern is the one 
   * of the j-th block.
   * 
   * That means this->patterns[i] is always the first block to have the 
   * corresponding sparsity pattern.
   */
  std::vector<unsigned int> patterns;
  
  /** @brief pattern of each block
   * 
   * The size of this vector is the number of blocks. If 
   * this->pattern_of_blocks[i] = j, then i >= j and the i-th block has the same
   * sparsity pattern as the this->patterns[j]-th block.
   */
  std::vector<unsigned int> pattern_of_blocks;
  
 public:
  /** @brief construct a BlockPattern for a given problem type 
   * 
   * This is probably the constructor you want to use for standard problems. 
   */
  BlockPattern(const Problem_type, unsigned int space_dimension, 
               bool mass_matrix = false);
  
  /** @brief construct a BlockPattern without Finite Element spaces
   * 
   * The vector fe_spaces contains only one nullptr. The vectors row_spaces and 
   * column_spaces are filled with zeros. The problem type is set to dummy.
   * The vector patterns only has one zero and the vector pattern_of_block is
   * filled with zeros.
   */
  BlockPattern(unsigned int n_rows, unsigned int n_cols);
  
  //Declaration of special member functions - rule of zero

  //! Default copy constructor. Performs deep copy.
  BlockPattern(const BlockPattern&) = default;

  //! Default move constructor.
  BlockPattern(BlockPattern&&) = default;

  //! Default copy assignment operator. Performs deep copy.
  BlockPattern& operator=(const BlockPattern&) = default;

  //! Default move assignment operator
  BlockPattern& operator=(BlockPattern&&) = default;

  //! Default destructor.
  ~BlockPattern() = default;

  /** @brief return true if the b-th block is on the diagonal, otherwise false
   * 
   * This is needed to put ones on the diagonal only of the correct matrix 
   * blocks.
   */
  bool diagonal_block(unsigned int b) const;
  
  /* getter functions */
  unsigned int get_n_spaces() const { return this->n_spaces; }
  unsigned int n_blocks() const { return n_block_rows * n_block_columns; }
  unsigned int n_rows() const {return n_block_rows; } // number of block rows
  unsigned int n_cols() const {return n_block_columns;}//number of block columns
  /** @brief return the number of \em different sparsity patterns */
  unsigned int n_patterns() const { return patterns.size(); }
  
  /** @brief print some information on this BlockPattern */
  void info(size_t verbose) const;
};

#endif // __BLOCKPATTERN__
