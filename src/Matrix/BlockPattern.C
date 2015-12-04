#include <algorithm>
#include <MooNMD_Io.h>
#include <BlockPattern.h>


/* ************************************************************************* */
BlockPattern::BlockPattern(const Problem_type t, unsigned int space_dimension,
                           bool mass_matrix)
 : type(t), space_dim(space_dimension), n_spaces(0), n_block_rows(0), 
   n_block_columns(0), row_spaces(0), column_spaces(0), 
   patterns(0), pattern_of_blocks(0)
{
  // some checks, to make sure the input makes sense
  if(space_dim != 2 && space_dim != 3)
   ErrThrow("unsupported space dimension ", space_dim);
  // check if all space were defined on the same Mesh
  
  switch(this->type)
  {
    case Problem_type::ConvDiffReac:
      this->n_spaces = 1;
      this->n_block_rows = 1;
      this->n_block_columns = 1;
      this->row_spaces = {0};
      this->column_spaces = {0};
      this->patterns = {0}; // only one pattern
      this->pattern_of_blocks = {0};
      break;
    case Problem_type::NavierStokes:
    {
      this->n_spaces = 2;
      this->n_block_rows = space_dim + 1;
      this->n_block_columns = space_dim + 1;
      if(space_dim == 2)
      {
        // ( a11 a12 a13 )   ( A11 A12 BT1 )
        // ( a21 a22 a23 ) = ( A21 A22 BT2 )
        // ( a31 a32 a33 )   ( B1  B2  C   )
        // a11, a12, a21, a22, a33, a13, a23, a31, a32
        // 0 is velocity space, 1 is pressure
        this->row_spaces =    {0, 0, 0, 0, 0, 0, 1, 1, 1};
        this->column_spaces = {0, 0, 1, 0, 0, 1, 0, 0, 1};
        // four different patterns:
        // velocity-velocity, velocity-pressure, 
        // pressure-velocity, pressure-pressure,
        this->patterns = {0, 2, 6, 8}; 
        // the first four blocks are the velocity-velocity blocks
        // the fourth is the pressure-pressure block
        // the sixth and seventh are the velocity-pressure blocks
        // the last two are the pressure-velocity blocks
        this->pattern_of_blocks = {0, 0, 1, 0, 0, 1, 2, 2, 3};
      }
      else if(space_dim == 3)
      {
        // ( a11 a12 a13 a14 )   ( A11 A12 A13 BT1 )
        // ( a21 a22 a23 a24 ) = ( A21 A22 A23 BT2 )
        // ( a31 a32 a33 a34 )   ( A31 A32 A33 BT3 )
        // ( a41 a42 a43 a44 )   ( B1  B2  B3  C   )
        // 0 is velocity space, 1 is pressure
        this->row_spaces =    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1};
        this->column_spaces = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
        // four different patterns:
        // velocity-velocity, velocity-pressure, 
        // pressure-velocity, pressure-pressure,
        this->patterns = {0, 3, 12, 15};
        // the first nine blocks are the velocity-velocity blocks
        // the 9th is the pressure-pressure block
        // the 10th, 11th and 12th are the velocity-pressure blocks
        // the last three are the pressure-velocity blocks
        this->pattern_of_blocks = {0,0,0,1 ,0,0,0,1, 0,0,0,1, 2,2,2,3};
      }
      /*if(mass_matrix) // time dependent case
      {
        ErrThrow("time dependent case not yet implemented for Navier--Stokes");
      }*/
      break;
    }
    case Problem_type::Darcy:
      this->n_spaces = 2;
      this->n_block_rows = 2;
      this->n_block_columns = 2;
      // 0 is velocity space, 1 is pressure
      // a11, a22, a12, a21
      this->row_spaces = {0, 1, 0, 1};
      this->column_spaces = {0, 1, 1, 0};
      this->patterns = {0, 1, 2, 3}; // four different patterns
      this->pattern_of_blocks = {0, 1, 2, 3}; // each block is different
      break;
    default:
      ErrThrow("unknown problem type ", t);
      break;
  }
}

/* ************************************************************************* */
BlockPattern::BlockPattern(unsigned int n_rows, unsigned int n_cols)
 : type(Problem_type::dummy), space_dim(0),
   n_spaces(std::max(n_rows, n_cols)), n_block_rows(n_rows),
   n_block_columns(n_cols), row_spaces(n_rows*n_cols, 0), 
   column_spaces(n_rows*n_cols, 0), patterns({0}), 
   pattern_of_blocks(n_rows*n_cols, 0)
{
  
}

/* ************************************************************************* */
unsigned int BlockPattern::test_space_index_of_block(unsigned int block) const
{
  if(block >= this->n_blocks())
  {
    ErrThrow("This BlockPattern only has ", this->n_blocks(), ", not ", block);
  }
  return this->row_spaces[block];
}
  
/* ************************************************************************* */
unsigned int BlockPattern::ansatz_space_index_of_block(unsigned int block) const
{
  if(block >= this->n_blocks())
  {
    ErrThrow("This BlockPattern only has ", this->n_blocks(), ", not ", block);
  }
  return this->column_spaces[block];
}

/* ************************************************************************* */
unsigned int BlockPattern::test_space_index_of_row(unsigned int r) const
{
  if(r >= this->n_rows())
  {
    ErrThrow("There are only ", this->n_rows(), " block rows, cannot get ",
             "the index of the test space of the ", r, "-th row");
  }
  return this->row_spaces[r*this->n_cols()];
}
  
/* ************************************************************************* */
unsigned int BlockPattern::ansatz_space_index_of_column(unsigned int c) const
{
  if(c >= this->n_cols())
  {
    ErrThrow("There are only ", this->n_cols(), " block columns, cannot get ", 
             "the index of the ansatz space of the ", c, "-th column");
  }
  return this->column_spaces[c];
}

/* ************************************************************************* */
unsigned int BlockPattern::get_space_dim() const
{
  return this->space_dim;
}

/* ************************************************************************* */
unsigned int BlockPattern::get_test_space_of_pattern(unsigned int p) const
{
  if(p >= this->n_patterns())
    ErrThrow("There are only ", this->n_patterns(), 
             " different sparsity patterns in this BlockPattern, ", p);
  return this->row_spaces[this->patterns[p]];
}

/* ************************************************************************* */
unsigned int BlockPattern::get_ansatz_space_of_pattern(unsigned int p) const
{
  if(p >= this->n_patterns())
    ErrThrow("There are only ", this->n_patterns(), 
             " different sparsity patterns in this BlockPattern, ", p);
  return this->column_spaces[this->patterns[p]];
}

/* ************************************************************************* */
unsigned int BlockPattern::pattern_of_block(unsigned int b) const
{
  if(b >= this->n_blocks())
  {
    ErrThrow("This BlockPattern only has ", this->n_blocks(), " blocks, not ",
             b);
  }
  return this->pattern_of_blocks[b];
}

/* ************************************************************************* */
bool BlockPattern::diagonal_block(unsigned int b) const
{
  if(b == 0) return true;
  if(b >= this->n_blocks())
  {
    ErrThrow("unable to tell if given block is on diagonal, ",
             "index out of range");
  }
  // the i-th diagonal entry is the i(m+1)-block where m is the number of 
  // columns.
  unsigned int diagonal = b/(this->n_block_columns+1); // floor function
  if(diagonal * (this->n_block_columns+1) == b)
    return true;
  else
    return false;
}

/* ************************************************************************* */
void BlockPattern::info(size_t verbose) const
{
  if(this->type != Problem_type::dummy)
  {
    OutPut("BlockPattern for a " << this->type << " problem\n");
  }
  else
  {
    OutPut("BlockPattern without a particular problem type\n");
  }
  
  OutPut(" number of blocks: " << this->n_blocks() << endl);
  OutPut(" number of block rows and columns: " << this->n_rows() << " and " <<
         this->n_cols() << endl);
  
  if(verbose > 0 && this->type != Problem_type::dummy)
  {
    OutPut(" number of FE_spaces " << this->n_spaces << endl);
    OutPut(" number of patterns " << this->patterns.size() << endl);
  }
  if(verbose > 1 && this->type != Problem_type::dummy)
  {
    for(unsigned int b = 0, n_b = this->n_blocks(); b < n_b; ++b)
    {
      OutPut(" block " << b << " has row space " << this->row_spaces.at(b) << 
             " and column space " << this->column_spaces.at(b) << 
             " and uses the pattern number " << this->pattern_of_blocks.at(b)
             << endl);
    }
  }
}

/* ************************************************************************* */
std::ostream& operator<<(std::ostream& out, const Problem_type value)
{
  const char* s = 0;
#define PROCESS_VAL(p) case(Problem_type::p): s = #p; break;
  switch(value)
  {
    PROCESS_VAL(ConvDiffReac);
    PROCESS_VAL(NavierStokes);
    PROCESS_VAL(Darcy);
    PROCESS_VAL(dummy);
    default: s = "unknown Problem_type "; break;
  }
#undef PROCESS_VAL
  return out << s;
}

/* ************************************************************************* */
std::string to_string(const Problem_type value)
{
#define PROCESS_VAL(p) case(Problem_type::p): return std::string(#p); break;
  switch(value)
  {
    PROCESS_VAL(ConvDiffReac);
    PROCESS_VAL(NavierStokes);
    PROCESS_VAL(Darcy);
    PROCESS_VAL(dummy);
    default: return std::string("unknown Problem_type "); break;
  }
#undef PROCESS_VAL
}

