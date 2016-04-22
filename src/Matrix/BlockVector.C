#include <BlockVector.h>
#include <BlockMatrix.h>
#include <BlockFEMatrix.h>
#include <stdlib.h>
#include <Constants.h>
#include <LinAlg.h>

/** ************************************************************************ */
BlockVector::BlockVector() : entries(), lengths(), actives()
{
  Output::print<3>("Constructor of BlockVector with no arguments");
}

/** ************************************************************************ */
BlockVector::BlockVector(unsigned int length)
 : entries(length, 0.0), lengths(1, length), actives(1, length)
{
  Output::print<3>("Constructor of BlockVector with length ", length);
}

/** ************************************************************************ */
BlockVector::BlockVector(int length)
 : BlockVector(length > 0 ? (unsigned int)length : 0)
{
  if(length < 0)
    ErrThrow("cannot construct BlockVector with negative number of entries");
}

BlockVector::BlockVector(const BlockMatrix& mat, bool result)
{
  Output::print<3>("Constructor of BlockVector using a BlockMatrix");
  // the total length of this vector
  size_t total_length = result ? mat.get_n_total_rows() : mat.get_n_total_columns();
  // number of blocks in this BlockVector
  size_t n_blocks = result ? mat.get_n_cell_rows() : mat.get_n_cell_columns();

  entries.resize(total_length, 0.0);
  // set all entries to zero
  std::fill(entries.begin(), entries.end(), 0.0);
  lengths.resize(n_blocks, 0);
  actives.resize(n_blocks, 0);

  //determine the length of each vector block
  for(size_t b = 0; b < n_blocks; b++)
  {
    if (result)
    {//vector on rhs
      lengths[b] = mat.get_n_rows_in_cell(b,0);
    }
    else
    {//vector as factor
      lengths[b] = mat.get_n_columns_in_cell(0,b);
    }
   // all entries are active
    actives[b] = lengths[b];
  }
}

BlockVector::BlockVector(const BlockFEMatrix& mat, bool result)
: BlockVector(static_cast<const BlockMatrix&>(mat), result)
{
  //now set actives correctly
  for(size_t b = 0; b < n_blocks(); ++b)
  {//go thorugh all blocks of the vector
    if( result )
    { // take actives from row space
      actives[b] = mat.get_n_row_actives( b );
    }
    else
    { //take actives from column space
      actives[b] = mat.get_n_column_actives( b );
    }
  }
}

/** ************************************************************************ */
unsigned int BlockVector::offset(unsigned int b) const
{
  if(b >= this->n_blocks())
    ErrThrow("trying to access block ", b, ", but there are only ", 
             this->n_blocks(), " blocks in this BlockVector");
  return std::accumulate(lengths.cbegin(), lengths.cbegin() + b, 0);
}

/** ************************************************************************ */
void BlockVector::reset()
{
  std::fill(this->entries.begin(), this->entries.end(), 0.0);
}

/** ************************************************************************ */
void BlockVector::ResetActive()
{
  unsigned int n_blocks = this->n_blocks();
  // loop over all blocks
  auto it = this->entries.begin();
  for(unsigned int i = 0; i < n_blocks; i++)
  {
    std::fill(it, it + this->active(i), 0.0);
    std::advance(it, this->length(i)); // it += this->length(i);
  }
}

/** ************************************************************************ */
void BlockVector::ResetNonActive()
{
  unsigned int n_blocks = this->n_blocks();
  // loop over all blocks
  auto it = this->entries.begin();
  for(unsigned int i = 0; i < n_blocks; i++)
  {
    std::fill(it + this->active(i), it + this->length(i), 0.0);
    std::advance(it, this->length(i)); // it += this->length(i);
  }
}

/** ************************************************************************ */
void BlockVector::scale(const double a, const unsigned int i)
{
  if((unsigned int)i < lengths.size())
    Dscal(lengths.at(i), a, this->block(i)); // scale i-th subvector
  else
    ErrThrow("trying to scale subvector ", i, " which does not exist");
}

/** ************************************************************************ */
void BlockVector::scaleActive(const double a)
{
  if(a==0)
    this->reset();
  else if(a != 1.0)
    for(unsigned int i=0; i<this->n_blocks(); ++i)  
    {
      Dscal(actives.at(i), a, this->block(i));
    }
}

/** ************************************************************************ */
void BlockVector::scale(const double a)
{
  if(a == 0.0)
    this->reset();
  else if(a != 1.0)
    Dscal(this->length(), a, this->get_entries());
  // else if a == 1.0 not action is taken
}


/** ************************************************************************ */
void BlockVector::add_scaled(const BlockVector& r, double factor)
{
  const unsigned int l = this->length(); // length of this BlockVector
  if(r.length() != l)
  {
    ErrThrow("unable to add two BlockVectors of different lengths\t", 
             l, "\t", r.length());
  }
  else if(this->n_blocks() != r.n_blocks())
  {
    Output::print("WARNING: BlockVector::operator+=\n adding to BlockVectors "
                  "with the same length but different numbers of blocks");
  }
  Daxpy(l, factor, r.get_entries(), this->get_entries());
}

/** ************************************************************************ */
void BlockVector::addScaledActive(const BlockVector& r, double factor)
{
  if(this->n_blocks() != r.n_blocks())
    ErrThrow("number of blocks in two vectors must be same");
  
  for(unsigned int i=0; i<this->n_blocks(); ++i)
  {
    Daxpy(this->actives.at(i), factor, r.block(i), this->block(i));
  }
}

/** ************************************************************************ */
void BlockVector::addScaledNonActive(const BlockVector& r, double factor)
{
  if(this->n_blocks() != r.n_blocks())
    ErrThrow("number of blocks in two vectors must be same");

  for(unsigned int i=0; i<this->n_blocks(); ++i)
  {
    if(this->actives.at(i) != r.active(i))
    {
      ErrThrow("Number of actives in block ",i," do not match!");
    }
    // do the daxpy for the nonactives only (includes pointer arithmetic...)
    int n_non_actives = this->lengths.at(i)- this->actives.at(i);
    Daxpy(n_non_actives , factor, r.block(i) + this->actives.at(i), this->block(i) + this->actives.at(i));
  }
}

/** ************************************************************************ */
void BlockVector::copy(const double * x, const int i)
{
  if(i < 0)
    *this = x; // copy entire vector
  else if((unsigned int)i < this->n_blocks())
    std::copy(x, x + this->length(i), this->block(i));
    //memcpy(this->block(i), x, this->length(i)*sizeof(double)); // in string.h
  else
  {
    ErrThrow("trying to copy to subvector ", i, " which does not exist.");
  }
}

/** ************************************************************************ */
void BlockVector::copy_nonactive(const BlockVector& r)
{
  const unsigned int l = this->length(); // length of this BlockVector
  if(r.length() != l)
  {
    ErrThrow("unable to copy from one BlockVector to another one of different ",
             "length\t", l, "\t", r.length());
  }
  else if(this->n_blocks() != r.n_blocks())
  {
    Output::print("WARNING: BlockVector::operator+=\n adding to BlockVectors "
                  "with the same length but different numbers of blocks\n");
  }
  
  for(unsigned int b = 0, n_b = this->n_blocks(); b < n_b; ++b)
  {
    if(this->lengths[b] > this->active(b))
      std::copy(r.block(b) + r.active(b), r.block(b) + r.length(b),
                this->block(b) + this->active(b));
      // in string.h
      //memcpy(this->block(b) + this->active(b), r.block(b) + this->active(b),
      //       (this->lengths[b] - this->active(b))*sizeof(double));
  }
}

/** ************************************************************************ */
void BlockVector::add(const double* x, const int i, double a)
{
  if(i < 0) 
  {
    const unsigned int l = length(); // length of this BlockVector
    Daxpy(l, a, x, this->get_entries());
  }
  else if((unsigned int)i < this->lengths.size()) 
  {
    double * b = this->block(i);
    for(unsigned int j = 0; j < this->length(i); j++)
    {
      b[j] += a*x[j];
    }
  }
  else 
  {
    ErrThrow("trying to add to subvector ", i, " which does not exist.");
  }
}


/** ************************************************************************ */
double BlockVector::norm() const
{
  return Dnorm(this->length(), this->get_entries());
}

/** ************************************************************************ */
void BlockVector::print(const std::string name, const int iB) const
{
  if(iB < 0)
  { // print full BlockVector
    for (unsigned int i = 0, l = this->length(); i < l; i++) 
      Output::print(name , "(" , i+1 , ")= " , this->at(i) , ";");
  }
  else if((unsigned int)iB < lengths.size())
  {
    const double *myBlock = this->block(iB);
    for (unsigned int i = 0; i < this->lengths[iB]; i++) 
      Output::print(name ,"(" , i+1 , ")= " , myBlock[i] ,";");
  }
  else
    ErrThrow("trying to print subvector ", iB, " which does not exist");
}

/** ************************************************************************ */
void BlockVector::write(std::string filename) const
{
  std::ofstream vectorfile;
  vectorfile.open(filename.c_str());

  //write the header line - array format, real values
  vectorfile << "%%MatrixMarket matrix array real general \n";

  //write general matrix information
  vectorfile << length() << "\t" << 1 << "\t" << "\n";

  //loop through matrix and print each entry
  for(auto it : entries)
  {
    vectorfile << it << "\n";
  }


  vectorfile.close();
}

/** ************************************************************************ */
void BlockVector::info()
{
  using namespace Output;
  Output::print(" | info to BlockVector :");
  unsigned int nb = n_blocks();
  Output::print(" | -- total length ", length(), "\tnumber of blocks ", nb);
  if(nb > 1)
  {
    for(unsigned int i = 0; i < nb; i++)
    {
      if(this->lengths[i] != this->actives[i])
      {
        Output::print(" | ---- block ", i, " has length ", this->lengths[i],
                      " (", this->actives[i], " active)");
      }
      else
      {
        Output::print(" | ---- block ", i, " has length ", this->lengths[i]);
      }
    }
  }
}

/** ************************************************************************ */
BlockVector& BlockVector::operator=(const double *r)
{
  if(this->get_entries() == r) 
    // both are the same, no copying necessary
    return *this;
  std::copy(r, r + this->length(), this->entries.begin());
  //memcpy(this->get_entries(), r, this->length()*sizeof(double));// in string.h
  return *this;
}

/** ************************************************************************ */
BlockVector& BlockVector::operator=(const double a)
{
  std::fill(this->entries.begin(), this->entries.end(), a);
  return *this;
}

/** ************************************************************************ */
BlockVector& BlockVector::operator*=(const double a)
{
  this->scale(a);
  return *this;
}

/** ************************************************************************ */
BlockVector& BlockVector::operator+=(const BlockVector& r)
{
  this->add_scaled(r, 1.0);
  return *this;
}

/** ************************************************************************ */
BlockVector& BlockVector::operator-=(const BlockVector& r)
{
  this->add_scaled(r, -1.0);
  return *this;
}

/** ************************************************************************ */
double dot(const BlockVector& a, const BlockVector& b)
{
  unsigned int l = a.length();
  if(b.length() != l)
  {
    ErrThrow("unable to compute dot product of two BlockVectors of different ", 
             "lengths\t", l, "\t", b.length());
  }
  else if(a.n_blocks() != b.n_blocks())
  {
    Output::print("WARNING: dot(BlockVector a, BlockVector b)\n computing the "
                  "dot product of two vectors with the same length but "
                  "different numbers of blocks");
  }
  return Ddot(l, a.get_entries(), b.get_entries());
}

/** ************************************************************************ */
void BlockVector::copy_structure(const BlockVector& r)
{
  entries.resize(r.length(), 0.0);
  lengths.resize(r.n_blocks());
  actives.resize(r.n_blocks());
  
  for(unsigned int i = 0; i < r.n_blocks(); i++)
  {
    lengths[i] = r.length(i);
    actives[i] = r.active(i);
  }
}

/** ************************************************************************ */
void BlockVector::write_fo_file(std::string filename)
{
  std::ofstream dat(filename);
  if(!dat)
  {
    ErrMsg("cannot open file '" << filename << "' to save data. return");
    return;
  }
  dat << this->length() << endl;
  dat.write((char *) this->get_entries(), sizeof(double) * this->length());
  dat.close();
}

/** ************************************************************************ */
void BlockVector::read_from_file(std::string filename)
{
  std::ifstream dat(filename);
  if(!dat)
  {
    ErrThrow("cannot open file '", filename, "' to read data");
  }
  
  // check if the size of this vector and the number of doubles in the file
  // coincide:
  unsigned int length_read;
  std::string line;
  std::getline(dat, line);
  std::istringstream parser( std::string( line.begin(), line.end() ) );
  parser >> length_read >> std::ws;
  if(!parser || parser.get() != EOF)
  {
    ErrThrow("formatting error, the first line in the file should contain ",
             "only a number idicating the number of entries to be read");
  }
  else if(this->length() == 0)
  {
    // the vector is not yet filled with any data, allocate the arrays:
    this->entries.resize(length_read, 0.0);
    lengths.resize(1, length_read);
    actives.resize(1, length_read);
  }
  else if(this->length() != length_read)
  {
    ErrThrow("unexpected size to be read. Expected: ", this->length(),
             "\tin file: ", length_read);
  }
  // now we are sure that (this->length == length_read) 
  
  dat.read((char *)this->get_entries(), sizeof(double) * length_read);
  dat.close();
}

/** ************************************************************************ */
double* BlockVector::block(const unsigned int i)
{
  return this->get_entries() + this->offset(i);
}

/** ************************************************************************ */
const double* BlockVector::block(const unsigned int i) const
{
  return this->get_entries() + this->offset(i);
}

/** ************************************************************************ */
double& BlockVector::at(const unsigned int i)
{
  try
  {
    return entries.at(i);
  }
  catch(...)
    ErrThrow("index out of bounds");
}

/** ************************************************************************ */
const double& BlockVector::at(const unsigned int i) const
{
  try
  {
    return entries.at(i);
  }
  catch(...)
  {
    ErrThrow("index out of bounds");
  }
}
