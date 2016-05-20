// =======================================================================
// @(#)Matrix.C        1.2 11/20/98
// 
// Class:       TMatrix
//
// Purpose:     store a  matrix (ansatz != test space)
//
// Author:      Gunar Matthies
//
// History:     26.08.1998 start implementation
//
// =======================================================================

#include <Matrix.h>
#include <string.h>
#include <Constants.h>
#include <LinAlg.h>
#include <stdlib.h>

#include <MooNMD_Io.h>
#include <fstream>

TMatrix::TMatrix(std::shared_ptr<TStructure> structure)
 : structure(structure), entries(this->structure->GetN_Entries(), 0.)
{
}

TMatrix::TMatrix(int nRows, int nCols)
 : TMatrix(std::make_shared<TStructure>(nRows,nCols))
{
  ;
}

void TMatrix::reset()
{
  memset(this->GetEntries(), 0., this->structure->GetN_Entries()*SizeOfDouble);
}

void TMatrix::setEntries(std::vector<double> entries)
{
  if(this->entries.size() != entries.size())
  {
    ErrThrow("resetting the entries of a matrix with the wrong number of ",
             "entries ", this->entries.size(), " != ", entries.size());
  }
  this->entries = entries;
}

void TMatrix::write(std::string filename) const
{
  std::ofstream matrixfile;
  matrixfile.open(filename.c_str());

  //write the header line - coordinate format, real values, no symmetry used
  matrixfile << "%%MatrixMarket matrix coordinate real general \n";

  //write general matrix information
  matrixfile << GetN_Rows() << "\t" << GetN_Columns() << "\t" << GetN_Entries() << "\n";

  //loop through matrix and print row - column - entry
  // for each entry in the sparsity structure
  const int* RowPtr = structure->GetRowPtr();
  const int* KCol = structure->GetKCol();
  int begin, end, pos=0;

  for (int i=0; i<structure->GetN_Rows(); ++i)
  {
    begin = RowPtr[i];
    end   = RowPtr[i+1];

    for (int j=begin; j<end; ++j)
    {
      // shift row and col by +1 (fortran style)
      matrixfile << i +1 << "\t" << KCol[pos] + 1 << "\t" << entries[pos] << "\n";
      ++pos;
    }
  }


  matrixfile.close();
}

void TMatrix::Print(const char *name) const
{
  const int* RowPtr = structure->GetRowPtr();
  const int* KCol = structure->GetKCol();
  int begin, end, pos=0;
  
  for (int i=0; i<structure->GetN_Rows(); ++i)
  {
    begin = RowPtr[i];
    end   = RowPtr[i+1];
    
    for (int j=begin; j<end; ++j)
    {
      Output::print<1>(name, "(", i, ",", KCol[pos], ") = ",
		       entries[pos],";");
      ++pos;
    }
  }
}

void TMatrix::PrintFull(std::string name, int fieldWidth) const
{
  const int* rowPtr = structure->GetRowPtr();
  const int* KCol = structure->GetKCol();

  cout << endl << name << " = " << endl;
  for (int curRow = 0; curRow < structure->GetN_Rows(); curRow++) 
  {
    int rowEnd = rowPtr[curRow+1];
    int posKCol = rowPtr[curRow];
    for (int curCol = 0; curCol < structure->GetN_Columns(); curCol++) 
    {
      if(curCol == KCol[posKCol] && posKCol < rowEnd)
      {
        cout << setw(fieldWidth) << entries[posKCol] << ", ";
        posKCol++;
      }
      else 
      {
        cout << setw(fieldWidth) << 0.0 << ", ";
      }
    }
    cout << endl;
  }
  cout << endl;
}

// add val to a matrix element
// return an error if the entry is not in the sparse structure
void TMatrix::add(int i,int j, double val)
{
  if(val != 0.0)
    this->get(i, j) += val;
}

void TMatrix::add(int i, std::map<int,double> vals, double factor)
{
  if(i < 0 || i > this->GetN_Rows())
  {
    ErrThrow("This matrix does not have a row ", i, 
             ".\nThe dimension of this matrix is ", this->GetN_Rows(), " x ",
             this->GetN_Columns());
  }
  const int* RowPtr = structure->GetRowPtr();
  const int* KCol = structure->GetKCol();
  std::map<int,double>::iterator it = vals.begin();
  for (int m=RowPtr[i];m < RowPtr[i+1] && it != vals.end(); m++) 
  {
    if (KCol[m] == it->first) 
    {
      entries[m] += factor*it->second;
      ++it;
    }
  }
  if(it != vals.end())
  {
    ErrThrow("Error in TMatrix::add. There are entries in 'vals' which are ",
             "not in the sparse structure. row ", i, ", column ", it->first,
             ".\nExit\n");
  }
}

void TMatrix::add(std::map<int, std::map<int,double> > vals, double factor)
{
  // add every row to the matrix
  std::map<int,std::map<int,double> >::iterator it;
  for(it = vals.begin(); it != vals.end(); ++it)
    add(it->first, it->second, factor);
}
  
  

// set val of a matrix element
// return an error if the entry is not in the sparse structure
void TMatrix::set(int i,int j, double val)
{
  int index_of_entry = this->structure->index_of_entry(i, j);
  if(index_of_entry == -1)
  {
    ErrThrow("Entry ", i,",",j," not in sparsity structure.");
  }
  this->entries[index_of_entry] = val;
}

// get val of a matrix element
// return an error if the entry is not in the sparse structure
const double& TMatrix::get(int i,int j) const
{
  // index of the entry (i,j) within the vector Entries
  int index = this->structure->index_of_entry(i, j);
  if(index >= 0 )
    return this->entries[index];
  ErrThrow("could not find the entry (", i, ",", j,
           ") in the sparsity structure");
}

double& TMatrix::get(int i,int j)
{
  // index of the entry (i,j) within the vector Entries
  int index = this->structure->index_of_entry(i, j);
  if(index >= 0 )
    return this->entries[index];
  ErrThrow("could not find the entry (", i, ",", j,
           ") in the sparsity structure");
}

double & TMatrix::operator()(const int i, const int j)
{
  return this->get(i, j);
}

const double & TMatrix::operator()(const int i, const int j) const
{
  return this->get(i, j);
}

double TMatrix::GetNorm(int p) const
{
  double result = 0.0;
  switch(p)
  {
    case -2:
      result = Dnorm(this->GetN_Entries(), this->GetEntries());
      break;
    case -1:
    {
      const int * rows = this->GetRowPtr();
      for(int row=0; row<this->GetN_Rows(); row++)
      {
        double row_sum = 0.0;
        //#pragma omp parallel for
        for(int i=rows[row]; i<rows[row+1]; i++)
        {
          row_sum += fabs(entries[i]);
        }
        if(row_sum>result && row_sum!=1.0)
          result = row_sum;
      }
      break;
    }
    case 0:
      for(int i=0; i<this->GetN_Entries(); i++)
      {
        double a = fabs(entries[i]);
        if(a > result)
          result = a;
      }
      break;
    case 1:
      ErrThrow("maximum absolute column sum norm of a matrix not yet ",
               "implemented!");
      break;
    case 2:
      ErrThrow("spectral norm of a matrix not yet implemented!");
      break;
    default:
      ErrThrow("undefined norm of a matrix!");
      break;
  }
  return result;
}


double* operator*(const TMatrix & A,const double* x)
{
  const double *AEntries = A.GetEntries();
  const int *ARowPtr = A.GetRowPtr();
  const int *AColIndex = A.GetKCol();

  int nrows = A.GetN_Rows();
  
  double *y = new double[nrows];
  double value;
  int index;

  for(int i=0;i<nrows;i++) 
  {
    value = 0;
//#pragma omp parallel for
    for (int j=ARowPtr[i]; j<ARowPtr[i+1]; j++) 
    {

      index = AColIndex[j];
      value += AEntries[j] * x[index];
    }
    y[i] = value;
  }
  return y;
}

TMatrix & TMatrix::operator+=(const TMatrix* A)
{
  if(this->GetStructure() != A->GetStructure()) // compare objects
  {
    ErrThrow("TMatrix::operator+= : the two matrices do not match.");
  }
  
  int n_entries = this->GetN_Entries();
  const double *AEntries = A->GetEntries();
  for(int i = 0; i < n_entries; ++i)
  {
    this->entries[i] += AEntries[i];
  }
  return *this;
}



TMatrix & TMatrix::operator-=(const TMatrix* A)
{
  if(this->GetStructure() != A->GetStructure()) // compare objects
  {
    ErrThrow("TMatrix::operator-= : the two matrices do not match.");
  }
  
  int n_entries = this->GetN_Entries();
  const double *AEntries = A->GetEntries();
  for(int i = 0; i < n_entries; ++i)
  {
    this->entries[i] -= AEntries[i];
  }
  return *this;
}


void TMatrix::multiply(const double * const x, double * const y, double a) const
{
  if(a == 0.0)
    return;
  
  const int *rowPtr = GetRowPtr();
  const int *colIndex = GetKCol();

  int nrows = GetN_Rows();
  
  double value;
  int i, j, end;
  for(i = 0; i < nrows; i++) 
  {
    value = 0;
    end = rowPtr[i+1];
    for (j = rowPtr[i]; j < end; j++) 
    {
      value += entries[j] * x[colIndex[j]];
    }
    y[i] += a * value;
  }
}

void TMatrix::transpose_multiply(const double * const x, double * const y, double a)
      const
{
  if(a == 0.0)
    return;
  
  const int *rowPtr = GetRowPtr();
  const int *colIndex = GetKCol();
  
  int nrows = GetN_Rows();
  
  double value = 0.;
  int i, j, end;
  for(i = 0; i < nrows; i++)
  {
    value = a*x[i];
    end = rowPtr[i + 1];
    for(j = rowPtr[i]; j < end; j++)
    {
      // Entries[j] is the (i,colIndex[j])-th entry of this matrix
      y[colIndex[j]] += entries[j] * value;
    }
  }
}

TMatrix* TMatrix::multiply(const TMatrix * const B, double a) const
{
  const int n_A_rows = this->GetN_Rows();   // = n_C_rows
  const int n_A_cols = this->GetN_Columns();
  const int n_B_rows = B->GetN_Rows();
  
  if(n_A_cols != n_B_rows)
  {
    ErrThrow("dimention mismatch during matrix-matrix multiplication");
  }
  const int * const a_rows = this->GetRowPtr();
  const int * const a_cols = this->GetKCol();
  
  const TStructure & strucB = B->GetStructure();
  const double * const b_entries = B->GetEntries();
  
  std::shared_ptr<TStructure> struc_c = get_product_structure(this->GetStructure(), strucB);
  const int * c_rows = struc_c->GetRowPtr();
  const int * c_cols = struc_c->GetKCol();
  TMatrix * c = new TMatrix(struc_c);
  double * c_entries = c->GetEntries();
  
  // fill the entries
  // loop over all rows in C
//#pragma omp parallel for
  for(int row = 0; row < n_A_rows; row++)
  {
    
    // loop over all entries in this row in C
    for(int col = c_rows[row]; col < c_rows[row + 1]; col++)
    {
      // multiply 'this row of A' x 'this column of B'
      // loop over all entries in this row in A
      for(int i = a_rows[row]; i < a_rows[row+1]; i++)
      {
        int ib = strucB.index_of_entry(a_cols[i], c_cols[col]);
        if(ib != -1)
        {
          c_entries[col] += entries[i] * b_entries[ib];
        }
      }
    }
  }
  return c;
}

TMatrix* TMatrix::multiply_with_transpose_from_right() const{

  // put up a unity diagonal matrix as a vector
  std::vector<double> unityScaling(structure->GetN_Columns(), 1.0);

  // let the other implementations do the work.
  return multiply_with_transpose_from_right(unityScaling);
}

TMatrix* TMatrix::multiply_with_transpose_from_right(
  const std::vector<double>& diagonalScaling) const
{
  // put up the structure by call to the specific structure generating method
  TStructure * productStructure = 
    structure->get_structure_of_product_with_transpose_from_right();

  // construct the matrix and return a pointer
  TMatrix * ret =  multiply_with_transpose_from_right(diagonalScaling, 
                                                      *productStructure);
  delete productStructure; // this has been (deep) copied in the above method
  return ret;
}

TMatrix* TMatrix::multiply_with_transpose_from_right(
  const std::vector<double>& diagonalScaling, const TStructure& knownStructure) 
const
{
  //check if the dimensions match
  if((int)diagonalScaling.size() != structure->GetN_Columns())
  {
    ErrThrow("Dimension mismatch! ", diagonalScaling.size(), "  ",
             structure->GetN_Columns());
  }

  const size_t nProductRows = structure->GetN_Rows();

  // copy construct the product's TStructure
  std::shared_ptr<TStructure> productStructure = 
    std::make_shared<TStructure>(knownStructure);
  // create new matrix with this structure
  TMatrix* product = new TMatrix(productStructure);
  
  int * productRowPtr = productStructure->GetRowPtr();
  double * productEntries = product->GetEntries();

  // fill the entries
  // loop over all rows in the product
  for(size_t row = 0; row < nProductRows; row++)
  {
    // loop over all entries in this row in the product
    for(int iEntries = productRowPtr[row]; iEntries < productRowPtr[row + 1];
        ++iEntries)
    {
      // hold row1 and row2 to fix ideas
      int row1 = row;
      int row2 = productStructure->GetKCol()[iEntries];

      //store begin and end indices for the columns and entries segments
      int beginRow1 = structure->GetRowPtr()[row1];
      int endRow1 = structure->GetRowPtr()[row1+1];
      int beginRow2 = structure->GetRowPtr()[row2];
      int endRow2 = structure->GetRowPtr()[row2+1];

      // work on two segments of column array
      const int* row1ColBegin = &structure->GetKCol()[beginRow1];
      const double* row1EntriesBegin = &entries[beginRow1];
      const int row1SegmentSize = endRow1 - beginRow1;
      const int* row2ColBegin = &structure->GetKCol()[beginRow2];
      const double* row2EntriesBegin = &entries[beginRow2];
      const int row2SegmentSize = endRow2 - beginRow2;

      //initialize the value to be written later on
      double temp = 0;
      //better use index for control of the loop!
      size_t index1 = 0;
      size_t index2 = 0;
      while ((int)index1 < row1SegmentSize && (int)index2 < row2SegmentSize)
      {
        if (row1ColBegin[index1] > row2ColBegin[index2])
        {
          index2++;
        }
        else if (row2ColBegin[index2] > row1ColBegin[index1])
        {
          index1++;
        }
        else
        {
          //we found a pair of indices with equal entry in KCol: account for 
          // the scaling matrix!
          temp +=  row1EntriesBegin[index1]
                 * row2EntriesBegin[index2]
                 * diagonalScaling[row1ColBegin[index1]];
          index1++;
          index2++;
        }
      }
      //set the current entry to be the freshly calculated vector product
      productEntries[iEntries]= temp;
    } //end loop over all entries in this row
  } //end loop over all rows of the product

  // return the pointer
  return product;
}

std::shared_ptr< TMatrix > 
  TMatrix::multiply_with_transpose_from_right(const TMatrix& B) const
{
  // dimension check
  if(B.GetN_Rows() != this->GetN_Columns() 
    || B.GetN_Columns() !=this->GetN_Columns())
  {
    ErrThrow("Dimension mismatch  ", B.GetN_Rows(), "  ", this->GetN_Columns());
  }
  
  // construct a product structure
  std::shared_ptr<TStructure> productStructure(
    structure->get_structure_of_product_with_transpose_from_right(B.GetStructure()));
  
  // lambda function which returns the entry in the BA^T matrix  
  auto return_BAT_entry = [this, B](int i, int j)
  {
    const int* row1ColBegin = &B.GetKCol()[B.GetRowPtr()[i]];
    const int row1ColSize = B.GetRowPtr()[i+1] - B.GetRowPtr()[i];
    const double *entriesB = B.GetEntries();
    
    const int* row2ColBegin = &this->GetKCol()[this->GetRowPtr()[j]];
    const int row2ColSize = this->GetRowPtr()[j+1] 
                               - this->GetRowPtr()[j];
    
    size_t indexB = 0;
    size_t indexA = 0;
    
    double entry_in_BAT_product = 0;
    while ((int)indexB < row1ColSize && (int)indexA < row2ColSize)
    {
      if (row1ColBegin[indexB] > row2ColBegin[indexA])
      {
        indexA++;
      }
      else if (row2ColBegin[indexA] > row1ColBegin[indexB])
      { 
        indexB++;
      }
      else
      {        
        entry_in_BAT_product +=  entriesB[indexB+B.GetRowPtr()[i]]
                               * entries[indexA+this->GetRowPtr()[j]];
        indexB++;
        indexA++;
      }
    }
    return entry_in_BAT_product;
  };
  
  // number of rows in product structure
  const size_t nProductRows = productStructure->GetN_Rows();
  int * productRowPtr = productStructure->GetRowPtr();
  
  // create product matrix 
  std::shared_ptr<TMatrix> productMatrix 
          = std::make_shared<TMatrix>(productStructure);
  // entries in the product matrix
  double * productEntries = productMatrix->GetEntries();
  
  // fill the entries
  // loop over all rows in the product
  for(unsigned int row=0; row<nProductRows; row++)
  {
    unsigned int begin = productRowPtr[row];
    unsigned int end = productRowPtr[row+1];
    // loop over entries in "this row" in the product 
    for(unsigned int iEntries = begin; iEntries<end; iEntries++)
    {
      int rowABAT = row;
      int colABAT = productStructure->GetKCol()[iEntries];
      
      // store begin and end indices for columns and entries segments
      int beginRowA = this->GetStructure().GetRowPtr()[rowABAT];
      int endRowA   = this->GetStructure().GetRowPtr()[rowABAT+1];
      
      double entry_ABAT = 0;
      for(int k=beginRowA; k<endRowA; ++k)
      {
        // compute the entry in the BA^T
        double BAT_entry = return_BAT_entry(this->GetKCol()[k], colABAT);
        double a_entry = entries[k];
        entry_ABAT += a_entry * BAT_entry;
      }
      // set the current entry
      productEntries[iEntries] = entry_ABAT;
    } //endfor loop over all entries in this row
  }// endfor loop over all entries in the product 
  // return point of this matrix
  return productMatrix;
}

TMatrix* TMatrix::GetTransposed() const
{
  // get transposed structure
  std::shared_ptr<TStructure> structureT = structure->GetTransposed();
  const int * rowsT= structureT->GetRowPtr();
  const int * colsT= structureT->GetKCol();
  const int * rows = this->GetRowPtr();
  const int * cols = this->GetKCol();
  
  // transpose the entries:
  TMatrix * mT = new TMatrix(structureT);
  double *entriesT = mT->GetEntries();
  // loop over all rows of the original matrix
  for(int i=0; i<this->GetN_Rows(); i++)
  {
    // loop over all entries in this row of the original matrix
    for(int j=rows[i]; j<rows[i+1]; j++)
    {
      // cols[j] is the column of this entry in the original matrix,
      // it corresponds to a row of the transposed matrix
      // look for the column index in that row in the transposed matrix 
      // which equals this (non-transposed) row index
      for(int k=rowsT[cols[j]]; k<rowsT[cols[j]+1] ;k++)
      {
        if(i==colsT[k])
        {
          entriesT[k] = entries[j];
          continue; // entry found
        }
      }
    }
  }
  
  return mT;
}

void TMatrix::remove_zeros(double tol)
{
  if(!this->structure.unique()) //could be adapted in time: work on a copy
    ErrThrow("Cannot remove zeroes in a structure shared by multiple matrices!");

  if(this->structure->GetHangingN_Entries() != 0) //if you need hanging nodes, fix this
    ErrThrow("Matrix structure has hanging node entries. "
        "TMatrix::remove_zeros will not yet reset their number and arrays correctly.");

  if(tol < 0)
    tol = this->GetNorm(0) * 1e-15; // largest (in magnitude) entry

  int *row_ptr = structure->GetRowPtr();
  int n_rows = structure->GetN_Rows();
  int* kcol = structure->GetKCol();
  //"entries" will be treated as the std::vector it is

  int row_begin_old = row_ptr[0];

  for(int i = 0 ; i<n_rows ; ++i)
  {
    int row_end_old = row_ptr[i+1];
    std::vector<int> kcol_part_new;
    std::vector<double> entries_part_new;
    kcol_part_new.reserve( row_end_old -row_begin_old );
    entries_part_new.reserve( row_end_old -row_begin_old );

    for(int j = row_begin_old; j < row_end_old; ++j)
    {
      if(!(fabs(entries.at(j) - 0) < tol)) //entry does not count as zero
      {
        kcol_part_new.push_back(kcol[j]);
        entries_part_new.push_back(entries.at(j));
      }
    } //new kcol and entries part for row i are filled

    //copy new parts into old arrays
    int row_begin_new = row_ptr[i]; //was updated in loop step i-1
    size_t n_entries_row_new = kcol_part_new.size();
    if(n_entries_row_new != 0) //this is not a zero row
    {
      memcpy(&kcol[row_begin_new], &kcol_part_new.at(0), n_entries_row_new*SizeOfInt);
      memcpy(&entries.at(row_begin_new), &entries_part_new.at(0), n_entries_row_new*SizeOfDouble); //stl here!!
    }
    //update row_ptr
    row_begin_old = row_ptr[i+1];
    row_ptr[i+1] = row_begin_new + n_entries_row_new;
  }

  //Re-fit the TStructure (nEntries and columns array)
  structure->reset_n_entries();
  //Re-fit the entries array (throw out all trailing (nearly) zeroes)
  entries.resize(structure->GetN_Entries());

}


void TMatrix::add_scaled(const TMatrix& m, double factor)
{
  if(this->GetStructure() != m.GetStructure()) // compare objects
  {
    ErrThrow("TMatrix::add : the two matrices do not match.");
  }
  Daxpy(this->GetN_Entries(), factor, m.GetEntries(), this->GetEntries());
}

void TMatrix::scale(double factor)
{
  Dscal(this->GetN_Entries(), factor, this->GetEntries());
}

void TMatrix::scale(const double * const factor, bool from_left)
{
  const int *rowPtr = GetRowPtr();
  const int *colIndex = GetKCol();
  
  if(from_left)
  {
    for(int i = 0, nrows = GetN_Rows(); i < nrows; i++) 
    {
      int end = rowPtr[i+1];
      // scale entire row with the same factor
      for(int j = rowPtr[i]; j < end; j++) 
      {
        entries[j] *= factor[i];
      }
    }
  }
  else
  {
    for(int i = 0, nrows = GetN_Rows(); i < nrows; i++) 
    {
      int end = rowPtr[i+1];
      for(int j = rowPtr[i]; j < end; j++) 
      {
        // scale columnwise
        entries[j] *= factor[colIndex[j]];
      }
    }
  }
}

TMatrix & TMatrix::operator*=(const double a)
{
  this->scale(a);
  return *this;
}


void TMatrix::reorderMatrix() 
{
  // make a deep copy of the structure in case of other matrices sharing it
  this->copyOwnStructure();
  int l, begin, end;
  double value;
  
  const int* Row = structure->GetRowPtr();
  int* KCol = structure->GetKCol();
  
  //OutPut("\tPARDISO: reordering of the columns will be performed"<<endl);
  //OutPut("\tPARDISO:   no back ordering implemented !!!"<<endl);
  
  for(int i=0;i<structure->GetN_Rows();i++) 
  {
    begin=Row[i];
    end=Row[i+1];
    for(int j=begin;j<end;j++) 
    {
      for(int k=j+1;k<end;k++) 
      {
        if(KCol[j] > KCol[k]) 
        {
          l = KCol[j];       value = entries[j];
          KCol[j] = KCol[k]; entries[j] = entries[k];
          KCol[k] = l;       entries[k] = value;
        }
      }
    }
  }
}

/** ************************************************************************* */
void TMatrix::copyOwnStructure()
{
  if(this->structure.unique())
  {
    // no one else knows this->structure
    return;
  }
  // create a deep copy of the own structure
  std::shared_ptr<TStructure> new_structure(new TStructure(*this->structure));
  this->structure = new_structure;
}


/** ************************************************************************* */
void TMatrix::info(size_t verbose) const
{
  OutPut(" TMatrix M with " << this->GetN_Rows() << " rows and "
         << this->GetN_Columns() << " columns\n");
  if(verbose > 2)
  {
    this->Print("M");
  }
}
