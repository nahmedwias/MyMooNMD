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

TMatrix::TMatrix(std::shared_ptr<TStructure> structure)
 : structure(structure), entries(this->structure->GetN_Entries(), 0.)
{
}

void TMatrix::reset()
{
  memset(this->GetEntries(), 0., this->structure->GetN_Entries()*SizeOfDouble);
}

void TMatrix::setEntries(double* entries)
{
  std::copy(entries, entries + this->GetN_Entries(), this->entries.begin());
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
  this->entries[this->structure->index_of_entry(i, j)] = val;
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


void TMatrix::multiply(const double * const x, double *y, double a) const
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

void TMatrix::transpose_multiply(const double * const x, double *y, double a)
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
  if(tol < 0)
    tol = this->GetNorm(0) * 1e-15; // largest (in magnitude) entry
  // we want to call this->changeRows(new_entries, true)
  std::map<int,std::map<int,double> > new_entries;
  const int *rows = structure->GetRowPtr();
  const int *cols = structure->GetKCol();
  int n_rows = structure->GetN_Rows();
  int n_removed = 0; // number of entries to be removed
  // loop over all rows of this matrix
  for(int row=0; row<n_rows; row++)
  {
    new_entries[row]; // empty row
    for(int col = rows[row]; col < rows[row + 1]; col++)
    {
      if(fabs(entries[col]) > tol)
        (new_entries[row])[cols[col]] = entries[col];
    }
    n_removed += rows[row + 1] - rows[row] - new_entries[row].size();
  }
  if(n_removed != 0)
  {
    OutPut("TMatrix::remove_zeros: tol " << tol << "\tn_removed " << n_removed
            << "\tratio " << (double)n_removed/(rows[n_rows]) << endl);
    this->changeRows(new_entries);
  }
  else
    OutPut("TMatrix::remove_zeros: no removable entries\n");
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

void TMatrix::changeRows(std::map<int,std::map<int,double> > entries)
{
  if(entries.size() == 0)
    return; // nothing needs to be done
  
  const int *oldRows = structure->GetRowPtr();
  const int *oldCols = structure->GetKCol();
  
  // find out how many entries there are after all changes are applied, i.e. 
  // how many entries are deleted/created
  int offset = 0;
  for(std::map<int,std::map<int,double> >::iterator it=entries.begin(); 
       it!=entries.end(); ++it)
  {
    int row = it->first;
    offset -= oldRows[row+1]-oldRows[row];// number of entries in old structure
    offset += (it->second).size();        // number of entries in new structure
  }
  
  int n_rows = structure->GetN_Rows();// new number of rows = old number of rows
  // new number of columns = old number of columns
  int n_cols = structure->GetN_Columns(); 
  int n_entries = structure->GetN_Entries() + offset; // new number of entries
  int *columns = new int[n_entries];  // new pointer to columns
  int *rows = new int[n_rows+1];      // new row pointer
  rows[0] = 0;
  
  // create new array to store the entries
  std::vector<double> new_entries(n_entries);
  
  // fill the arrays 'rows', 'columns' and 'new_entries'
  for(int row=0; row<n_rows; row++)
  {
    std::map<int,std::map<int,double> >::iterator it = entries.find(row);
    if(it == entries.end())
    {
      // this row stays unchanged
      // number of (old) entries in this row
      unsigned int n_old_entries = oldRows[row+1] - oldRows[row];
      // copy pointer to columns in this row
      memcpy(columns+rows[row], oldCols+oldRows[row], n_old_entries*SizeOfInt);
      // update row pointer
      rows[row+1] = rows[row] + n_old_entries;
      // copy entries
      memcpy(&new_entries[0]+rows[row], this->GetEntries()+oldRows[row],
             n_old_entries*SizeOfDouble);
    }
    else
    {
      // this row will be replaced
      std::map<int,double> newRow = it->second;
      // loop over all new entries in this row
      int columnIndex=0;
      for(std::map<int,double>::iterator it2 = newRow.begin(); 
          it2 != newRow.end(); ++it2)
      {
        int colInd = it2->first; // column index of new entry
        double entry = it2->second; // value of new entry
        columns[columnIndex+rows[row]] = colInd;
        new_entries[columnIndex+rows[row]] = entry;
        columnIndex++;
      }
      rows[row+1] = rows[row] + newRow.size();
      //if(newRow.size() != columnIndex)
      //  OutPut("ERROR: wrong number of columns in this row "<< newRow.size()
      //      << "\t" << columnIndex << "\t" << row << endl);
    }
  }
  
  // change Structure of this matrix
  this->structure = std::make_shared<TStructure>(n_rows, n_cols, n_entries, 
                                                 columns, rows);
  
  this->entries = new_entries;
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
