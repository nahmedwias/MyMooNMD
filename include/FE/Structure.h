/**
 * @class TStructure
 * @brief represent the sparsity structure of a matrix
 * 
 * Matrices in ParMooN are stored in the compressed row storage (CRS) format.
 * The structure essentially stores the number of rows (TStructure::nRows), the 
 * number of columns (TStructure::nColumns), the number of entries 
 * (TStructure::nEntries) and two vectors (TStructure::rows and 
 * TStructure::columns). 
 * 
 * The first vector TStructure::rows has length TStructure::nRows+1, indicating
 * how many entries there are in all previous rows. It always is
 *     TStructure::rows[0] = 0;
 *     TStructure::rows[TStructure::nRows] = TStructure::nEntries;
 * The number of entries in a given row `i` is given by the difference 
 *     unsigned int nEntriesInRow = TStructure::rows[i+1] - TStructure::rows[i];
 * 
 * The second vector TStructure::columns has length TStructure::nEntries and 
 * represents the column of each entry.  
 * 
 * To get the `j`-th entry in the `i`-th row (not the `j`-th entry overall) it
 * is
 *     unsigned int indexOfEntry = TStructure::rows[i] + j; 
 *     unsigned int columnOfEntry = TStructure::columns[indexOfEntry];
 * 
 * Now the entry (`i`, `columnOfEntry`) in a matrix with this structure is the
 * `indexOfEntry`-th entry in the list of entries. 
 * 
 * \todo How do we handle Hanging nodes? We need a test case!
 */

#ifndef __STRUCTURE__
#define __STRUCTURE__

#include <FESpace1D.h>
#include <FESpace2D.h>
#include <FESpace3D.h>
#include <memory>

class TStructure
{
  private:
    /** @brief number of rows */
    int nRows;

    /** @brief number columns */
    int nColumns;

    /** @brief number of matrix entries */
    int nEntries;

    /** @brief vector of column information
     * 
     * This vector has length TStructure::nEntries and stores the column index
     * for each entry.
     */
    int *columns;
    
    /** @brief vector storing the number of entries in all previous rows
     * 
     * This vector tells you the global indices of the entries for each row.
     * The `j`-th entry in the `i`-th row is at the `TStructure[i]+j`-th 
     * position in the global entries vector. So this also gives you the right 
     * index if you are interested in the column of that entry, 
     * `TStructure::columns[TStructure[i]+j]`.
     */
    int *rows;

    /** @brief number of active rows
     * 
     * If this structure has been created using finite element spaces where the
     * test space has a Dirichlet (non-active) boundary, then this number is
     * smaller than TStructure::nRows. Otherwise it is the number of rows. In 
     * ParMooN the non-active degrees of freedom (those which are on a 
     * Dirichlet boundary) are always the last rows. These rows typically only
     * have one entry (on the diagonal).  
     */
    int ActiveBound;

    /** @brief ordering of the column entries
     * 
     * 0 - no special ordering 
     * 1 - increasing ordering (like used for umfpack)
     * 2 - diagonal entry first, then increasing ordering (for AMG)
     */
    int ColOrder;

    // entries related to hanging nodes
    /** @brief number of matrix entries in hanging nodes part */
    int nHangingEntries;

    /** @brief in which column is the current entry (hanging nodes part) */
    int *hangingColums;

    /** @brief index in hangingColums where each row starts */
    int *HangingRows;
    
    
    /// @name ansatz spaces
    /// @brief the ansatz space (pre-image space)
    /// @details Exactly one of these pointers is not a nullptr.
    //@{
    const TFESpace1D* AnsatzSpace1D;
    const TFESpace2D* AnsatzSpace2D;
    const TFESpace3D* AnsatzSpace3D;
    //@}

    /// @name test spaces
    /// @brief the test space (image space)
    /// @details Exactly one of these pointers is not a nullptr.
    //@{
    const TFESpace1D* TestSpace1D;
    const TFESpace2D* TestSpace2D;
    const TFESpace3D* TestSpace3D;
    //@}
    
    
    /** @brief sort all rows in increasing order */
    void Sort();
    
    /** @brief sort one row to increasing order */
    void SortRow(int *BeginPtr, int *AfterEndPtr);
    
    /** @brief sort column numbers
     * 
     * diag is first element, other numbers are increasing
     */
    void SortDiagFirst();
    
  public:
    /** @brief default constructor sets everything to 0/nullptr */
    TStructure();
    
    /// @name construct square structures using one finite element space
    /// @brief ansatz and test space is the same 
    //@{
    TStructure(const TFESpace1D * space);
    TStructure(const TFESpace2D * Space);
    TStructure(const TFESpace3D * space);
    //@}
    
    /// @name construct rectangular structures using two finite element spaces
    /// @brief test and ansatz space are possibly different
    /// 
    /// A matrix using this structure represents a linear map from the ansatz 
    /// to the test space.
    //@{
    TStructure(const TFESpace2D * testspace, const TFESpace2D * ansatzspace);
    TStructure(const TFESpace3D * testspace, const TFESpace3D * ansatzspace);
    //@}
    
    /**
     * @brief construct rectangular structures using two finite element spaces 
     * defined on different grids
     * 
     * These grids must be part of the same hierarchy and the corresponding 
     * levels must be given. 
     */
    TStructure(const TFESpace2D * testspace, int test_level, 
               const TFESpace2D * ansatzspace, int ansatz_level);
    

    /**
     * @brief generate a square structure, all arrays are already defined
     * 
     * Note that the pointers will be deleted from within this class. No deep 
     * copy is done.
     * 
     * @param n number of rows/columns
     * @param N_entries number of entries in this structure
     * @param col_ptr the new TStructure::columns vector
     * @param row_ptr the new TStructure::rows vector
     */
    TStructure(int n, int N_entries, int *col_ptr, int *row_ptr);

    /**
     * @brief generate a structure, all arrays are already defined
     * 
     * Note that the pointers will be deleted from within this class. No deep 
     * copy is done.
     * 
     * @param nRows number of rows
     * @param nCols number of columns
     * @param N_entries number of entries in this structure
     * @param col_ptr the new TStructure::columns vector
     * @param row_ptr the new TStructure::rows vector
     */
    TStructure(int nRows, int nCols, int N_entries, int *col_ptr, int *row_ptr);

    /** @brief Generates an empty `n`*`n` Structure with no entries */
    explicit TStructure(int n);
    
    /** @brief Generates an empty `nRows`*`nCols` Structure with no entries
     */
    TStructure(int nRows, int nCols);

    
    /** @brief return if this structure is square */
    bool isSquare() const
    { return nRows == nColumns; }
    
    /** return TestSpace */
    const TFESpace1D *GetTestSpace1D() const
    {
      return TestSpace1D;
    }
    
    /** return AnsatzSpace */
    const TFESpace1D *GetAnsatzSpace1D() const
    {
      return AnsatzSpace1D;
    }
        
    /** return TestSpace */
    const TFESpace2D *GetTestSpace2D() const
    {
      return TestSpace2D;
    }

    /** return AnsatzSpace */
    const TFESpace2D *GetAnsatzSpace2D() const
    {
      return AnsatzSpace2D;
    }
    
    /** return TestSpace */
    const TFESpace3D *GetTestSpace3D() const
    {
      return TestSpace3D;
    }
    
    /** return AnsatzSpace */
    const TFESpace3D *GetAnsatzSpace3D() const
    {
      return AnsatzSpace3D;
    }
    
    
    /** return TestSpace */
    const TFESpace *GetTestSpace() const
    {
      if(TestSpace1D)
      {
        return TestSpace1D;
      }
      else if(TestSpace2D)
      {
        return TestSpace2D;
      }
      else
      {
        return TestSpace3D;
      }
    }
    
    const TFESpace *GetAnsatzSpace() const
    {
      if(AnsatzSpace1D)
      {
        return AnsatzSpace1D;
      }
      else if(AnsatzSpace2D)
      {
        return AnsatzSpace2D;
      }
      else
      {
        return AnsatzSpace3D;
      }
    }
    
    
    /** return FESpace */
    const TFESpace1D *GetFESpace1D() const
    { if(this->isSquare()) return TestSpace1D;
    else ErrThrow("accessing FESpace for non-square matrix, but which one?");}
    /** return FESpace */
    const TFESpace2D *GetFESpace2D() const
    { if(this->isSquare()) return TestSpace2D;
        else ErrThrow("accessing FESpace for non-square matrix, but which one?");}
    /** return FESpace */
    const TFESpace3D *GetFESpace3D() const
    { if(this->isSquare()) return TestSpace3D;
        else ErrThrow("accessing FESpace for non-square matrix, but which one?");}
    
    
    /**
     * @brief return number of active rows
     * 
     * @return between 0 and TStructure::nRows 
     */
    int GetActiveBound() const
    {
      return ActiveBound;
    }
    
    unsigned int getNActiveEntries() const;
    
    /** @brief return ordering of columns, see also TStructure::colOrder */
    int GetColOrder() const
    {
      return ColOrder;
    }
    
    
    /** @brief return number of rows */
    int GetN_Rows() const
    { return nRows; }

    /** @brief return number of columns */
    int GetN_Columns() const
    { return nColumns; }
    
    /** @brief return number of entries */
    int GetN_Entries() const
    { return nEntries; }

    /** @brief return number of entries (hanging nodes part) */
    int GetHangingN_Entries() const
    { return nHangingEntries; }

    /** @brief return vector columns */
    int *GetKCol() const
    { return columns; }

    /** @brief return array hangingColums */
    int *GetHangingKCol() const
    { return hangingColums; }

    /** return array RowPtr */
    int *GetRowPtr() const
    { return rows; }

    /** return array HangingRowPtr */
    int *GetHangingRowPtr() const
    { return HangingRows; }
    
    /**
     * @brief find the index of a given entry
     * 
     * If the (i,j)-th entry is not in the sparsity pattern, -1 is returned. 
     * This is how this function can be used to check whether an entry is in the
     * sparsity pattern.
     * 
     * @param i row of entry to check
     * @param j column of entry to check
     */ 
    int index_of_entry(const int i, const int j) const;
    
    /** 
     * @brief return a new structure for a transposed matrix
     * 
     * If this structure has been created using finite element spaces, then 
     * this might be different from directly constructing a structure with 
     * test and anstz space exchanged. This can happen due to the non-active 
     * degrees of freedom.
     * 
     * This function returns a true (algebraic) transposed structure.
     */
    std::shared_ptr<TStructure> GetTransposed() const;
    
    
    /** @brief copy constructor */
    TStructure(const TStructure&);

    /** @brief destructor: delete all used arrays */
    ~TStructure();
    
    /**
     * @brief return a structure for the matrix-matrix-product A*B
     * 
     * if A and B are matrices with structures 'strucA' and 'strucB', this
     * function computes a structure for the product C = A*B
     * 
     * @param strucA structure of left factor
     * @param strucB structure of right factor
     */
    friend std::shared_ptr<TStructure> get_product_structure(
        TStructure const & strucA, TStructure const & strucB);
    
    /** @brief Comparision Operator 
     * 
     * It is not explicitly checked if the arrays are the same, only the 
     * integers are compared.
     */
    friend bool operator==(const TStructure &lhs, const TStructure &rhs);
    friend bool operator!=(const TStructure &lhs, const TStructure &rhs);
    
    
#ifdef __MORTAR__ // \todo can this mortar stuff be removed?
  protected:
    int *AnsatzMortarSpaceGlobNo;
    int *TestMortarSpaceGlobNo;

    int *AnsatzNonMortarSpaceGlobNo;
    int *TestNonMortarSpaceGlobNo;
  public:
    /** generate the matrix Structure2D, one space with 1D and the other
     with 2D collection */
    TStructure(TFESpace1D *testspace, TFESpace2D *ansatzspace);
    
    /** generate the matrix Structure2D, one space with 1D and the other
     with 2D collection */
    TStructure(TFESpace1D *testspace, TFESpace2D *ansatzspace,
        int **ansatzcelljoints);

    /** generate the matrix Structure2D, one space with 1D and the other
     with 2D collection */
    TStructure(TFESpace1D *testspace, TFESpace2D *ansatzspace,
        TNonMortarData *NonMortarFEData);

    TStructure(TFESpace2D *testspace, TFESpace1D *ansatzspace,
        TNonMortarData *NonMortarFEData);
#endif
    
};

#endif
