// =======================================================================
// @(#)Structure.h        1.3 09/14/99
// 
// Class:       TStructure
//
// Purpose:     build and store a matrix structure
//
// Author:      Gunar Matthies
//
// History:     24.11.97 start implementation
//
//              start of reimplementation 26.08.1998 (Gunar Matthies)
//
// =======================================================================

#ifndef __STRUCTURE__
#define __STRUCTURE__

#include <FESpace1D.h>
#include <FESpace2D.h>
#include <FESpace3D.h>

class TStructure
{
  protected:
    /** number of rows */
    int N_Rows;

    /** number columns */
    int N_Columns;

    /** number of matrix entries */
    int N_Entries;

    /** number of matrix entries in hanging nodes part */
    int HangingN_Entries;

    /** in which column is the current entry */
    int *KCol;

    /** in which column is the current entry (hanging nodes part */
    int *HangingKCol;

    /** index in KCol where each row starts */
    int *RowPtr;

    /** index in HangingKCol where each row starts */
    int *HangingRowPtr;
    
    
    // Structure2D
    /** Ansatzspace */
    const TFESpace1D *AnsatzSpace1D;
    const TFESpace2D *AnsatzSpace2D;

    /** Testspace */
    const TFESpace1D *TestSpace1D;
    const TFESpace2D *TestSpace2D;

    int *AnsatzMortarSpaceGlobNo;
    int *TestMortarSpaceGlobNo;

    int *AnsatzNonMortarSpaceGlobNo;
    int *TestNonMortarSpaceGlobNo;
    
    
    // Structure3D
    /** Ansatzspace */
    //TFESpace2D *AnsatzSpace2D;
    const TFESpace3D *AnsatzSpace3D;

    /** Testspace */
    //TFESpace2D *TestSpace2D;
    const TFESpace3D *TestSpace3D;

  public:
    /** generate the matrix structure, both space with 2D collection */
    TStructure();
    
    /** generate a (square) matrix structure, all arrays are already defined */
    TStructure(int n, int N_entries, int *col_ptr, int *row_ptr);

    /** generate the matrix structure, all arrays are already defined */
    TStructure(int nRows, int nCols, int N_entries, int *col_ptr, int *row_ptr);
    
    /** Generates an empty nRows*nCols Structure for a Zero-Matrix */
    TStructure(int nRows, int nCols);
    
    
    /** generate the matrix Structure2D, both space with 2D collection */
    TStructure(const TFESpace2D *testspace, const TFESpace2D *ansatzspace);
    
    /** generate the matrix structure, both spaces are 2D */
    /** both spaces are defined on different grids */
    TStructure(TFESpace2D *testspace, int test_level, TFESpace2D *ansatzspace,
               int ansatz_level);
#ifdef __MORTAR__
    /** generate the matrix Structure2D, one space with 1D and the other
     with 2D collection */
    TStructure2D(TFESpace1D *testspace, TFESpace2D *ansatzspace);
#endif
    
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
    
    /** generate the matrix Structure3D, both space with 3D collection */
    TStructure(const TFESpace3D *testspace, const TFESpace3D *ansatzspace);
    
    
    
    
    
    
    
    // structure2d
    /** return AnsatzSpace */
    const TFESpace2D *GetAnsatzSpace2D() const
    {
      return AnsatzSpace2D;
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
    
    /** return TestSpace */
    const TFESpace2D *GetTestSpace2D() const
    {
      return TestSpace2D;
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

    
    
    
    //structure 3d
    
    
    
    
    /** @brief copy constructor */
    TStructure(const TStructure&);

    /** destructor: free all used arrays */
    ~TStructure();

    /** return number of rows */
    int GetN_Rows() const
    { return N_Rows; }

    /** return number of columns */
    int GetN_Columns() const
    { return N_Columns; }
    
    /** return number of matrix entries */
    int GetN_Entries() const
    { return N_Entries; }

    /** return number of matrix entries (hanging nodes part) */
    int GetHangingN_Entries() const
    { return HangingN_Entries; }

    /** return array KCol */
    int *GetKCol() const
    { return KCol; }

    /** return array HangingKCol */
    int *GetHangingKCol() const
    { return HangingKCol; }

    /** return array RowPtr */
    int *GetRowPtr() const
    { return RowPtr; }

    /** return array HangingRowPtr */
    int *GetHangingRowPtr() const
    { return HangingRowPtr; }
    
    /** @brief set member variables. Careful, this can produce inconsistencies! */
    void setN_Rows(int n) { N_Rows = n; }
    void setN_Columns(int n) { N_Columns = n; }
    void setN_Entries(int n) { N_Entries = n; }
    void setKCol(int * p) { KCol = p; }
    void setRowPtr(int * p) { RowPtr = p; }

    /** sort one row */
    void SortRow(int *BeginPtr, int *AfterEndPtr);
    
    /** sort rows */
    void Sort();
    
    
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
    
    /** return a new structure for a transposed matrix 
     * If this is an object of a derived class (e.g. TStructure2D, 
     * TSquareStructure), then the number of active degrees of freedom is not 
     * taken into account. The returned TMatrix is really the algebraic 
     * transposed matrix.
     * */
    TStructure* GetTransposed();
    
    /**
     * @brief return a structure for the matrix-matrix-product A*B
     * 
     * if A and B are matrices with structures 'strucA' and 'strucB', this
     * function computes a structure for the product C = A*B
     * 
     * @param strucA structure of left factor
     * @param strucB structure of right factor
     */
    friend TStructure* get_product_structure(TStructure const * const strucA,
                                             TStructure const * const strucB);
    
    /** @brief Comparision Operator 
     * 
     * It is not explicitly checked if the arrays are the same, only the 
     * integers are compared.
     */
    friend bool operator==(const TStructure &lhs, const TStructure &rhs);
    friend bool operator!=(const TStructure &lhs, const TStructure &rhs);
};

#endif
