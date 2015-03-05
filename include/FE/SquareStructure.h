// =======================================================================
// @(#)SquareStructure.h        1.6 09/14/99
// 
// Class:       TSquareStructure
//
// Purpose:     build and store a structure for a square matrix
//
// Author:      Gunar Matthies
//
// History:     06.08.1998 start implementation
//
// =======================================================================

#ifndef __SQUARESTRUCTURE__
#define __SQUARESTRUCTURE__

class TSquareStructure
{
  protected:
    /** number of rows */
    int N_Rows;

    /** number of active rows */
    int ActiveBound;

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

    /** ordering of the column entries */
    /** 0 - no special ordering */
    /** 1 - increasing ordering (like used for umfpack) */
    /** 2 - diagonal entry first, then increasing ordering */
    int ColOrder;

    /** sort an integer array */
    void IntSort(int *BeginPtr, int *AfterEndPtr);

  public:
    /** generate the matrix structure, only one space needed */
    TSquareStructure();

    /** destructor: free all used arrays */
    ~TSquareStructure();

    /** generate the matrix structure, all arrays are already defined */
    TSquareStructure(int n, int N_entries, int *col_ptr,
      int *row_ptr);

    /** return number of rows */
    int GetN_Rows()
    { return N_Rows; }

    /** return number of columns */
    int GetN_Columns()
    { return N_Columns; }
    
    /** return ActiveBound */
    int GetActiveBound()
    { return ActiveBound; }

    /** return number of matrix entries */
    int GetN_Entries()
    { return N_Entries; }

    /** return number of matrix entries (hanging nodes part) */
    int GetHangingN_Entries()
    { return HangingN_Entries; }

    /** return array KCol */
    int *GetKCol()
    { return KCol; }

    /** return array HangingKCol */
    int *GetHangingKCol()
    { return HangingKCol; }

    /** return array RowPtr */
    int *GetRowPtr()
    { return RowPtr; }

    /** return array HangingRowPtr */
    int *GetHangingRowPtr()
    { return HangingRowPtr; }

    /** return ordering of columns */
    int GetColOrder()
    { return ColOrder;}

    /** sort column numbers in each row, increasing indices */
    void Sort();

    /** sort column numbers: diag is first element, other numbers are
    increasing */
    void SortDiagFirst();
    
    /** sort column numbers in increasing in all rows including Dirichlet DOF row*/
//     void Sort_ForDirectSolver();
};

#endif
