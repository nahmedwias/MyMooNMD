// =======================================================================
// @(#)SquareMatrix.h        1.3 11/20/98
// 
// Class:       TSquareMatrix
//
// Purpose:     store a square matrix (ansatz = test space)
//
// Author:      Gunar Matthies
//
// History:     10.08.1998 start implementation
//
// =======================================================================

#ifndef __SQUAREMATRIX__
#define __SQUAREMATRIX__

#include <SquareStructure.h>

class TSquareMatrix
{
  protected:
    /** number rof rows */
    int N_Rows;

    /** number columns */
    int N_Columns;

    /** number of matrix entries */
    int N_Entries;

    /** number of matrix entries for hanging nodes data */
    int HangingN_Entries;

    /** in which column is the current entry */
    int *KCol;

    /** KCol for hanging node data */
    int *HangingKCol;

    /** index in KCol where each row starts */
    int *RowPtr;

    /** RowPtr for hanging node data */
    int *HangingRowPtr;
    
    /** matrix elements in an array */
    double *Entries;

    /** bound for active degrees of freedom */
    int ActiveBound;

    /** bound for hanging nodes */
    int HangingBound;

    /** ordering of the column entries */
    /** 0 - no special ordering */
    /** 1 - increasing ordering (like used for umfpack) */
    /** 2 - diagonal entry first, then increasing ordering */
    int ColOrder;
  
  public:
    /** generate the matrix */
    TSquareMatrix(TSquareStructure *structure);

    /** destructor: free Entries array */
    ~TSquareMatrix();

    /** reset matrix entries to zero */
    void Reset();

    /** return number of rows */
    int GetN_Rows()
    { return N_Rows; }

    /** return number of columns */
    int GetN_Columns()
    { return N_Columns; }
    
    /** return number of matrix entries */
    int GetN_Entries()
    { return N_Entries; }

    /** return number of matrix entries for hanging node data */
    int GetHangingN_Entries()
    { return HangingN_Entries; }

    /** return array KCol */
    int *GetKCol()
    { return KCol; }

    /** return array HangingKCol */
    int *GetHangingKCol()
    { return HangingKCol; }

    /** return array HangingRowPtr */
    int *GetHangingRowPtr()
    { return HangingRowPtr; }

    /** return array RowPtr */
    int *GetRowPtr()
    { return RowPtr; }

    /** return matrix entries */
    double *GetEntries()
    { return Entries; }

    /** determine renumbering */
    void ReNumbering(int* &Numbers);

    /** return ActiveBound */
    int GetActiveBound()
    { return ActiveBound; }

    /** return HangingBound */
    int GetHangingBound()
    { return HangingBound; }

    /** return ordering of columns */
    int GetColOrder()
    { return ColOrder;}

    /** write matrix into file */
    int Write(const char *filename);
    
    /** print matrix */
    void Print();
};

#endif
