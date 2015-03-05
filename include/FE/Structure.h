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

  public:
    /** generate the matrix structure, both space with 2D collection */
    TStructure();

    /** destructor: free all used arrays */
    ~TStructure();

    /** return number of rows */
    int GetN_Rows()
    { return N_Rows; }

    /** return number of columns */
    int GetN_Columns()
    { return N_Columns; }
    
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

    /** sort one row */
    void SortRow(int *BeginPtr, int *AfterEndPtr);
    
    /** sort rows */
    void Sort();

};

#endif
