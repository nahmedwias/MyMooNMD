// =======================================================================
// @(#)Matrix.h        1.2 11/20/98
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

#ifndef __MATRIX__
#define __MATRIX__

#include <Structure.h>

class TMatrix
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

  public:
    /** generate the matrix */
    TMatrix(TStructure *structure);

    /** destructor: free Entries array */
    ~TMatrix();

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

    /** write matrix into file */
    int Write(const char *filename);
};

#endif
