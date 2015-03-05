// =======================================================================
// @(#)Structure.C        1.10 11/24/99
// 
// Class:       TStrcuture
//
// Purpose:     build and store a matrix structure
//
// Author:      Gunar Matthies
//
// History:     24.11.97 start implementation (Gunar Matthies)
//
//              04.08.1998 start reimplementation (Gunar Matthies)
//              02.09.1999 adding connections over mortar edges
//                         use __ADD_LINK__ (Volker Behns)
//
// =======================================================================

#include <Structure.h>

/** generate the matrix structure, both spaces are 2D */
TStructure::TStructure()
{
}

/** sort one row [BeginPtr, AfterEndPtr) */
void TStructure::SortRow(int *BeginPtr, int *AfterEndPtr)
{
  int *IPtr, *JPtr, T;

  for(IPtr=BeginPtr;IPtr<AfterEndPtr;IPtr++)
  {
    for(JPtr=IPtr+1;JPtr<AfterEndPtr;JPtr++)
    {
      if( *IPtr > *JPtr )
      {
        T = *IPtr;
        *IPtr = *JPtr;
        *JPtr = T;
      }
    } // endfor JPtr
  } // endfor IPtr
}

/** sort numbers within each row */
void TStructure::Sort()
{
  int i,j,k;
  int end, begin;

  end = 0;
  for(i=0;i<N_Rows;i++)
  {
    begin = end;
    end = RowPtr[i+1];
    SortRow(KCol+begin, KCol+end);
  } // endfor i
}

/** destructor */
TStructure::~TStructure()
{
  delete [] KCol;
  delete [] RowPtr;
}
