// =======================================================================
// @(#)SquareStructure.C        1.6 09/17/99
// 
// Class:       TSquareStructure
//
// Purpose:     build and store a structure for a square matrix in 2d
//
// Author:      Gunar Matthies
//
// History:     06.08.1998 start implementation
//
//              02.09.1999 adding connections over mortar edges
//                         use __ADD_LINK__ (Volker Behns)
//
// =======================================================================

#include <SquareStructure.h>
#include <Database.h>
#include <string.h>

/** generate the matrix structure, both space are 2D */
TSquareStructure::TSquareStructure()
{
} 

/** generate the matrix structure, all arrays are already defined */
TSquareStructure::TSquareStructure(int n, int N_entries, int *col_ptr,
				   int *row_ptr)
{
}

TSquareStructure::~TSquareStructure()
{
  delete [] KCol;
  delete [] RowPtr;
}

/** sort an integer array [BeginPtr, AfterEndPtr) */
void TSquareStructure::IntSort(int *BeginPtr, int *AfterEndPtr)
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

/** sort column numbers: numbers are increasing */
/** this is the ordering used by umfpack */

void TSquareStructure::Sort()
{
  int i,j;
  int end, begin;

  for(i=0;i<ActiveBound;i++)
  {
    begin = RowPtr[i];
    end = RowPtr[i+1];
    IntSort(KCol+begin, KCol+end);
  } // endfor i
  ColOrder = 1;

  // in this case umfpack needs an ordering before the call
  if (HangingN_Entries>0)
      ColOrder = -1;

}

/** sort column numbers: diag is first element, other numbers are
    increasing */
void TSquareStructure::SortDiagFirst()
{
  int i,j,k;
  int end, begin;

  end = 0;
  for(i=0;i<ActiveBound;i++)
  {
    begin = end;
    end = RowPtr[i+1];
    k = KCol[begin];
    for(j=RowPtr[i];j<end;j++)
    {
      if(KCol[j] == i)
      {
        // diag entry
        KCol[begin] = i;
        KCol[j] = k;
        break;
      } // endif
    } // endfor j
    IntSort(KCol+begin+1, KCol+end);
  } // endfor i
  
  ColOrder = 2;
}

/** sort column numbers in increasing in all rows including Dirichlet DOF row*/
// void TSquareStructure::Sort_ForDirectSolver()
// {
//   int i,j,k;
//   int end, begin;
// 
//   end = 0;
//   for(i=0;i<N_Rows;i++)
//   {
//     begin = end;
//     end = RowPtr[i+1];
//     IntSort(KCol+begin, KCol+end);
//   } // endfor i
// }
// 
