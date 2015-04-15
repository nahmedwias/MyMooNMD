// =======================================================================
// @(#)ParSolver.C
//
// Class:      ParSolver 
// Purpose:    Super class for all Parallel Solver classes 
//
// Author:     Sashikumaar Ganesan (19.10.10)
//
// History:    Start of implementation 19.10.10 (Sashikumaar Ganesan)
//
// =======================================================================

#include <ParSolver.h>

// insertion sort,
void TParSolver::Sort(int *Array, int length)
{
 int i, j, a, tmp;
 
 for(i=0; i<length-1; i++)
  {
   for(j=i+1; j<length; j++)
    if(Array[i]>Array[j])
     {
      tmp = Array[i];
      Array[i]=Array[j];
      Array[j]=tmp;
     }
  }

}

int TParSolver::GetIndex(int *Array, int Length, int pos)
{
  int l=0, r=Length, m=(r+l)/2;
  int Mid;

  Mid=Array[m];
  while(Mid!=pos)
  {
    if(Mid<pos)
    { l=m; }
    else
    { r=m; }

    m=(r+l)/2;
    Mid=Array[m];
  }
  return m;
}
