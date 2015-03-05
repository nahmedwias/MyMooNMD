// =======================================================================
// @(#)SquareMatrix.C        1.2 11/20/98
// 
// Class:       TSquareMatrix
//
// Purpose:     store a square matrix (ansatz = test space) in 
//
// Author:      Gunar Matthies
//
// History:     10.08.1998 start implementation
//
// =======================================================================

#include <SquareMatrix.h>
#include <Constants.h>
#include <string.h>
#include <MooNMD_Io.h>

TSquareMatrix::TSquareMatrix(TSquareStructure *squarestructure)
{
}

void TSquareMatrix::Reset()
{
  memset(Entries, 0, N_Entries*SizeOfDouble);
}

TSquareMatrix::~TSquareMatrix()
{
  if (Entries) delete [] Entries;
}

// find a renumbering of the DOFs from the matrix entries
void TSquareMatrix::ReNumbering(int* &Numbers)
{
  int i,j,k,l,N_;
  int *Inflow, N_Active;
  int begin, end, beginJ, endJ;
  int Index, State, CurrentNumber;
  double aij, aji;

  N_Active = ActiveBound;
  Inflow = new int[N_Active];
  memset(Inflow, 0, N_Active*SizeOfInt);
  Numbers = new int[N_Active];
  memset(Numbers, -1, N_Active*SizeOfInt);

  for(i=0;i<N_Active;i++)
  {
    begin = RowPtr[i];
    end = RowPtr[i+1];
    for(j=begin+1;j<end;j++)
    {
      k = KCol[j];
      if( (k>i) && (k<N_Active))
      {
        aij = Entries[j];
        beginJ = RowPtr[k];
        endJ = RowPtr[k+1];
        for(l=beginJ+1;l<endJ;l++)
        {
          if( KCol[l] == i )
          {
            aji = Entries[l];
            if(aij>aji) Inflow[k]++;
            if(aji>aij) Inflow[i]++;
          }
        }
      } // endif
    } // endfor j
  } // endfor i

  Index = 0; 
  for(i=0;i<N_Active;i++)
  {
    cout << i << "    " << Inflow[i] << endl;
    if(Inflow[i] == 0)
    {
      Numbers[Index] = i;
      Index++;
    }
  } // endfor i

  State = 0;
  while(State<N_Active)
  {
    CurrentNumber = Numbers[State];
    begin = RowPtr[CurrentNumber];
    end = RowPtr[CurrentNumber+1];
    for(j=begin+1;j<end;j++)
    {
      k = KCol[j];
      if(k<N_Active)
      {
        aij = Entries[j];
        beginJ = RowPtr[k];
        endJ = RowPtr[k+1];
        for(l=beginJ+1;l<endJ;l++)
        {
          if( KCol[l] == CurrentNumber )
          {
            aji = Entries[l];
            if(aij>aji)
            {
              if(Inflow[k] > 0)
              {
                Inflow[k]--;
                if(Inflow[k] == 0)
                {
                  Numbers[Index] = k;
                  Index++;
                }
              }
            }
          }
        } // endfor l
      } // endif
    } // endfor j
    State++;
  } // endwhile

  for(i=0;i<N_Active;i++)
  {
    cout << i << "   " << Numbers[i] << endl;
  }

  delete Inflow;
}

/** write matrix into file */
int TSquareMatrix::Write(const char *filename)
{
  int header[3];

  std::ofstream dat(filename);
  if(!dat)
  {
    cerr << "cannot open file '" << filename << "' for output" << endl;
    return -1;
  }

  header[0] = N_Rows;
  header[1] = N_Columns;
  header[2] = N_Entries;

  dat.write((char *)header, sizeof(int)*3);
  dat.write((char *)RowPtr, sizeof(int)*(N_Rows+1));
  dat.write((char *)KCol, sizeof(int)*N_Entries);
  dat.write((char *)Entries, sizeof(double)*N_Entries);

  dat.close();
  
  return 0;
}

void TSquareMatrix::Print()
{
  int begin, end, pos=0;
  
  for (int i=0;i<N_Rows;++i)
  {
    begin = RowPtr[i];
    end   = RowPtr[i+1];
    
    for (int j=begin;j<end;++j)
    {
      cout << "a(" << i << "," << KCol[pos] << ") = " << Entries[pos] << endl;
      ++pos;
    }
  }
}

