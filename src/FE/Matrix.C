// =======================================================================
// @(#)Matrix.C        1.2 11/20/98
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

#include <Matrix.h>
#include <string.h>
#include <Constants.h>

#include <MooNMD_Io.h>

TMatrix::TMatrix(TStructure *structure)
{
}

void TMatrix::Reset()
{
  memset(Entries, 0, N_Entries*SizeOfDouble);
}

TMatrix::~TMatrix()
{
  delete Entries;
}

/** write matrix into file */
int TMatrix::Write(const char *filename)
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
