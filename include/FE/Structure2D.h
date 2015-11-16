// =======================================================================
// @(#)Structure2D.h        1.3 09/14/99
// 
// Class:       TStructure2D
//
// Purpose:     build and store a matrix Structure2D
//
// Author:      Gunar Matthies
//
// History:     24.11.97 start implementation
//
//              start of reimplementation 26.08.1998 (Gunar Matthies)
//
// =======================================================================

#ifndef __STRUCTURE2D__
#define __STRUCTURE2D__

#include <FESpace1D.h>
#include <FESpace2D.h>
#include <Structure.h>

class TStructure2D : public TStructure
{
  protected:

  public:
    /** generate the matrix Structure2D, both space with 2D collection */
    TStructure2D(const TFESpace2D *testspace, const TFESpace2D *ansatzspace)
     : TStructure(testspace, ansatzspace) {}


    /** generate the matrix structure, both spaces are 2D */
    /** both spaces are defined on different grids */
    TStructure2D(TFESpace2D *testspace, int test_level, 
                 TFESpace2D *ansatzspace, int ansatz_level)
     : TStructure(testspace, test_level, ansatzspace, ansatz_level) {}
    #ifdef __MORTAR__
    /** generate the matrix Structure2D, one space with 1D and the other
        with 2D collection */
    TStructure2D(TFESpace1D *testspace, TFESpace2D *ansatzspace)
     : TStructure(testspace, ansatzspace) {}
    #endif

    /** generate the matrix Structure2D, one space with 1D and the other
        with 2D collection */
     TStructure2D(TFESpace1D *testspace, TFESpace2D *ansatzspace, int **ansatzcelljoints)
      : TStructure(testspace, ansatzspace, ansatzcelljoints) {}
     
     
     TStructure2D(TFESpace2D *testspace, TFESpace1D *ansatzspace, TNonMortarData *NonMortarFEData)
     : TStructure(testspace, ansatzspace, NonMortarFEData) {}

};

#endif
