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

#include <Matrix.h>
#include <Structure.h>

class TSquareMatrix : public TMatrix
{
  protected:
    /** bound for hanging nodes 
     * @todo is this structure->HangingN_Entries ?? */
    int HangingBound;

    /** ordering of the column entries */
    /** 0 - no special ordering */
    /** 1 - increasing ordering (like used for umfpack) */
    /** 2 - diagonal entry first, then increasing ordering */
    int ColOrder;
  
    /** generate the matrix, called from derived classes */
    TSquareMatrix(std::shared_ptr<TStructure> structure);
  public:

    /** destructor: free Entries array */
    ~TSquareMatrix() = default;
    
    TSquareMatrix(const TSquareMatrix & m) = default;

    /** reset all entries in active rows */
    void ResetActive();
    
    /** @brief set zeros in nonactive rows. 
     * 
     * This is e.g. for the off-diagonal blocks in a Stokes matrix 
     */
    void resetNonActive();

    /** determine renumbering */
    void ReNumbering(int* &Numbers) const;

    /** return ActiveBound */
    int GetActiveBound() const
    { return structure->GetActiveBound(); }
    
    /** return ordering of columns */
    int GetColOrder() const
    { return structure->GetColOrder(); }
    
    /** write matrix into file */
    int Write(const char *filename);
    
    /** print matrix */
    void Print();
};

#endif
