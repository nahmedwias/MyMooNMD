/**
 * @class FEMatrix
 * @brief matrix which represents the mapping between two finite element spaces
 * 
 * 
 * This class essentially implements constructors which take finite element
 * spaces. Also it saves which spaces were used during construction. Otherwise
 * it behaves just like its base class Matrix.
 * 
 * Note that there is no constructor taking a TStructure together with the 
 * spaces it was created with. If you want to create multiple matrices with the
 * same structure object, use the copy constructor which will do just that.
 * 
 * @ruleof0
 * 
 * @todo write a test for this class
 */

#ifndef __FEMATRIX__
#define __FEMATRIX__

#include <FESpace1D.h>
#include <FESpace2D.h>
#include <FESpace3D.h>
#include <Matrix.h>

class FEMatrix : public TMatrix
{
  public:
    /// @brief there is no default constructor, this makes no sense
    FEMatrix() = delete;
    
    /// @name construct square structures using one finite element space
    /// @brief ansatz and test space is the same 
    //@{
    FEMatrix(const TFESpace1D * space);
    FEMatrix(const TFESpace2D * space);
    #ifdef __3D__
    FEMatrix(const TFESpace3D * space);
    #endif // 3D
    //@}
    
    /// @name construct rectangular structures using two finite element spaces
    /// @brief test and ansatz space are possibly different
    /// 
    /// A matrix using this structure represents a linear map from the ansatz 
    /// to the test space.
    /// @param[in] is_empty If true, this TMatrix is created with an empty structure.
    //@{
    FEMatrix(const TFESpace2D * testspace, const TFESpace2D * ansatzspace,  bool is_empty = false);
    #ifdef __3D__
    FEMatrix(const TFESpace3D * testspace, const TFESpace3D * ansatzspace);
    #endif // 3D
    //@}
    
    //! Default copy constructor. Performs shallow copy.
    FEMatrix(const FEMatrix&) = default;
    
    //! Default move constructor.
    FEMatrix(FEMatrix&&) = default;
    
    //! Default copy assignment operator. Performs shallow copy.
    FEMatrix& operator=(const FEMatrix&) = default;
    
    //! Default move assignment operator
    FEMatrix& operator=(FEMatrix&&) = default;
    
    //! Default destructor.
    ~FEMatrix() = default;
    
    
    /** @brief reset all entries in active rows to zero */
    void resetActive();

    /** @brief set zeros in nonactive rows. 
     * 
     * This is e.g. for the off-diagonal blocks in a Stokes matrix 
     */
    void resetNonActive();

    /** @brief scale this matrix by a factor
     * 
     * Only rows corresponding to active d.o.f are scaled. Other rows remain
     * unscaled.
     */
    void scaleActive(double factor = 1.0);

    /** @brief adding a scaled matrix to this matrix
     * 
     * This is only done for those rows which correspond to active degrees of 
     * freedom.
     * 
     * The summation is index-wise, i.e. A(i,j) += factor*m(i.j), where A is 
     * this matrix. 
     * 
     * Note that this only works if the sparsity structure is the same for this
     * matrix and m.
     */
    void addActive(const FEMatrix& m, double factor = 1.0);
    
    /** @brief compute y = y + a * Ax 
     *
     * add the matrix-vector product "Ax", scaled by "a", to y: only active
     * "A" is this matrix.
     * 
     * This function can be used to compute the residual r = b - Ax.
     *
     * @param x the vector which is multiplied by this matrix
     * @param y result of matrix-vector-multiplication and scaling
     * @param factor optional scaling   factor, default to 1.0
     */
    void multiplyActive(const double *x, double *y, double factor = 1.0)
      const;
    
    
    /** @brief return the number of active rows */
    int GetActiveBound() const;
    
    
    /** @brief return 1D test space */
    const TFESpace1D *GetTestSpace1D() const;
    
    /** @brief return 1D ansatz space */
    const TFESpace1D *GetAnsatzSpace1D() const;
    
    /** @brief return 2D test space */
    const TFESpace2D *GetTestSpace2D() const;
    
    /** @brief return 2D ansatz space */
    const TFESpace2D *GetAnsatzSpace2D() const;
    
    #ifdef __3D__
    /** @brief return 3D test space */
    const TFESpace3D *GetTestSpace3D() const;
    
    /** @brief return 3D ansatz space */
    const TFESpace3D *GetAnsatzSpace3D() const;
    #endif // 3D
    
    /** @brief return test space */
    const TFESpace *GetTestSpace() const;
    
    /** @brief return ansatz space */
    const TFESpace *GetAnsatzSpace() const;
    
    /** return FESpace */
    const TFESpace1D *GetFESpace1D() const;
    /** return FESpace */
    const TFESpace2D *GetFESpace2D() const;
    #ifdef __3D__
    /** return FESpace */
    const TFESpace3D *GetFESpace3D() const;
    #endif // 3D
    
  private:
    /// @name ansatz spaces
    /// @brief the ansatz space (pre-image space)
    /// @details Exactly one of these pointers is not a nullptr.
    /// @todo make this a share_ptr
    //@{
    const TFESpace1D* AnsatzSpace1D;
    const TFESpace2D* AnsatzSpace2D;
    const TFESpace3D* AnsatzSpace3D;
    //@}
    
    /// @name test spaces
    /// @brief the test space (image space)
    /// @details Exactly one of these pointers is not a nullptr.
    /// @todo make this a share_ptr
    //@{
    const TFESpace1D* TestSpace1D;
    const TFESpace2D* TestSpace2D;
    const TFESpace3D* TestSpace3D;
    //@}
};

#endif //__FEMATRIX__
