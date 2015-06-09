/** ************************************************************************ 
*
* @class     SystemMat2D
* @brief     stores the information of a 2D system matrix
* 
* This is a base class. Derived classes for differnet problem types exist and 
* should be used rather than using this class directly.
* 
* @author    Ulrich Wilbrandt,
* @date      17.03.15
* 
* @history Ulrich Wilbrandt, June 2015, added some linear algebra routines
* 
* ************************************************************************  */
#ifndef __SYSTEMMAT2D__
#define __SYSTEMMAT2D__

#include <SquareMatrix2D.h>
#include <Matrix2D.h>
#include <vector>

class SystemMat2D
{
  protected:

    /** @brief vector of finite element spaces 
     * 
     * Each space appears exactly once in this vector. Typically this vector has
     * size 1 (scalar) or two (saddle-point problems).
     */
    std::vector<TFESpace2D*> fe_spaces;
    
    /** @brief vector of the square matrices */
    std::vector<TSquareMatrix2D*> sq_matrices;
  
    /** @brief vector of the rectangular matrices */
    std::vector<TMatrix2D*> rect_matrices;
    
    /** @brief method for resudual calculation */
    DefectProc *defect;
    
  public:
    /** @brief standard constructor
     * 
     * This initializes the vectors with nullptr. The vectors are filled within
     * the derived classes. 
     */
    SystemMat2D(unsigned int n_spaces, unsigned int n_sq_matrices, 
                unsigned int n_rect_matrices)
     : fe_spaces(n_spaces, nullptr), sq_matrices(n_sq_matrices, nullptr),
       rect_matrices(n_rect_matrices, nullptr), defect(nullptr)
    {};
    
    virtual ~SystemMat2D();
    
    /** @brief adding a scaled matrix to this matrix
     * 
     * The summation is index-wise, i.e. M(i,j) += factor*A(i.j), where M is 
     * this matrix.
     * 
     * Only the rows corresponding to active degrees of freedom are scaled, all
     * other rows are left unchanged.
     * 
     * Note that this only works if the sparsity structure is the same for this
     * matrix and A.
     */
    void add_active(const SystemMat2D &A, double factor = 1.0);
    
    /** @brief adding a scaled matrix to this matrix
     * 
     * The summation is index-wise, i.e. M(i,j) += factor*A(i.j), where M is 
     * this matrix. 
     * 
     * Note that this only works if the sparsity structure is the same for this
     * matrix and A.
     */
    void add(const SystemMat2D &A, double factor = 1.0);

    /** @brief scale this matrix by a factor
     * 
     * Only rows corresponding to active degrees of freedom are scaled. Other 
     * rows remain unchanged.
     */
    void scale_active(double factor);
    
    /** @brief scale this matrix
     * 
     * That means for each submatrix all entries are scaled.
     */
    void scale(double factor);
    
    /** @brief compute y = factor* A*x 
     *
     * write the matrix-vector product "Ax" scaled by a factor to y. 
     * "A" is this matrix. Both "A" and "x" remain unchanged. If the factor is
     * 0.0, the vector y is only reset without performing the actual 
     * multiplication.
     *
     * @param x the vector which is multiplied by this matrix
     * @param y result of matrix-vector-multiplication and scaling
     * @param factor optional scaling factor, default to 1.0
     */
    virtual void apply(const double *x, double *y, double factor = 1.0) const;
    
    /** @brief compute y = y + a * Ax 
     *
     * add the matrix-vector product "Ax", scaled by "a", to y.
     * "A" is this matrix. Both "A" and "x" remain unchanged. If the factor is
     * 0.0, nothing is done.
     * 
     * This function can be used to compute the residual r = b - Ax.
     *
     * @param x the vector which is multiplied by this matrix
     * @param y result of matrix-vector-multiplication and scaling
     * @param factor optional scaling   factor, default to 1.0
     */
    virtual void apply_scaled_add(const double *x, double *y, 
                                  double factor = 1.0) const;
    
    // some getters
    const TFESpace2D* get_space(unsigned int i = 0) const 
    { return fe_spaces.at(i); }
    
    const TSquareMatrix2D* get_square_matrix(unsigned int i = 0) const
    { return sq_matrices.at(i); }
    TSquareMatrix2D* get_square_matrix(unsigned int i = 0) 
    { return sq_matrices.at(i); }
    
    const TMatrix2D* get_rectangular_matrix(unsigned int i = 0) const
    { return rect_matrices.at(i); }
    TMatrix2D* get_rectangular_matrix(unsigned int i = 0)
    { return rect_matrices.at(i); }
    
    unsigned int n_blocks() const; // total number of blocks
    virtual unsigned int n_rows() const; // number of block rows
    virtual unsigned int n_cols() const; // number of block columns
    virtual unsigned int n_total_rows() const;//total number of rows(>n_rows)
    virtual unsigned int n_total_cols() const;//total number of columns(>n_cols)
    unsigned int n_total_entries() const; // total number of entries
};

#endif // __SYSTEMMAT2D__
