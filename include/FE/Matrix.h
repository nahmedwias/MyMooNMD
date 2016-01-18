/** ************************************************************************ 
*
* @class     TMatrix
* @brief     store a sparse matrix
* 
* Sparse matrices in ParMooN are stored in <em>compressed row storage</em>. 
* Basically the matrix consists of a structure and entries and can perform 
* algebraic operations.
* 
* To create a TMatrix either copy or move an existing one or use a TStructure 
* to construct a new one.
* 
* Objects ob type TMatrix are not supposed to be copied. If you want that 
* copy/move construct a new TMatrix.
* 
* \todo apply naming convention to all methods
* 
* @ruleof0
* 
 ************************************************************************  */

#ifndef __MATRIX__
#define __MATRIX__

#include <Structure.h>
#include <string>
#include <map>
#include <vector>

class TMatrix
{
  protected:
    /** @brief structure of the matrix
     * 
     * This stores the information where this matrix has entries. Many objects
     * of type TMatrix can have the same TStructure, therefore we use a 
     * `std::shared_ptr` here.
     */
    std::shared_ptr<TStructure> structure;
    
    /** @brief matrix entries
     * 
     * Its size is determined by the TMatrix::structure.
     */
    std::vector<double> entries;
    
    /**
     * @brief replace the structure by a copy
     * 
     * If you plan on calling a method which changes the structure, such as
     * TMatrix::remove_zeros, TMatrix::changeRows or TMatrix::reorderMatrix,
     * this method has to be called before. Otherwise other matrices sharing 
     * the same TStructure will become corrupted. 
     * 
     * This method makes a deep copy of the structure. Then this matrix only 
     * knows its copy and is the only one knowing this copy. Then this matrix 
     * can safely call any structure changing methods.
     * 
     * Note that this method does nothing, if this matrix is the only one 
     * sharing its structure.
     */
    void copyOwnStructure();
    
  public:
    /** @brief generate the matrix, initialize entries with zeros */
    TMatrix(std::shared_ptr<TStructure> structure);

    /**
     * @brief Generates an empty `nRows`*`nCols` Matrix with no entries
     */
    TMatrix(int nRows, int nCols);

    /// @brief Default copy constructor
    ///
    /// Performs a deep copy of the entries, but not of the structure.
    TMatrix(const TMatrix&) = default;

    /// @brief Default move constructor.
    TMatrix(TMatrix&&) = default;

    /// @brief no copy assignment operator to avoid accidental copies
    TMatrix & operator=(const TMatrix& A) = delete;
    
    /// @brief no move assignment operator to avoid accidental moves
    TMatrix& operator=(TMatrix&&) = delete;

    /// @brief Default destructor.
    virtual ~TMatrix() = default;
    
    
    /// @brief reset all matrix entries to zero
    void reset();
    
    /// @brief return number of rows
    int GetN_Rows() const
    { return structure->GetN_Rows(); }
    
    /// @brief return number of columns
    int GetN_Columns() const
    { return structure->GetN_Columns(); }
    
    /// @brief return number of matrix entries
    int GetN_Entries() const
    { return structure->GetN_Entries(); }
    
    /// @brief return the column pointer in the TStructure of this matrix
    const int *GetKCol() const
    { return structure->GetKCol(); }
    /** @brief return the column pointer in the TStructure of this matrix
     * 
     * This version should never be used. It only exists because some other
     * parts of the software do not respect the const keyword (like AMG) or are
     * not well implemented (changing the structure of a matrix).
     */
    int *GetKCol()
    { return structure->GetKCol(); }
    
    /// @brief return the row pointer in the TStructure of this matrix
    const int *GetRowPtr() const
    { return structure->GetRowPtr(); }
    
    /** @brief return the row pointer in the TStructure of this matrix
     * 
     * This version should never be used. It only exists because some other
     * parts of the software do not respect the const keyword (like AMG) or are
     * not well implemented (changing the structure of a matrix).
     */
    int *GetRowPtr()
    { return structure->GetRowPtr(); }
    
    /// @brief return number of matrix entries for hanging node data
    int GetHangingN_Entries() const
    { return structure->GetHangingN_Entries(); }
    
    /// @brief return the column pointer corresponding to hanging nodes
    const int *GetHangingKCol() const
    { return structure->GetHangingKCol(); }
    
    /// @brief return the row pointer corresponding to hanging nodes
    const int *GetHangingRowPtr() const
    { return structure->GetHangingRowPtr(); }
    
    /// @brief return structure
    const TStructure& GetStructure() const
    { return *structure; }
    
    /// @brief return matrix entries as a pointer to const double
    const double *GetEntries() const
    { return &entries[0]; }
    
    /// @brief return matrix entries as a pointer to double
    double *GetEntries()
    { return &entries[0]; }
    
    /** @brief return the norm of the matrix 
     * 
     * The parameter \p p determines which norm to compute. Choose \p as  
     * -2 for Frobenius norm
     * -1 for maximum absolute row sum
     *  0 for maximum entry
     *  1 for maximum absolute column sum (not yet implemented)
     *  2 for euclidean norm, (not yet implemented)
     */
    double GetNorm(int p=-1) const;

    /// @brief write matrix into file
    int Write(const char *filename) const;
    
    /// @brief Print matrix into the shell
    void Print(const char *name = "a") const;
    
    /** @brief print the full matrix, including all zeros
     * 
     * This is only meaningful for very small matrices.
     */
    void PrintFull(std::string name="", int fieldWidth=4) const;

    /// @brief add a value at selected entry
    void add(int i,int j, double val);
    /** @brief add values in row 'i' given by the map 'vals', multiplied by 
     * 'factor'
     * 
     * This should be faster than adding all values in 'vals' individually
     */
    void add(int i, std::map<int,double> vals, double factor = 1.0);
    /** @brief add values `vals[i][j]` to this matrix at the positions `(i,j)`
     * for all `i,j` defined in the map `vals`
     */
    void add(std::map<int, std::map<int,double> > vals, double factor = 1.0);
    /// @brief set a value at selected entry
    void set(int i, int j, double val);
    /// @brief get a value at selected entry
    const double& get(int i,int j) const;
    /// @brief  get a value at selected entry (you may change that value)
    double& get(int i,int j);
    
    /// @brief reset the entries, the given vector must be of the same size
    void setEntries(std::vector<double> entries);
    
    /** @brief reorders the Matrix to comply with direct solvers. 
     * 
     * @warning This changes the structure of the matrix
     */
    void reorderMatrix();
    
    /// @brief return ordering of columns, see TStructure::ColOrder
    int GetColOrder() const
    { return structure->GetColOrder(); }
    
    /** @brief return a new TMatrix which is the transposed of this matrix 
     * 
     * If this is an object of a derived class (e.g. TMatrix2D, TSquareMatrix),
     * then the number of active degrees of freedom is not taken into account. 
     * The returned TMatrix is really the algebraic transposed matrix.
     * */
    TMatrix* GetTransposed() const;
    
    /** 
     * @brief replace several rows in the matrix with new entries.
     * 
     * Replace rows by new ones. This creates a new structure for the sparsity 
     * pattern of the matrix. Therefore reallocation is necessary. 
     * 
     * If there are no rows to change, i.e. if entries.size()==0, nothing is 
     * done.
     * 
     * This will create a new structure for this matrix. The structure 
     * previously belonging to this matrix is not changed. So other matrices 
     * are not affected.
     * 
     * @param entries for every row a map of columns-to-entries map
     */
    void changeRows(std::map<int,std::map<int,double> > entries);
    
    
    /** @brief compute y = A*x   (Matrix-Vector-Multiplication)
     *
     * Note that 'y' is created here and it is up to the user to delete it. 
     */
    friend double* operator*(const TMatrix & A, const double* x);

    /** @brief add another matrix to this one
     * 
     * This is of course only possible if the corresponding structures are the
     * same. 
     */
    virtual TMatrix & operator+=(const TMatrix * A);
    /** @brief substract another matrix to this one
     * 
     * This is of course only possible if the corresponding structures are the
     * same. 
     */
    virtual TMatrix & operator-=(const TMatrix * A);
    /** @brief add another matrix to this one
     * 
     * This is of course only possible if the corresponding structures are the
     * same. This method exists only for convenience and uses the same method 
     * with a pointer to A instead of the reference.
     */
    virtual TMatrix & operator+=(const TMatrix & A) 
    { *this += &A; return *this; }
    /** @brief scale matrix by a factor */
    virtual TMatrix & operator*=(const double a);
    
    /** @brief compute y += a * A*x
     *
     * 'A' is this TMatrix and 'x', and 'y' are given vectors. The scalar 'a'
     * is a scaling factor.
     * 
     * @param x array representing the vector which is multiplied by this matrix
     * @param y array representing the vector to which a * A*x is added
     * @param a scaling factor, default is 1.0
     */
    void multiply(const double * const x, double * const y, double a = 1.0) const;
    
    /** @brief compute y += a * A^T * x
     *
     * 'A^T' is the transposed of this TMatrix and 'x', and 'y' are given 
     * vectors. The scalar 'a' is a scaling factor. 
     * 
     * @param x array representing the vector which is multiplied by this matrix
     * @param y array representing the vector to which a * A*x is added
     * @param a scaling factor, default is 1.0
     */
    void transpose_multiply(const double * const x, double * const y, double a = 1.0)
      const;
    
    /**
     * @brief compute matrix-matrix product C = a*A*B, 
     * 
     * 'A' is this matrix, 'a' is a scalar factor, 'B' is given. Then matrix
     * 'C' is created during this function and the user is responsible to 
     * delete C.
     * 
     * Note that this is rather slow.
     * 
     * @param B matrix to be multiplied (from right) to this matrix
     * @param a scaling factor, default is 1.0
     */ 
    TMatrix* multiply(const TMatrix& B, double a = 1.0) const;
    
    /** @brief adding a scaled matrix to this matrix
     * 
     * The summation is index-wise, i.e. A(i,j) += factor*m(i.j), where A is 
     * this matrix. 
     * 
     * Note that this only works if the sparsity structure is the same for this
     * matrix and m.
     */
    void add_scaled(const TMatrix& m, double factor = 1.0);
    
    /// @brief scale all entries of this matrix
    void scale(double factor);
    
    /**
     * @brief scale a matrix using a vector
     * 
     * think of this as multipling this matrix with a diagonal matrix from the 
     * left (if the second argument is true). The parameter 'factor' are the 
     * entries of that diagonal matrix.
     * 
     * The i-th row (colum if from_left is false) is scaled by the i-th entry in
     * 'factor'.
     * 
     * The array 'factor' must be of size this->GetN_Rows() if 'from_left' is 
     * true or of size this->GetN_Columns() otherwise.
     * 
     * If all entries in 'factor' are the same, you can use operator*= as well.
     * 
     * @param factor array of scaling factors
     * @param from_left scale rows (true) or columns (false)
     */
    void scale(const double * const factor, bool from_left = true);
    
    /**
     * @brief remove all entries from sparsity structure where a zero is stored
     * 
     * This changes the sparsity structure of this matrix. Afterwards all stored
     * entries are nonzero. This can help if a lot of zeros are explicitly 
     * stored in the matrix. Afterwards matrix operations should be faster.
     * 
     * if tol is greater or equal to zero, all entries with magnitude smaller 
     * than tol are removed. If tol is smaller than zero, tol is set such that 
     * the ratio of the largest over the smallest entry in the resulting matrix
     * is smaller than 10^{15}.
     */
    void remove_zeros(double tol = 0.0);
    
    /** @brief get/set a specific matrix entry 
     * 
     * This will give an error if that entry is not in the sparsity structure
     */
    double & operator()(const int i, const int j);
    /** @brief get a specific matrix entry 
     * 
     * This will give an error if that entry is not in the sparsity structure
     */
    const double & operator()(const int i, const int j) const;
    
    /// @brief print some information on this TMatrix
    void info(size_t verbose) const;

    ///! @return true if this matrix is square
    bool is_square() const
    {
      return structure->isSquare();
    }
};

#endif
