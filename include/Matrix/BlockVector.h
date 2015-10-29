/** ************************************************************************ 
*
* @name       BlockVector
* @brief      Store a block-vector 
*
*             Store block-vector which consists of several subvectors 
*             (diffenrent length is possible). Some algebraic operations on the
*             subvectors are implemented
*
* @author     Ulrich Wilbrandt
* @date       01.11.2013
*
****************************************************************************/

#ifndef __BLOCKVECTOR__
#define __BLOCKVECTOR__

class BlockMatrix; // forward declaration

#include <numeric>
#include <BlockMatrix.h>

class BlockVector
{
  protected:
    /** @brief vector of (pointers to) blocks. blocks[0] is the entire vector */
    std::vector<double> entries;
    
    /** @brief the lengths of each block 
     * 
     * The size of this vector also determines the number of blocks.
     */
    std::vector<unsigned int> lengths;
    
    /** @brief the number of active degrees of freedom for each block 
     * 
     * The size of this vector is the same as that of 'lengths'.
     */
    std::vector<unsigned int> actives;
    
    /** @brief return the accumulated length of all blocks i with i < b
     * 
     * This tells you the start index of a given block in the entries vector.
     */
    unsigned int offset(unsigned int b) const;
  public:
    
    /** standard constructor */
    BlockVector();
    
    /** constructor for a BlockVector consisting of a single block of length
     * 'l' filled with zeros.
     */
    BlockVector(unsigned int l);
    
    /** constructor for a BlockVector consisting of a single block of length
     * 'l' filled with zeros.
     *
     * @note This constructor is for backwards compatibility, to avoid
     * explicit casts of integers to unsigned int. Use carefully.
     */
    BlockVector(int l);

    /** @brief constructor for a BlockVector suitable for a given BlockMatrix
     * 
     * If 'image' is set to true, the vector will be in the image space of the 
     * BlockMatrix 'mat'. Otherwise it will be in the pre-image space. It
     * is initialized with zeros.
     * 
     * The default behavior, 'image' set to false, means this BlockVector will 
     * have as many blocks as 'mat' has blocks in each row. The length of 
     * this BlockVector will be the number of (total) columns of the given 
     * BlockMatrix 'mat'. 
     * 
     * For a linear system \f$Ax=b\f$, the solution \f$x\f$ can be created using
     * \f$A\f$ and 'image' set to false. With 'image' set to true a BlockVector 
     * suitable for \f$b\f$ is created.
     * 
     * In contrast to the templated constructor, all entries are always active.
     */
    BlockVector(const BlockMatrix& mat, bool image = false);
    
    /** @brief constructor for a BlockVector suitable for a given BlockMatrix
     * 
     * If 'image' is set to true, the vector will be in the image space of the 
     * BlockMatrix 'mat'. Otherwise it will be in the pre-image space. It
     * is initialized with zeros.
     * 
     * The default behavior, 'image' set to false, means this BlockVector will 
     * have as many blocks as 'mat' has blocks in each row. The length of 
     * this BlockVector will be the number of (total) columns of the given 
     * BlockMatrix 'mat'. 
     * 
     * For a linear system \f$Ax=b\f$, the solution \f$x\f$ can be created using
     * \f$A\f$ and 'image' set to false. With 'image' set to true a BlockVector 
     * suitable for \f$b\f$ is created.
     * 
     * In contrast to the non-templated constructor, this sets the the active
     * entries appropriatly. The class 'BM' should be a derived class from 
     * BlockMatrix.
     */
    template<class BM>
    BlockVector(const BM& mat, bool image = false)
    {
      this->copy_structure<BM>(mat, image);
    }

    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    BlockVector(const BlockVector&) = default;

    //! Default move constructor.
    BlockVector(BlockVector&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    BlockVector& operator=(const BlockVector&) = default;

    //! Default move assignment operator
    BlockVector& operator=(BlockVector&&) = default;

    //! Default destructor.
    ~BlockVector() = default;
    
    /**
     * @brief Set all entries to zero
     *
     * Reset all entries of all subvectors to zero.
     */
    void reset();
    
    /**
     * @brief Set all active entries to zero
     *
     * Reset all active entries of all subvectors to zero leaving only Dirichlet
     * values unchanged.
     * 
     */
    void ResetActive();
    
    /**
     * @brief Set all non-active entries to zero
     *
     * Reset all non-active entries of all subvectors to zero leaving all  
     * non-Dirichlet values unchanged.
     */
    void ResetNonActive();
    
    /**
     * @brief Scale a subvector 
     *
     * The subvector with index i is scaled by a. If i is negative, the entire
     * BlockVector is scaled
     *
     * @param a factor by which the subvector is scaled
     * @param i index of target subvector
     *
     */
    void scale(const double a, const unsigned int i);
    
    /**
     * @brief Scale the entire vector
     */
    void scale(const double a);
    
    /**
     * @brief add scaled vector to this
     *
     * @param r BlockVector which is added to this
     * @param factor factor with which r is multiplied
     *
     */
    void add_scaled(const BlockVector& r, double factor);
    
    
    
    /** @brief copy the structure of another BlockVector, 
     * 
     * No values are copied. If there are old values, they are deleted! New 
     * values are set to 0.
     * 
     * @param r BlockVector from which the structure is copied to this one
     */
    void copy_structure(const BlockVector & r);
    
    /** @brief change the structure of a BlockVector so that it is suitable for
     *         a given BlockMatrix
     * 
     * This is exactly what the constructor 
     * BlockVector(const BlockMatrix& mat, bool image = false) does.
     * 
     * This deletes possibly existing data in this BlockVector. After this 
     * method is called, all values are set to zero. All entries are active.
     */
    void copy_structure(const BlockMatrix & mat, bool image);
    
    /** @brief change the structure of a BlockVector so that it is suitable for
     *         a given BlockMatrix
     * 
     * This is essentially the method 
     * copy_structure(const BlockMatrix& mat, bool image), however this method 
     * is called using a derived class of BlockMatrix, which is used to set the
     * active entries in this BlockVector.
     * 
     * This deletes possibly existing data in this BlockVector. After this 
     * method is called, all values are set to zero. This sets active entries
     * appropriatly.
     * 
     * The class BM is supposed to be a derived class of BlockMatrix and it must
     * implement a method
     *  const TFESpace * get_space_of_block(unsigned int block, bool test) const
     * which returns the test space (or ansatz space if 'test' is false) of the
     * b-th block.
     */
    template<class BM>
    void copy_structure(const BM & mat, bool image);
    
    /**
     * @brief add (scaled) values into a subvector, "this += a*x"
     *
     * The values given by 'x' are added to the subvector with index 'i'. If 
     * 'i' is negative, the entire BlockVector is added. The scalar 'a' is a 
     * scaling factor.
     *
     * @param[in] x values which are copied to the target subvector
     * @param i index of target subvector, -1 means entire vector
     * @param a factor by which x is multiplied
     *
     */
    void add(const double * x, const int i = -1, double a = 1.0);

    /**
     * @brief copy values into a subvector
     *
     * The values given by x are added to the subvector with index 'i'.
     * If 'i' is negative, the entire BlockVector is added (as for operator=)
     *
     * @param x values which are copied to the target subvector
     * @param i index of target subvector
     *
     */
    void copy(const double * x, const int i = -1);
    
    /** @brief copy the nonactive values of each block of r to this */
    void copy_nonactive(const BlockVector& r);
    
    /**
     * @brief compute norm of this BlockVector
     * 
     * compute the 2-norm (square root of sum of squares)
     * Note: possibly implement other norms as well
     */
    double norm() const;
    
    /**
     * @brief Print subvector iB to console in Matlab format
     *
     * If iB < 0, the entire block vector is printed to console. Otherwise the 
     * subvector of index iB is printed.
     * Format: name(i)= val       (i>0)
     *
     * @param name name for output
     * @param iB index of subvector
     *
     */
    void print(const std::string name = "rhs", const int iB = -1) const;
    
    /**
     * @brief Print some information without explicitly printing values
     */
    void info();
    
    /** 
     * @brief write the entire BlockVector to a file
     * 
     * The file will start with a line containing only one number, indicating
     * the length of the vector. The following lines are just written and are
     * then not human readable. 
     * 
     * A file created with this method can be read into a BlockVector using 
     * BlockVector::read_from_file(std::string);
     */
    void write_fo_file(std::string filename);
    
    /** 
     * @brief read data into this BlockVector from a file
     * 
     * The file has to start with a line containing only one number, indicating
     * the length of the vector. The following lines are just read in and are
     * not human readable. 
     * 
     * A file created with BlockVector::write_fo_file(std::string) can be read 
     * using this method.
     */
    void read_from_file(std::string filename);
    
    /** getters */
    
    unsigned int length() const 
    { return std::accumulate(lengths.begin(),lengths.end(),0); }
    
    unsigned int length(const int i) const
    { return lengths.at(i); }
    
    /** @brief return the number of active entries for a given block i */
    unsigned int active(const int i) const
    { return actives.at(i); }
    
    unsigned int n_blocks() const 
    { return lengths.size(); }
    
    const double* block(const unsigned int i) const;
    
    double* block(const unsigned int i);
    
    const double* get_entries() const
    { return &(this->entries.at(0)); }
    
    double* get_entries()
    { return &(this->entries.at(0)); }
    
    /* *********************************************************************
     * overloading standard operators
     * we avoid operators +,-,*,/ because they would require to alocate new
     * memory, which is often undesirable and unnecessary.
     */
    
    /** array acces:
     * For a BlockVector v and unsigned int i use either 'v.at(i)', 'v[i]' or 
     * 'v(i)' to acces the i-th element. All three methods provide the same 
     * functionality. The function 'at()' also checks if the given index is not 
     * too large. Furthermore for each method two versions are implemented. One
     * returns an lvalue so that a value can be set at the given index 
     * (e.g. a[i]=5;). The other leaves this BlockVector constant 
     * (e.g. in cout << a[i]);
     */
    double & at(const unsigned int i);
    const double & at(const unsigned int i) const;
    double& operator[](const unsigned int i)
    { return entries[i]; }
    const double& operator[](const unsigned int i) const
    { return entries[i]; }
    double& operator()(const unsigned int i)
    { return entries[i]; }
    const double& operator()(const unsigned int i) const
    { return entries[i]; }
    
    
    // copy, r should be as long as entire BlockVector
    BlockVector& operator=(const double *r);
    BlockVector& operator=(const double a); // set all values to a
    BlockVector& operator*=(const double a); // multiply by a scalar
    // mulitply matrix from left (i.e. x <-- A*x), only square A supported
    BlockVector& operator*=(const BlockMatrix& A);
    BlockVector& operator+=(const BlockVector& r); // add other BlockVector
    BlockVector& operator-=(const BlockVector& r);// substract other BlockVector
    
    /** ******************************************************************** */
    
    /** friends */
    
    /** 
     * @brief compute dot product of two BlockVectors
     *
     * it is checked if the two BlockVectors are of equal length
     *
     * @note we put 'friend' here, so that we don't need a static function call
     *       in e.g. a template iterative solver
     * @param a,b the two BlockVectors
     */
    friend double dot(const BlockVector& a, const BlockVector& b);
    friend double norm(const BlockVector& a) { return a.norm(); };
};



#endif

