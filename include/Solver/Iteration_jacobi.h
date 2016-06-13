#ifndef __ITERATION_JACOBI__
#define __ITERATION_JACOBI__

#include <IterativeMethod.h>
#include <vector>

template <class LinearOperator, class Vector>
class Iteration_jacobi : public IterativeMethod<LinearOperator, Vector>,
                         public Preconditioner<Vector>
{
  public:
    /** constructor, this uses NoPreconditioner<Vector> ! */
    Iteration_jacobi(const LinearOperator & mat);
    
    void apply(const Vector & z, Vector & r) const override;
    
    /** destructor */
    virtual ~Iteration_jacobi() = default;
    
    /** iterate routine */
    std::pair<unsigned int, double> iterate(const LinearOperator & A, 
                                            const Vector & rhs,
                                            Vector & solution);
    
    /** @brief update after the matrix has changed */
    void update(const LinearOperator& A);
    void update(){ this->update(this->linear_operator); };
    
    const LinearOperator& get_operator() const;

  protected:
    // the diagonal entries
    std::vector<double> diagonal_entries;
    // the linear operator (usually the matrix)
    const LinearOperator & linear_operator;
};

#endif // __ITERATION_JACOBI__
