#ifndef __ITERATION_GMRES__
#define __ITERATION_GMRES__

#include <IterativeMethod.h>

// not nice: pollute global namespace
enum class gmres_type {left, right, flexible};

template <class LinearOperator, class Vector>
class Iteration_gmres : public IterativeMethod<LinearOperator, Vector>
{
  public:
    /** constructor */
    Iteration_gmres(std::shared_ptr<Preconditioner<Vector>> p,
                    gmres_type t = gmres_type::flexible);
    
    /** destructor */
    virtual ~Iteration_gmres() = default;
    
    /** iterate routine */
    std::pair<unsigned int, double> iterate(const LinearOperator & A, 
                                            const Vector & rhs,
                                            Vector & solution);
    
  protected:
    
    gmres_type type;
    
    std::pair<unsigned int, double> left_gmres(const LinearOperator & A, 
                                               const Vector & rhs,
                                               Vector & solution);
    std::pair<unsigned int, double> right_gmres(const LinearOperator & A, 
                                                const Vector & rhs,
                                                Vector & solution);
    std::pair<unsigned int, double> flexible_gmres(const LinearOperator & A, 
                                                   const Vector & rhs,
                                                   Vector & solution);
};


#endif // __ITERATION_GMRES__