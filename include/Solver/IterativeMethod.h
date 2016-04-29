#ifndef __ITERATIONMETHOD__
#define __ITERATIONMETHOD__

#include <utility> // std::pair
#include <memory> // std::shared_ptr
#include <Preconditioner.h>

/** @brief an abstract base class to describe iterative methods for solving
 */
template <class LinearOperator, class Vector>
class IterativeMethod
{
  protected:
    /** @brief preconditioner */
    std::shared_ptr<Preconditioner<Vector>> prec;
    
    /** @brief just a name for nicer output */
    std::string name;

    /** @brief absolute tolerance for stopping 
     * 
     * This is the desired norm of the residual when calling 
     * IterativeMethod::iterate. After IterativeMethod::iterate it is the 
     * achieved norm of the residual. That's why there is a getter and a setter
     * method for this member variable.
     */
    double residual_tolerance;

    /** @brief reduction factor of residual for stopping
     * 
     * Stop the iteration if the current residual is smaller than the initial
     * residual multiplied with this factor.
     */
    double residual_reduction;

    /** @brief limit divergence for stopping
     *
     * Stop the iteration if the current residual is larger than the initial
     * residual multiplied with this factor.
     */
    double divergence_factor;

    /** @brief maximal number of iterations */
    unsigned int max_n_iterations;

    /** @brief minimal number of iterations */
    unsigned int min_n_iterations;

    /** @brief number of iterations until restart (for gmres) */
    unsigned int restart;
    
    /** @brief damping factor (typically between 0 and 1) */
    double damping;
    
    // a small helper function to reduce code duplication and to assure the same
    // stopping criterion for all iterative methods
    bool converged(double current_residual, double initial_residual,
                   unsigned int current_iterate)
    {
      bool residual_small_enough = current_residual < this->residual_tolerance;
      bool sufficient_reduction = 
        current_residual < this->residual_reduction * initial_residual;
      bool enough_iterations = 
        current_iterate > this->min_n_iterations || this->min_n_iterations == 0;
      return (residual_small_enough || sufficient_reduction)
              && enough_iterations;
    }
    
  public:
    /** @brief constructor */
    IterativeMethod(std::shared_ptr<Preconditioner<Vector>> prec,
                    std::string name = "")
      : prec(prec), name(name), residual_tolerance(1.e-8),
        residual_reduction(0.), divergence_factor(1.5), max_n_iterations(100),
        min_n_iterations(0), restart(10)
    {
    }
    
    /** destructor */
    virtual ~IterativeMethod() = default;

    /** @brief iterate 
     * 
     * @return the number of performed iterations and the residual as a pair
     */
    virtual std::pair<unsigned int, double> iterate(const LinearOperator & A,
                                                    const Vector & rhs,
                                                    Vector & solution) = 0;

    /** return absolute tolerance for stopping */
    double get_residual_tolerance() const
    { return residual_tolerance; };
    
    /// @brief set all parameters related to stopping criteria
    void set_stopping_parameters(unsigned int max_it, unsigned int min_it,
                                 double tolerance, double reduction,
                                 double divergence, double new_damping,
                                 unsigned int rest = 10)
    {
      residual_tolerance = tolerance;
      residual_reduction = reduction;
      divergence_factor = divergence;
      max_n_iterations = max_it;
      min_n_iterations = min_it;
      if(min_n_iterations > max_n_iterations)
        min_n_iterations = max_n_iterations; // do exactly this many iterations
      restart = rest;
      if(restart > max_n_iterations)
        restart = max_n_iterations;
      damping = new_damping;
    }
    
    /** set absolute tolerance for stopping */
    void set_residual_tolerance(double new_residual_tolerance)
    { residual_tolerance = new_residual_tolerance; }

    /** set relative tolerance for stopping */
    void set_residual_reduction(double new_residual_reduction)
    { residual_reduction = new_residual_reduction; }

    /** set tolerance for divergence */
    void set_divergence_factor(double new_divergence_factor)
    { divergence_factor = new_divergence_factor; }

    /** set maximal number of iterations */
    void set_max_n_iterations(unsigned int max_it)
    {
      max_n_iterations = max_it;
      if(min_n_iterations > max_n_iterations)
        min_n_iterations = max_n_iterations;// do exactly max_it many iterations
      if(restart > max_n_iterations)
        restart = max_n_iterations;
    }

    /** set minimal number of iterations */
    void set_min_n_iterations(unsigned int min_it)
    { 
      min_n_iterations = min_it;
      if(min_n_iterations > max_n_iterations)
        max_n_iterations = min_n_iterations;// do exactly min_it many iterations
    }

    /** set restart */
    void set_restart(unsigned int new_restart)
    {
      restart = new_restart;
      if(restart > max_n_iterations)
        restart = max_n_iterations;
    }
    
    void set_damping(double new_damping)
    {
      damping = new_damping;
    }
    
    std::string get_name() const
    { 
      return name;
    }
};

#endif // __ITERATIONMETHOD__