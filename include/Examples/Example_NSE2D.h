/** ************************************************************************ 
*
* @class Example_NSE2D
* @brief store all functions needed to describe a (Navier-)Stokes example 
* 
* Depending on the value of TDatabase::ParamDB->EXAMPLE, the standard 
* constructor of this class will fill the vectors (in Example2D) with pointers
* to the functions needed to fully describe a particular example.
* 
* @author    Ulrich Wilbrandt, 
* @date      14.03.15
* 
* @ruleof0
* 
 ************************************************************************  */

#ifndef __EXAMPLE_NSE2D__
#define __EXAMPLE_NSE2D__

#include<Example2D.h>
#include <functional>
#include <FEFunction2D.h>

class Example_NSE2D : public Example2D
{
  /** @brief a function which can be called as a post processing step 
   *
   * The arguments are the two velocity components and the pressure function. 
   * This is for stationary problems
   */
  std::function<void(TFEFunction2D *, TFEFunction2D *,TFEFunction2D *)>
    post_processing;

  /** @brief A function which can be called as a post processing step
   *    in the time dependent case
   *
   * The arguments are the two velocity components and the pressure function
   * at the last time step, plus the two velocity components at the
   * last-but-one time step.
   */
  std::function<void(TFEFunction2D *, TFEFunction2D *,TFEFunction2D *,TFEFunction2D *,TFEFunction2D *)>
    post_processing_time;
    
  public:
    /** @brief default constructor
     * 
     * This intializes a (Navier-)Stokes example in 2D. It is chosen according
     * to TDatabase::ParamDB->EXAMPLE.
     */
    Example_NSE2D();
    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example_NSE2D(std::vector <DoubleFunct2D*> exact,
                  std::vector <BoundCondFunct2D*> bc,
                  std::vector <BoundValueFunct2D*> bd, CoeffFct2D *coeffs)
      : Example2D(exact, bc, bd, coeffs) {};
  
    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    Example_NSE2D(const Example_NSE2D&) = default;

    //! Default move constructor.
    Example_NSE2D(Example_NSE2D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    Example_NSE2D& operator=(const Example_NSE2D&) = default;

    //! Default move assignment operator
    Example_NSE2D& operator=(Example_NSE2D&&) = default;

    //! Default destructor.
    ~Example_NSE2D() = default;

    void do_post_processing(TFEFunction2D *u1, TFEFunction2D *u2,TFEFunction2D *p) const
    { 
      if(this->post_processing)
        this->post_processing(u1, u2, p);
    }
  /*! @brief Example dependent post processing and output generating method.
   *
   * Checks if the std::function member post_processing_time was initialized
   * for the current example and if so calls the method to do whatever it was
   * intended to do.
   */
  void do_post_processing(TFEFunction2D *u1, TFEFunction2D *u2, TFEFunction2D *p,
                          TFEFunction2D *u1old, TFEFunction2D *u2old) const
  {
    if(this->post_processing_time)
      this->post_processing_time(u1, u2, p, u1old, u2old);
  }
};


#endif // __EXAMPLE_NSE2D__
