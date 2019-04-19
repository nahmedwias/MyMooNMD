/** ************************************************************************ 
*
* @class Example_TnseTcd2D
* @brief store all functions needed to describe coupled (Navier-)Stokes and convection--diffusion example 
* 
* The standard
* constructor of this class will fill the vectors (in Example) with pointers
* to the functions needed to fully describe a particular example.
* 
* @author 
* @date 
* 
* @ruleof0
* 
 ************************************************************************  */
#ifndef _Example_TnseTcd2D_
#define _Example_TnseTcd2D_

#include <Example_TimeNSE2D.h>
#include <Example_TimeCD2D.h>
#include <functional>
#include <Example_NonStationary2D.h>

template <int d> class MultiPhase; //forward declaration

class Example_TnseTcd2D : public Example_NonStationary2D
{
public:
  /** @brief default constructor
   * This intializes a Time dependent (Navier-)Stokes and convection-diffusion example in 2D. 
   * It is chosen according to example_code.
   */
  explicit Example_TnseTcd2D(const ParameterDatabase& user_input_parameter_db);
  /** @brief initialize your own example
   * 
   * Create an example with all vectors already defined.
   */
  Example_TnseTcd2D(const std::vector<DoubleFunct2D*>& exact,
                   const std::vector<BoundCondFunct2D*>& bc,
                   const std::vector<BoundValueFunct2D*>& bd, 
                   const CoeffFct2D& coeffs,                    
                   bool timedependentrhs, bool timedependentcoeffs, 
                   std::vector <DoubleFunct2D*> init_cond);
//       :  Example2D(const std::vector<DoubleFunct2D*>& exact,
//               const std::vector<BoundCondFunct2D*>& bc,
//               const std::vector<BoundValueFunct2D*>& bd, 
//               const CoeffFct2D& coeffs)

  
//---------------------------We need to change this aswell------------------------    
      /// Apply the function stored as post processing routine.
  void do_post_processing(MultiPhase<1> & tcd2d) const;
  void do_post_processing(MultiPhase<2>& tnse2d, double& val) const;
  void do_post_processing(MultiPhase<2>& tnse2d) const;

  /// Return kinematic viscosity, if set.
  double get_nu() const;
  

    //Declaration of special member functions - rule of zero
  //! Default copy constructor. Performs deep copy.
  Example_TnseTcd2D(const Example_TnseTcd2D&) = default;
  
  //! Default move constructor.
  Example_TnseTcd2D(Example_TnseTcd2D&&) = default;
  
  //! Default copy assignment operator. Performs deep copy.
  Example_TnseTcd2D& operator=(const Example_TnseTcd2D&) = default;
  
  //! Default move assignment operator
  Example_TnseTcd2D& operator=(Example_TnseTcd2D&&) = default;
  
  //! Default destructor.
  ~Example_TnseTcd2D() = default;

  private:
  /// Function doing the post processing for a stationary example.
  /// TODO put Time_CD2D argument const as soon as FEFunctions can be copied properly!
  //std::function<void(TNSE_TCD &)> post_processing_stat;
  
//    std::function<void(TNSE_TCD&, double& val)> post_processing_stat;
   /// TODO Function doing the post processing for a time dependent example.
 std::function<void(MultiPhase<2>&)> post_processing_stat_old;
  /// TODO Function doing the post processing for a time dependent example.

};
#endif // _Example_TnseTcd2D_
  
  
  
  
  

  
  
  
