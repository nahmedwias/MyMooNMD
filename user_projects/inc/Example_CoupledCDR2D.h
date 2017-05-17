#ifndef _EXAMPLE_COUPLEDCDR2D_
#define _EXAMPLE_COUPLEDCDR2D_

#include <memory>
#include <Example_CD2D.h>
//for the typedef CouplingFunction.
#include <Constants.h>


class Example_CoupledCDR2D : public Example2D
{


  public:

    //! @brief Constructor,chooses example according to the example code
    //! which is given in the input database.
    Example_CoupledCDR2D(const ParameterDatabase &);

    /**
     * @brief Get one of the underlying CD(R) examples without the coupling term.
     * @param[in] n The index of the equation whose example we request.
     * @return A constant reference to one particular decoupled example.
     */
    const Example_CD2D& getDecoupledExample(size_t n) const;

    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    Example_CoupledCDR2D(const Example_CoupledCDR2D&) = default;

    //! Default move constructor.
    Example_CoupledCDR2D(Example_CoupledCDR2D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    Example_CoupledCDR2D& operator=(const Example_CoupledCDR2D&) = default;

    //! Default move assignment operator
    Example_CoupledCDR2D& operator=(Example_CoupledCDR2D&&) = default;

    //! Default destructor.
    ~Example_CoupledCDR2D() = default;


    //Getter functions.

    //! Return function pointer to coefficient function of equation nr. equationIndex.
    CoeffFct2D* getCoeffFct(size_t equationIndex) const;

    //! Return function pointer to assembling function of equation nr. equationIndex.
    AssembleFctParam2D* getAssemblingFct(size_t equationIndex) const;

    //! Return function pointer to parameter function of equation nr. equationIndex.
    ParamFct* getParamFct(size_t equationIndex) const;

    /*! Get number of equations
     * @return The value of nEquations.*/
    size_t getNEquations() const{
      return nEquations_;
    }

    CoeffFct2D* get_coeffs() const override;




  protected:
    //! Number of contained equations.
    size_t nEquations_;

    //!@brief List of the coefficient function pointers used in the assembling of the uncoupled (convection-diffusion) parts.
    // NOTE: This replaces usage of problem_coefficients from Base class. Imho best of the bad alternatives
    std::vector<CoeffFct2D*> bilinCoeffs_;
    //!@brief List of the coefficient function pointers used in the assembling of the rhs part (linear-decoupled strategy).
    std::vector<AssembleFctParam2D*> rhsAssemblingFunctions_;
    //!@brief List of the parameter function pointers used in the assembling of the rhs part (linear-decoupled strategy).
    std::vector<ParamFct*> parameterFunctions_;

    /** @brief A vector which stores the individual examples without the coupling part.
     *  @note It is important to note that although each of these decoupled examples is
     *  equipped with an "exact solution", this is of little relevance - depending on
     *  the implementation of the used hard-coded example, this exact solution is either the
     *  exact solution of the equation in a coupled system or, if unknown, just zero or some other
     *  reference function.
     */
    std::vector<Example_CD2D> decoupledExamples_;


  private:
    //! @brief Generates decoupled CD2D examples and fills the data member decoupledExamples_ with them.
    void generateDecoupledExamples();


};


#endif
