/** ***************************************************************************
 *
 * @name   ConvectionDiffusion_AFC
 * @brief  Stores everything needed to solve a Convection Diffusion problem
 *         using Algebraic Flux Correction.
 *
 *         Store the exact interpolant, gamma vector (for BJK limiter) and
 *         diffusion matrix D. Also contains optimization tools such as Anderson
 *         acceleration and Dynamic Damping. This wraps up everything which 
 *         is necessary to solve an AFC problem in 2D and 3D.
 *
 * 
 * @author     Volker John, Ellen Schmeyer, Clemens Bartsch, Abhinav Jha
 * @date       2015/11/11
 *
 ******************************************************************************/

#ifndef __SYSTEM_CONVECTIONDIFFUSION_AFC_H__
#define __SYSTEM_CONVECTIONDIFFUSION_AFC_H__

#include "ConvectionDiffusion.h"
#include "anderson.h"

template<int d>
class LocalAssembling;

template<int d>
class ConvectionDiffusion_AFC : public ConvectionDiffusion<d> 
{
  protected:
    
    /** @brief sets algebraic_flux_correction to afc and space discretization
     *         to Galerkin*/
    void set_AFC_parameters();
    /** @brief stores interpolation of the exact solution. Useful in finding
     *         AFC errors */
    std::vector<double> exact_interpolant;
    /** @brief vector of weights for an AFC scheme */
    std::vector<double> afc_gamma;
    /** @brief entries of correction matrix for AFC schemes */
    std::vector<double> afc_matrix_D_entries;
    
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_CD = typename Template_names<d>::Example_CD;

    ConvectionDiffusion_AFC(const TDomain& domain,
                            const ParameterDatabase& param_db);  
    
    /** @brief assemble the system
     * 
     *  @param[in] iteration Required for Fixed_point_RHS scheme. It helps in 
     *                       assembling the matrix only once.
     *
     */
    void assemble(const int iteration);
    
    /** @brief solve the system 
     * 
     * @param[in] iteration Required for factorization of system matrix only 
     *                      once in case of Fixed_point_RHS and then store that 
     *                      for further iterations.
     * 
     * Returns whether or not the method is converged
     */
    bool solve(const int iteration);
    
    /** 
     * @brief measure AFC errors
     * 
     * The current AFC errors will be printed out. 
     * 
     * @param[in] i use to call output from base class
     */
    void output(int i = -1);
    
    ///Dynamic Damping from [JK08.CMAME]
    void dynamic_damping(const int iteration);
        
    ///Anderson acceleration from [WN11.SINUM]
    ///check the paper for notations.
    void anderson_acceleration_damping(int N_Unknowns, int iteration, 
     std::list<std::vector<double>> & solAnderson,
     std::list<std::vector<double>> & deltaAnderson);
    
    
    ~ConvectionDiffusion_AFC()=default; 
    
  private:
    /**
     * @brief: Apply an algebraic flux correction scheme to the assembled matrix.
     * Should be called within the assemble routine, after the actual assembling
     * has been performed with the INTERNAL_FULL_MATRIX_STRUCTURE switch on.
     *
     * @param[in] iteration required to copy the matrix A+D after first 
     *                      iteration.
     * @param[in] check_fpr Checks the iteration scheme
     * 
     */
    void do_algebraic_flux_correction(const int iteration, const int check_fpr);
    
         
    /**
     * @brief Parameters for Anderson acceleration
     * 
     * alphas_x_i: weights used in Anderson acceleration
     * solAnderson, deltaAnderson: storage of anderson previous updates
     * kAnderson: number of steps in anderson acceleration
     */
    BlockVector alphas_x_i;
    std::list<std::vector<double>>  solAnderson;
    std::list<std::vector<double>>  deltaAnderson;
    int kAnderson;
    
    /**
     * residual: norm of current residual r_k+1
     * residual_old: norm of previous residue r_k 
     * old_solution: Storing the old solution, will be used to compute new 
     *               iterate     
     * time_total: time taken to solve the problem
     * rejected_steps: rejected iteration steps in dynamic_damping
     */
     double residual, residual_old; 
     BlockVector old_solution;
     double time_total;
     int rejected_steps;

    /**
     * @brief Parameters for Fixed Point RHS (FPR)
     * 
     * matrix_copy: For AFC scheme while computing afc_matrix_D_entries we need
     *              the matrix A computed as all dofs were active dofs. For FPR 
     * the assembly is done only once and hence we need to store the A matrix 
     * after first computation.
     * 
     * is_not_afc_fixed_point_rhs: checks whether we have FPR and iteration>1
     * 
     * rhs_copy: Copy the RHS after the first iteration is done so that it can 
     *           be used for FPR
     */
    BlockFEMatrix matrix_copy;
    int is_not_afc_fixed_point_rhs;
    BlockVector rhs_copy;
};

#endif