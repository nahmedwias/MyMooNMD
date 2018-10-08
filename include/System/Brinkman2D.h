/** ***************************************************************************
 *
 * @name   Brinkman2D
 * @brief  store everything needed to solve a Brinkman problem
 *
 *         Store the matrix, right hand side, FE spaces, FE functions and the 
 *         solution vector of a Brinkman problem. This wraps up everything which 
 *         is necessary to solve a Brinkman problem in 2D.
 *
 * @author     Alfonso Caiazzo and Laura Blank
 * @date       18.07.2016
 *
 ******************************************************************************/

#ifndef Brinkman2D_H_
#define Brinkman2D_H_

#include <FEVectFunct2D.h>
#include <Example_Brinkman2D.h>
#include <MainUtilities.h> // FixedSizeQueue
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Solver.h>
#include <ParameterDatabase.h>
#include <DataWriter.h>
#include <utility>
#include <array>
#include <FEFunction2D.h>



class Brinkman2D
{
// This is only interesting if we consider different matrix types, see "todo" in Brinkman2D.C
   enum class Matrix{Type14, Type1, Type2, Type3, Type4};

public:
    /**
     * @brief a simple struct storing one set of residuals
     *
     * The full residual is the \f$\ell^2\f$-norm of the vector \f$Ax-b\f$
     * where \f$A\f$ is the current matrix and \f$b\f$ the right hand side. It
     * is composed of two parts, the impuls and the residual.
     *
     * If not default constructed it holds
     *     fullResidual*fullResidual = impulsResidual*impulsResidual
     *                                 +massResidual*massResidual
     */
    struct Residuals
    {
        /// @brief the impuls residual
        double impulsResidual;
        /// @brief the mass residual
        double massResidual;
        /// @brief the fulf residual
        double fullResidual;
        ///@brief standard constructor, initialize with large numbers
        Residuals();
        /// @brief constructor given the \e square of the impuls and mass
        /// residuals
        Residuals(double imR, double maR);
        /// @brief write out the three numbers to a stream.
        friend std::ostream& operator<<(std::ostream& s, const Residuals& n);
    };
    
//    BoundCondFunct2D *boundary_conditions[3];
   //BoundValueFunct2D **BoundaryValues;

    
 // Initialize/Declare Brinkman 2D database, called with the constructor (for the use outside Brinkman2D.C)
 static ParameterDatabase get_default_Brinkman2D_parameters();


    /** @brief constructor
     *
     * This constructor calls the other constructor creating an Example_Brinkman2D
     * object for you. See there for more documentation.
     */
    Brinkman2D(const TDomain& domain, const ParameterDatabase& param_db,
               int reference_id = -4711);
    
    /** @brief constructor
     *
     * The domain must have been refined a couple of times already if you want
     * to use multigrid. On the finest level the finite element spaces and
     * functions as well as matrices, solution and right hand side vectors are
     * initialized.
     *
     * The reference_id can be used if only the cells with the give reference_id
     * should be used. The default implies all cells.
     */
    Brinkman2D(const TDomain & domain, const ParameterDatabase& param_db,
               const Example_Brinkman2D & _example,
               unsigned int reference_id = -4711);
    
    /** @brief standard destructor */
    ~Brinkman2D();
    
    /** @brief assemble matrix,
     *
     * This assembles everything which is not related to the nonlinear term.
     * I.e. it assembles a Stokes matrix.
     */

// LB NEW 16.04.18 start
void assemble(size_t level = 4711, TFEFunction2D* coefficient_function = nullptr);
TFEFunction2D* u1;
TFEFunction2D* u2;
// LB NEW 16.04.18 end
/*
// LB OLD 16.04.18 start
    void assemble();
// LB OLD 16.04.18 end
*/
    /** @brief solve the system */
    void solve();
    
    /**
     * @brief Compute the defect Ax-b and store its norm in Brinkman2D::oldResiduals
     *
     * where A is the current matrix, x is the current solution and b is the
     * right hand side. Call this function after assembling the nonlinear
     * matrix with the current solution.
     */
    void computeNormsOfResiduals();
    
    /** @brief check if one of the stopping criteria is fulfilled
     *
     * either converged, maximun number of iterations reached, or slow
     * convergence
     *
     * @param iteration_counter current iterate
     */
    bool stopIt(unsigned int iteration_counter);
    
    /**
     * @brief measure errors and write pictures
     *
     * The current errors will be printed out. If desired, further output, e.g.,
     * vtk files are created.
     *
     * @param i suffix for output file name, -1 means no suffix
     */
    void output(int level, int i = -1);
    
    /// @name return computed errors
    ///
    /// You have to call Brinkman2D::output for any of these to return a
    /// meaningful value.
    //@{
    /// @brief return the computed L2 error of the velocity
    double getL2VelocityError() const;
    /// @brief return the computed L2 error of the divergence of the velocity
    double getL2DivergenceError() const;
    /// @brief return the computed H1-semi error of the velocity
    double getH1SemiVelocityError() const;
    /// @brief return the computed L2-error of the velocity at the boundary
    double getL2BoundaryError() const;
    /// @brief return the computed L2-error of the normal velocity at the boundary
    double getL2NormNormalComponentError() const;
    /// @brief return the computed L2 error of the pressure
    double getL2PressureError() const;
    /// @brief return the computed L2 error of the pressure
    double getH1SemiPressureError() const;
    //@}

    /////////////// Routines for periodic boundary conditions /////////////////

    // for the riverbed example with "nonhomogeneous periodic boundary conditions".
    // this map contains all pairs of velocity dofs which have to be identified
    std::map<int,int> periodic_dofs;


    /** find periodic boundaries dofs.·
     * This fills the map<int,int> 'periodic_dofs' such that a call to·
     * 'getPeriodicDOF(int)' now makes sense
     */
    void findPeriodicDOFs();
    /** enable periodic boundaries in some matrix which has as many rows as the
     * Stokes matrix. This is used for·
     * - S, the Stokes matrix itself  (default, i.e., no arguments)
     * - C, coupling matrix used for the direct solution
     * - E, representing (eta_f,u.n), which is added to the rhs during iteration
     *·
     * The method makePeriodicBoundary() (without arguments must be called first
     * If the second argument is true then there will be ones on the diagonal and·
     * -1 in the off-diagonal entry coupling with this row. If additionally the·
     * third argument is true, then there are zeros put where otherwise 1 and -1
     * are put (as described in the previous sentence). This enlarges the·
     * structure of the matrix 'mat', it is then possible to add two such·
     * matrices. This is currently implemented only due to periodic boundary
     * values (riverbed example).··
     */
		 void makePeriodicBoundary(std::shared_ptr<TMatrix> mat = nullptr,
				 bool stokesMat = false, bool p = false);

     void checkPeriodicDOFs();

     /** check if the given dof is a periodic dof. If no, -1 is returned. If yes,·
      * the dof to which this dof is coupled (periodicDOF) is returned. The row·
      * 'dof' should then be deleted. If mat==getMat().squareBlock(0), the row
      * 'dof' should be replaced by a 1 on the diagonal and -1 on 'periodicDOF'.
      */
			int getPeriodicDOF(int dof) const
      {
      	std::map<int,int>::const_iterator it = periodic_dofs.find(dof);
      	if (it == periodic_dofs.end()) return -1;
      	else return it->second;
      }

			/////////////// /////////////// /////////////// ///////////////

    // getters and setters
    //    const BlockMatrixNSE2D & get_matrix() const TODO
    //    { return this->systems.front().matrix; }
    //    BlockMatrixNSE2D & get_matrix()
    //    { return this->systems.front().matrix; }
    const BlockFEMatrix & get_matrix() const
    { return this->systems.front().matrix; }
    BlockFEMatrix & get_matrix()
    { return this->systems.front().matrix; }
    const BlockVector & get_rhs() const
    { return this->systems.front().rhs; }
    BlockVector & get_rhs()
    { return this->systems.front().rhs; }
    const TFEVectFunct2D & get_velocity() const
    { return this->systems.front().u; }
    TFEVectFunct2D & get_velocity()
    { return this->systems.front().u; }
    // try not to use this as it is not const
    TFEFunction2D *get_velocity_component(int i)
    { return (i==0) ? this->systems.front().u.GetComponent(0)
        : this->systems.front().u.GetComponent(1); }
    const TFEFunction2D & get_pressure() const
    { return this->systems.front().p; }
    TFEFunction2D & get_pressure()
    { return this->systems.front().p; }
    const TFESpace2D & get_velocity_space() const
    { return this->systems.front().velocity_space; }
    const TFESpace2D & get_pressure_space() const
    { return this->systems.front().pressure_space; }
    const BlockVector & get_solution() const
    { return this->systems.front().solution; }
    BlockVector & get_solution()
    { return this->systems.front().solution; }
    unsigned int get_size() const
    { return this->systems.front().solution.length(); }
    const Example_Brinkman2D & get_example() const
    { return example; }
    const ParameterDatabase & get_db() const
    { return brinkman2d_db; }
    /// @brief get the current residuals  (updated in Brinkman2D::normOfResidual)
    const Residuals& getResiduals() const;
    /// @brief get the current impuls residual (updated in Brinkman2D::normOfResidual)
    double getImpulsResidual() const;
    /// @brief get the current mass residual (updated in Brinkman2D::normOfResidual)
    double getMassResidual() const;
    /// @brief get the current residual (updated in Brinkman2D::normOfResidual)
    double getFullResidual() const;
    /// @brief return the computed errors
    std::array<double, int(8)> get_errors();


protected:
    
    /** @brief store a complete system on a particular grid
     *
     * This combines a matrix, rhs, solution, spaces and functions needed to
     * describe one Brinkman problem in 2D.
     */
    struct System_per_grid
    {
        /** @brief Finite Element space for the velocity */
        TFESpace2D velocity_space;
        /** @brief Finite Element space for the pressure */
        TFESpace2D pressure_space;
        /** The system matrix. */
        BlockFEMatrix matrix;
        /** @brief the right hand side vector */
        BlockVector rhs;
        /** @brief solution vector with two components. */
        BlockVector solution;
        /** @brief Finite Element function for velocity */
        TFEVectFunct2D u;
        /** @brief Finite Element function for pressure */
        TFEFunction2D p;
        
        /** @brief constructor */
        System_per_grid(const Example_Brinkman2D& example, TCollection& coll,
                        std:: pair <int,int> velocity_pressure_orders,
                        Brinkman2D::Matrix type);
        
        /**
         * Special member functions mostly deleted,
         * for struct takes ownership of the bad
         * classes TFEFunction2D, TFEVectFunct2D and TFESpace2D.
         */
        //! Delete copy constructor.
        System_per_grid(const System_per_grid&) = delete;
        
        //! Delete move constructor.
        System_per_grid(System_per_grid&&) = delete;
        
        //! Delete copy assignment operator.
        System_per_grid& operator=(const System_per_grid&) = delete;
        
        //! Delete move assignment operator.
        System_per_grid& operator=(System_per_grid&&) = delete;
        
        //! Default destructor.
        ~System_per_grid() = default;
    };
    
    /** @brief a local parameter database which constrols this class
     *
     * The database given to the constructor will be merged into this one. Only
     * parameters which are of interest to this class are stored (and the
     * default ParMooN parameters). Note that this usually does not include
     * other parameters such as solver parameters. Those are only in the
     * Solver object.
     */

     // LB NEW 19.04.18 start
    ParameterDatabase brinkman2d_db;
    // LB NEW 19.04.18 end
/*      // LB OLD 19.04.18 start
    ParameterDatabase brinkman2d_db;
    // LB OLD 19.04.18 end
*/
    
    /** @brief a solver object which will solve the linear system
     *
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
     */
    Solver<BlockFEMatrix, BlockVector> solver;

    /** @brief class for output handling (vtk and case files) */
    DataWriter2D outputWriter;
    
    /** @brief a complete system on each grid
     *
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<System_per_grid> systems;
    
    /** @brief Definition of the used example */
    const Example_Brinkman2D example;
    
    /** @brief an array to store defect, so that we don't have to reallocate
     *         so often
     */
    BlockVector defect;
    
    /** @brief stores the norms of the residuals of previous iterations.
     * The default length is 10
     */
    std::vector<double> norms_of_residuals;

    /**
     * @brief store the norms of residuals from previous iterations
     */
    FixedSizeQueue<10, Residuals> oldResiduals;

    /** @brief store the initial residual so that the nonlinear iteration can
     *         be stopped as soon as a desired reduction is achieved
     */
    double initial_residual;
    
    /** @brief set the velocity and pressure orders
     *
     * This function sets the corresponding velocity and
     * pressure orders. The pressure order is set if it is
     * not specified by the readin file. Default is -4711
     */
    void get_velocity_pressure_orders(std::pair <int,int>
                                      &velocity_pressure_orders);
    
    /** @brief set parameters in database
     *
     * This functions checks if the parameters in the database are meaningful
     * and resets them otherwise. The hope is that after calling this function
     * this class is fully functional.
     *
     * If some parameters are set to unsupported values, an error occurs and
     * throws an exception.
     */
    void set_parameters();

    // This function checks if the input Parameters are set in an 'admissible' way.↲
    // (Basicalley it chhecks the 'consistency' of symmetry and nonsymmetry)
    void check_input_parameters();

    /** @brief Errors to be accesed from outside the class
     * The array is filled during the function call Brinkman2D::output()
     * Currently, the errors store the L2 and H1 errors of the velocity
     * and pressure
     */
    std::array<double, int(8)> errors;
   



};


#endif /* Brinkman2D_H_ */
