/** ***************************************************************************
 *
 * @name   Brinkman3D
 * @brief  Store everything needed to solve a Brinkman problem
 *
 *         Store matrix, right hand side, FE spaces, FE functions and the
 *         solution vector of a Brinkman problem. This wraps up everything which
 *         is necessary to solve a Brinkman problem in 3D.
 *
 * @author     Alfonso Caiazzo and Laura Blank
 * @date       06.10.2016
 *
 ******************************************************************************/

#ifndef Brinkman3D_H_
#define Brinkman3D_H_

#include <FEVectFunct3D.h>
#include <Example_Brinkman3D.h>
#include <MainUtilities.h> // FixedSizeQueue
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Solver.h>
#include <ParameterDatabase.h>
#include <utility>
#include <array>
#include <FEFunction3D.h>
#include <Residuals.h>
#include <FESpace3D.h>
#include <DataWriter.h>
#include <vector>
#include <deque>
#include <list>


class Brinkman3D
{
    enum class Matrix{Type14, Type1, Type2, Type3, Type4};
    
protected:
    
    /** @brief store a complete system on a particular grid
     *
     * This combines a matrix, rhs, solution, spaces and functions needed to
     * describe one Darcy problem in 3D.
     */
    struct System_per_grid
    {
        /** @brief constructor
         * @param[in] example The current example.
         * @param[in] coll The collection of mesh cells of this grid.
         */
        System_per_grid(const Example_Brinkman3D& example, TCollection& coll,
                        std:: pair <int,int> velocity_pressure_orders,
                        Brinkman3D::Matrix type);
        /** @brief Finite Element space for the velocity */
        TFESpace3D velocity_space;
        /** @brief Finite Element space for the pressure */
        TFESpace3D pressure_space;
        /** The system matrix. */
        BlockFEMatrix matrix;
        /** @brief the right hand side vector */
        BlockVector rhs;
        /** @brief solution vector with two components. */
        BlockVector solution;
        /** @brief Finite Element function for velocity */
        TFEVectFunct3D u;
        /** @brief Finite Element function for pressure */
        TFEFunction3D p;
        
        /**
         * Special member functions mostly deleted,
         * for struct takes ownership of the bad
         * classes TFEFunction3D, TFEVectFunct3D and TFESpace3D.
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
    
    /** @brief a complete system on each grid
     *
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<System_per_grid> systems;
    
    /** @brief Definition of the used example */
    const Example_Brinkman3D example;
    
    /** @brief a local parameter database which constrols this class
     *
     * The database given to the constructor will be merged into this one. Only
     * parameters which are of interest to this class are stored (and the
     * defualt ParMooN parameters). Note that this usually does not include
     * other parameters such as solver parameters. Those are only in the
     * Solver object.
     */
    ParameterDatabase brinkman3d_db;
    
    /** @brief a solver object which will solve the linear system
     *
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
     */
    Solver<BlockFEMatrix, BlockVector> solver;
    
    /** @brief an array to store defect, so that we don't have to reallocate
     *         so often
     */
    BlockVector defect;
    
    //! @brief The initial residual. Stored so that the nonlinear iteration can
    //!        be stopped as soon as a desired reduction is achieved
    double initial_residual;
    
    ///@brief The norms of residuals from up to 10 previous iterations
    FixedSizeQueue<10, Residuals> norms_of_old_residuals;
    
    /** @brief Errors, held in ready to be accessed from outside the class
     * The array is filled during the function call NSE3D::output()
     * Currently, the errors store the L2 and H1-semi errors of the velocity
     * (errors.at(0) is L2 and errors.at(1) is H1-semi)
     * and the pressure (errors.at(2) is L2 and errors.at(3) is H1-semi).
     */
    std::array<double, int(4)> errors;
    
    /** @brief output object*/
    DataWriter3D outputWriter;
    
protected:
    
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
    
public:
    
    /** @brief Standard constructor of an Brinkman3D problem.
     *
     * @note The domain must have been refined a couple of times already if you want
     * to use multigrid. On the finest level the finite element spaces and
     * functions as well as matrices, solution and right hand side vectors are
     * initialized.
     *
     * @param collections
     * @param example The example to perform
     */
    Brinkman3D(std::list<TCollection* > collections,
               const ParameterDatabase& param_db,
               const Example_Brinkman3D& example,
               int reference_id = -4711);
    
    
    /** @brief constructor
     *
     * This constructor calls the other constructor creating an Example_Brinkman2D
     * object for you. See there for more documentation.
     */
    Brinkman3D(const TDomain& domain,
               const ParameterDatabase& param_db,
               const Example_Brinkman3D& example,
               int reference_id = -4711);
    
    
    /** @brief standard destructor */
    ~Brinkman3D();
    
    /** @brief assemble matrix */
    void assemble();
    
    /** @brief solve the system */
    void solve();
    
    /** @brief Use Petsc solver */
    void solve_with_Petsc(ParameterDatabase parmoon_db);

    
    /**
     * @brief Compute the defect Ax-b, and the residuals and store it all
     *
     * where A is the current matrix, x is the current solution and b is the
     * right hand side.
     */
    
    /** @brief compute the norms of the residuals */
    void compute_norm_of_residual();
    
    /** @brief check if one of the stopping criteria is fulfilled
     *
     * either converged, maximun number of iterations reached, or slow
     * convergence
     *
     * @param iteration_counter current iterate
     */
    void stopIt(unsigned int iteration_counter);
    
    /**
     * @brief measure errors and create pictures
     *
     * The current errors will be printed out. If desired, further output, e.g.,
     * vtk files are created.
     *
     * @param i suffix for output file name, -1 means no suffix
     */
    void output(int i = -1);
    
    /** @brief This function produces the output of some data like dofs and grid size on the console */
    void output_problem_size_info() const;
    
    /*******************************************************************************/
    // Declaration of special member functions - delete all but destructor.
    // This problem class will be used to request the whole process of
    // putting up and solving a single flow problem. It is, due to
    // its procedural nature, not supposed to be copied.
    
    //! Delete copy constructor. No copies allowed.
    Brinkman3D(const Brinkman3D&) = delete;
    
    //! Delete move constructor. No moves allowed.
    Brinkman3D(Brinkman3D&&) = delete;
    
    //! Delete copy assignment operator. No copies allowed.
    Brinkman3D& operator=(const Brinkman3D&) = delete;
    
    //! Default move assignment operator. No moves allowed.
    Brinkman3D& operator=(Brinkman3D&&) = delete;
    
    //    //! Default destructor. Does most likely cause memory leaks.
    //    ~Brinkman3D() = default;
    
    /********************************************************************************/
    // getters
    const TFEVectFunct3D& get_velocity() const
    { return this->systems.front().u; }
    
    TFEVectFunct3D& get_velocity()
    { return this->systems.front().u; }
    
    TFEFunction3D *get_velocity_component(int i);
    
    const TFEFunction3D& get_pressure() const
    { return this->systems.front().p; }
    
    TFEFunction3D& get_pressure()
    { return this->systems.front().p; }
    
    const TFESpace3D& get_velocity_space() const
    { return this->systems.front().velocity_space; }
    
    const TFESpace3D& get_pressure_space() const
    { return this->systems.front().pressure_space; }
    
    int get_size(){return this->systems.front().solution.length();}
    
    const ParameterDatabase & get_db() const
    { return brinkman3d_db; }
    
    /// @brief Get the current residuals  (updated in compute_residuals)
    const Residuals& get_residuals() const;
    /// @brief get the current impuls residual (updated in compute_residuals)
    double get_impuls_residual() const;
    /// @brief get the current mass residual (updated in compute_residuals)
    double get_mass_residual() const;
    /// @brief get the current residual (updated in compute_residuals)
    double get_full_residual() const;
    
    /// @brief return the computed errors (computed in output())
    std::array<double, int(4)> get_errors() const;
    
};


#endif /* Brinkman3D_H_ */
