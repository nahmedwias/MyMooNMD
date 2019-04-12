/** **************************************************************************** 
*
* @name       TimeConvectionDiffusionPOD
* @brief      Computation of POD basis - specific routines for TCD problems
*
* @author     Swetlana Giere & Alfonso Caiazzo
* @date       08.03.2017 (start of implementation), 15.1.2019 (restart)
*
*******************************************************************************/

#ifndef TIMECONVECTIONDIFFUSIONPOD_H
#define TIMECONVECTIONDIFFUSIONPOD_H

// boost
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>


#include "POD.h"

#include "BlockFEMatrix.h"
#include "BlockVector.h"
#ifdef __2D__
#include "Example_TimeCD2D.h"
#include "FEFunction2D.h"
#else
#include "Example_TimeCD3D.h"
#include "FEFunction3D.h"
#endif

#include "DataWriter.h"
#include "templateNames.h"
#include "TimeDiscretizations.h"
#include "AlgebraicFluxCorrection.h"
#include "LocalAssembling.h"

namespace ublas = boost::numeric::ublas;

template<int d>
class TimeConvectionDiffusionPOD : public POD
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_TimeCD = typename Template_names<d>::Example_TimeCD;

    /** @brief constructor
    * this constructor calls the other constructor creating an
    * Example_TimeCD object
    */
    TimeConvectionDiffusionPOD(TCollection&             coll,
                               const ParameterDatabase& param_db);

    /** @brief The standard constructor, can be used for multigrid and 
    * non-multigrid.
    *
    * @param[in] collection collection of cells
    * @param[in] param_db A parameter database with parameters concerning this
    *                     class or any of its members (fe space, solver,
    *                     assemble,...)
    * @param[in] example The example which is to be calculated.
    */
    TimeConvectionDiffusionPOD(TCollection&             coll,
                               const ParameterDatabase& param_db, 
                               const Example_TimeCD&    ex);

    /** TODO: adapt comment (if more than tcd)
    * @brief return a database with all parameters necessary for
    * time-dependent convection-diffusion (tcd) probems
    */
    static ParameterDatabase set_pod_basis_database(
                                             const ParameterDatabase& param_db);

    /**
    * @brief Compute POD basis
    *
    * Compute POD basis from snapshots with the compute_basis() routine
    * from the class POD. Within the function, appropriate gramian matrix
    * according to the 'pod_inner_product' parameter is assembled and
    * incorporated into computation of the POD basis.
    */
    void compute_pod_basis();

    /**
    * @brief Read POD basis from file
    *
    * Read POD basis from file. The functionality could be interesting if POD
    * basis is already computed and stored into file and one
    * only wants to write vtk files of the POD basis functions.
    */
    void read_basis();

    /// @brief write POD basis to file and write vtk files
    void output();

  protected:
    /** @brief Finite Element space */
    std::shared_ptr<FESpace> fe_space;

    /** @brief gramian matrix (needed for POD computation) */
    BlockFEMatrix gramian_matrix;

    BlockVector pod_mode;

    /** @brief Finite Element function */
    FEFunction fe_function;

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

    /** @brief Assemble gramian matrix
    *
    * Assemble gramian matrix, which describes the inner product, with respect
    * to which the POD basis will be computed.
    */
    void assemble_gramian();

    /** @brief print the problem infomation
    */
    void output_problem_size_info() const;

    /** @brief a local parameter database which controls this class
    *
    * The database given to the constructor will be merged into this one. Only
    * parameters which are of interest to this class are stored (and the
    * default ParMooN parameters). Note that this usually does not include
    * other parameters such as solver parameters. Those are only in the
    * Solver object.
    */
    ParameterDatabase db;

    /// @brief time stepping scheme object to access everything
    TimeDiscretization time_stepping_scheme;

    /** @brief Definition of the used example */
    const Example_TimeCD example;

    /** @brief class for output handling */
    DataWriter<d> outputWriter;
};

#endif // TIMECONVECTIONDIFFUSIONPOD_H
