#ifndef PRESROBUSTNAVIERSTOKES_H
#define PRESROBUSTNAVIERSTOKES_H

#include "../include/System/TimeNavierStokes.h"

#include <../include/Geometry/Domain.h>
#include <../include/Matrix/BlockFEMatrix.h>
#include <BlockFEMatrixPr.h>

#ifdef __2D__
#include "Example_NSPR_NSE2D.h"
#include "../include/FE/FEFunction2D.h"
#include "../include/FE/FEVectFunct2D.h"
#else
#endif
#include "../include/Matrix/BlockVector.h"
#include "../include/General/ParameterDatabase.h"
#include "../include/Solver/Solver.h"
#include "../include/General/templateNames.h"


/**
 * @todo write docs
 */
template <int d>
class PresRobustNavierStokes : public TimeNavierStokes<d>
{
public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FEVectFunct = typename Template_names<d>::FEVectFunct;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_NSE = typename Template_names<d>::Example_NSE;
    /**
     * Default constructor
     */
    PresRobustNavierStokes(const TDomain& domain, const ParameterDatabase& param_db, 
            const Example_NSPR_NSE2D& ex, std::tuple<int,int,int>_spaces_orders_);

    /** assemble function */
    void assemble_matrices();
    
protected:
    /** @brief Finite Element space for projection */
    FESpace projecti_space;
    /** @brief projection matrix 
     * the same structure as the matrix of NS system
     */
    BlockFEMatrixPr pr_mat;
    
    /** @brief 
     * Modified mass matrix after reconstruction
     **/
    BlockFEMatrix mass_reconst_;
    
    /** @brief right hand side for reconstruction*/
    BlockVector rhs_xh;
    /** @brief a local parameter database which controls this class
     * 
     *
     * The database given to the constructor will be merged into this one. Only
     * parameters which are of interest to this class are stored (and the
     * default ParMooN parameters). Note that this usually does not include
     * other parameters such as solver parameters. Those are only in the
     * Solver object.
     */
    ParameterDatabase db;
    
    /** @brief Definition of the used example. */
    const Example_NSPR_NSE2D example;
    
};

#endif // PRESROBUSTNAVIERSTOKES_H
