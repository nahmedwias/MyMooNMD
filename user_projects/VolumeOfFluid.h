#ifndef VOLUMEOFFLUID_H
#define VOLUMEOFFLUID_H

#include "vector"
#include "../include/General/templateNames.h"
#include "../include/Matrix/BlockFEMatrix.h"
#include "../include/Matrix/BlockVector.h"
#include "../../include/General/ParameterDatabase.h"

#ifdef __2D__
#include "../include/FE/FEFunction2D.h"
#include "../include/FE/FEVectFunct2D.h"
#else
#include "../include/FE/FEFunction3D.h"
#include "../include/FE/FEVectFunct3D.h"
#endif
/**
 * @todo write docs
 */

template <int d>
class VolumeOfFluid
{
public:
  using FEFunction = typename Template_names<d>::FEFunction;
  using FEVectFunct = typename Template_names<d>::FEVectFunct;
  using FESpace = typename Template_names<d>::FESpace;
  using Example_TimeNSE = typename Template_names<d>::Example_TimeNSE;
  
protected:
  /** @brief see the other constructor */
  VolumeOfFluid(const TDomain& domain, const ParameterDatabase& param_db);
  /** @brief constructor
     *
     * The domain must have been refined a couple of times already. On the
     * finest level the finite element spaces and functions as well as
     * matrices, solution and right hand side vectors are initialized.
     * 
     * @param domain the computational domain to get the grid(s)
     * @param param_db parameters controlling this class
     * @param example The example to use
     */
  VolumeOfFluid(const TDomain& domain, const ParameterDatabase& param_db,
                     const Example_TimeNSE& ex);
  
  
private:
  struct System_per_grid
  {
    /** @brief Finite element space for velocity*/
    std::shared_ptr<FESpace> velocity_space;
    /** @brief Finite element space for pressure*/
    std::shared_ptr<FESpace> pressure_space;
    /** @brief Finite element space for pressure*/
    std::shared_ptr<FESpace> temp_space;
    
    /** @brief system matrix                     |    [ A11  A12  A13  B1T ]
       *                       [ A11  A12  B1 ]    |    [ A21  A22  A23  B2T ]
       *                       [ A21  A22  B2 ]    |    [ A31  A32  A33  B3T ]
       *                       [ B3   B4   C  ]    |    [ B1   B2   B3   C   ]
      */
    BlockFEMatrix matrix_NS;
    /***/
    BlockFEMatrix matrix_temp;
    /** @brief mass matrix: this will be the standard mass matrix
       * for the standard Galerkin scheme. However, for the SUPG/RBVMS
       * schems, this includes all the terms which are related to the
       * time derivatives in the fully discrete scheme Eq: (45)
       * Ahmed, Rebollo, John and Rubino (2015)
       *  [ M11  M12  0 ]
       *  [ M21  M22  0 ]
       *  [ 0     0   0 ]
       * This additionaly created due to the struture of the
       * residual based Variational Multiscale Method:
       */
    BlockFEMatrix mass_matrixNS;
    /****/
    BlockFEMatrix mass_matrixtemp;
    
    /** @brief right hand side vector*/
    BlockVector rhs_NS;
    /** @brief solution vector*/
    BlockVector solution_NS;
    /** @brief Finite element function for velocity*/
    FEVectFunct u;
    /** @brief Finite element function for pressure*/
    FEFunction p;
    
    /***/
    BlockVector rhs_temp;
    BlockVector sol_temp;
    
    /** @brief Finite element function for temprature*/
    FEFunction temp;
    
    /** @brief constructor*/
    System_per_grid(const Example_TimeNSE& example, TCollection& coll,
                      std::tuple<int,int,int> order);
    
  };
  
};

#endif // VOLUMEOFFLUID_H
