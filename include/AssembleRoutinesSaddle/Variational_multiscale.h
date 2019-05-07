// =======================================================================
// Variational_multiScale.h
//
// Purpose:     routines for projection-based VMS
// =======================================================================

#ifndef __Variational_multiScale__
#define __Variational_multiScale__

#include <memory>
#include<FEMatrix.h>
#include<BlockFEMatrix.h>

#ifdef __2D__
#include <FEVectFunct2D.h>
#endif

#ifdef __3D__
#include<FEVectFunct3D.h>
#endif

#include "templateNames.h"
#include "vector"

/** creates a diagonal mass matrix */
template<int d>
void LumpMassMatrixToDiagonalMatrix(std::shared_ptr<FEMatrix> & matrix);

/**
 *  block=
 *  [ A11  A12  A13  B1T ]
 *  [ A21  A22  A23  B2T ]
 *  [ A31  A32  A33  B3T ]
 *  [ B1   B2   B3   C   ]
 */

/**
 * matrices_vms[] ={Gt11, GT22, GT33, G11, G22, G33, M};
 */
template<int d>
void VMS_ProjectionUpdateMatrices(std::vector< std::shared_ptr< FEMatrix > >& blocks,
                                  std::vector< std::shared_ptr< FEMatrix > > matrices_vms); 
#endif
