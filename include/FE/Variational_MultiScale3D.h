// =======================================================================
// Variational_MultiScale3D.h
//
// Purpose:     routines for projection-based VMS
//
// Author:       Volker John  2006/05/18
// 
// =======================================================================

#ifndef __Variational_MultiScale3D__
#define __Variational_MultiScale3D__

#include <memory>
#include<FEMatrix.h>
#include<BlockFEMatrix.h>
#include<FEVectFunct3D.h>

/** creates a diagonal mass matrix */
void LumpMassMatrixToDiagonalMatrix3D(std::shared_ptr<FEMatrix> & matrix);

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
void VMS_ProjectionUpdateMatrices3D(std::vector< std::shared_ptr< FEMatrix > >& blocks, 
                         std::array< std::shared_ptr< FEMatrix >, int(7) > matrices_vms);

//TODO:
void ComputeVMSProjection(const std::array< std::shared_ptr< FEMatrix >, int(4) > matrixG, 
                          const TFEVectFunct3D& velocity, TFEVectFunct3D& projection);
#endif
