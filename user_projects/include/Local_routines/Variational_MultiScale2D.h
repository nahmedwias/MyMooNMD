// =======================================================================
// Variational_MultiScale2D.h
//
// Purpose:     routines for projection-based VMS
//
// Author:       Najib Alia  2017/03/24
//
// =======================================================================


#ifndef _VARIATIONAL_MULTISCALE2D_H_
#define _VARIATIONAL_MULTISCALE2D_H_

#include <memory>
#include<FEMatrix.h>
#include<BlockFEMatrix.h>
#include<FEVectFunct2D.h>


/** creates a diagonal mass matrix */
void LumpMassMatrixToDiagonalMatrix2D(std::shared_ptr<FEMatrix> & matrix);

/**
 *  block=
 *  [ A11  A12  B1T ]
 *  [ A21  A22  B2T ]
 *  [ B1   B2   C   ]
 */

/**
 * matrices_vms[] ={GT11, GT22, G11, G22, M};
 */
void VMS_ProjectionUpdateMatrices2D(std::vector< std::shared_ptr< FEMatrix > >& blocks,
                         std::array< std::shared_ptr< FEMatrix >, int(5) > matrices_vms);

//TODO:
void ComputeVMSProjection2D(const std::array< std::shared_ptr< FEMatrix >, int(4) > matrixG,
                          const TFEVectFunct2D& velocity, TFEVectFunct2D& projection);



#endif /* _VARIATIONAL_MULTISCALE2D_H_ */
