﻿# CMakeLists.txt for subdirectory FE of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/FE"
                    "${CMAKE_SOURCE_DIR}/include/FE1D"
                    "${CMAKE_SOURCE_DIR}/include/FE2D"
                    "${CMAKE_SOURCE_DIR}/include/FE3D")

# Source files used in 2D and 3D.
#list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/ADICell.C") 
#list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/ADICell1D.C") 
#list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/ADISystem.C") 
#list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/ADISystem1D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/AlgebraicFluxCorrection.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/Assemble2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/AuxParam2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/BaseFunct1D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/BaseFunct2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/Bcgs.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/Blas1.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/Blas2.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/Bulk.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/Cg.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/ConvDiff.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/Convolution.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/DirectSolver.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/DiscreteForm2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FE1D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FE2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FE2DMapper.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FE2DMapper1Reg.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FEDatabase2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FEDesc1D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FEDesc2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FEFunction1D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FEFunction2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FEMatrix.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FEM_TVD_FCT.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FESpace.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FESpace1D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FESpace2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FgmresIte.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FixedPointIte.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/HNDesc.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/HangingNode.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/ItMethod.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/JacobiIte.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/LinAlg.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/LineAffin.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/LocalProjection.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/MGComponents2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/Matrix.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/MultiGridIte.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/MultiGridScaIte.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/NSE_MGLevel.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/NSE_MGLevel1.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/NSE_MGLevel2.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/NSE_MGLevel3.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/NSE_MGLevel4.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/NSE_MGLevel14.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/NSE_MultiGrid.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/NodalFunctional1D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/NodalFunctional2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/Output2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/PardisoSolver.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/QuadAffin.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/QuadIsoparametric.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/QuadBilinear.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/RationalLES.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/RefTrans1D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/RefTrans2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/Solver.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/SquareMatrix.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/SquareMatrix1D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/SquareMatrix2D.C")
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/SSORIte.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/Structure.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/TimeDiscRout.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/TriaAffin.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/TriaIsoparametric.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/VMS.C") 

# Source files only used in 2D
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/Assemble2D_edge_Oseen.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/Bulk_2d3d.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/CD2DErrorEstimator.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/ConvDiff2D.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/Drop.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/Darcy2DMixed.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/FEVectFunct2D.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/FreeSurface2D.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/LocalAssembling2D.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/MGLevel2D.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/MovingNavierStokes.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/MultiGrid2D.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/Matrix2D.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/NSE2DErrorEstimator.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/NSE2D_AxialSymm3D_FixPo.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/NSE2D_EquOrd_FixPo.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/NSE2D_FixPo.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/NSE2D_FixPoSkew.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/NSE2D_FixPoRot.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/NSE2D_Friction_FixPo.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/NSE2D_Newton.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/PDAE_Index2D_2.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/RFB.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/RKV_FDM.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/SSMUM.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/Superconvergence2D.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/TCD2D.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/TimeConvDiff2D.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/TNSE2D_FixPo.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/TNSE2D_FixPoRot.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/TNSE2D_FixPo_SSMUM.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/Upwind.C") 

# Source files only used in 3D
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/Assemble3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/Aux2D3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/AuxParam3D.C")  
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/BaseFunct3D.C")  
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/Bulk_3d4d.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/ConvDiff3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/DiscreteForm3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FE3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FE3DMapper.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FE3DMapper1Reg.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FEDatabase3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FEDesc3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FEFunction3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FESpace3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FEVectFunct3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FreeSurface3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/HexaAffin.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/HexaIsoparametric.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/HexaTrilinear.C")
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/LocalAssembling3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/MGComponents3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/MGLevel3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/Matrix3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/MultiGrid3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/NSE3D_FixPo.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/NSE3D_FixPoSkew.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/NSE3D_Friction_FixPo.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/NSE3D_Newton.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/NodalFunctional3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/Output3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/SquareMatrix3D.C")
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/Superconvergence3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/TCD3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/TNSE3D_FixPo.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/TNSE3D_Newton.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/RefTrans3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/RFB.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/TetraAffin.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/TetraIsoparametric.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/Upwind3D.C") 
#list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/Urea_3d4d.C") 


list(APPEND PARMOON_SOURCES_2D ${FE_SOURCES} ${FE_SOURCES_2D})
list(APPEND PARMOON_SOURCES_3D ${FE_SOURCES} ${FE_SOURCES_3D})
