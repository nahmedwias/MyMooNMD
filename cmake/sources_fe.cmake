# CMakeLists.txt for subdirectory FE of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
# 2016/09/23 Najib Alia: added my own files for user_projects

# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/FE"
                    "${CMAKE_SOURCE_DIR}/include/FE1D"
                    "${CMAKE_SOURCE_DIR}/include/FE2D"
                    "${CMAKE_SOURCE_DIR}/include/FE3D"
		    "${CMAKE_SOURCE_DIR}/user_projects/include/FE"
		    "${CMAKE_SOURCE_DIR}/user_projects/include/Local_routines")

# Source files used in 2D and 3D.
#list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/ADICell.C") 
#list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/ADICell1D.C") 
#list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/ADISystem.C")
#list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/ADISystem1D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/AuxParam2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/BaseFunct1D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/BaseFunct2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/Blas1.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/Blas2.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/Convolution.C") 
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
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FESpace.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FESpace1D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/FESpace2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/HNDesc.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/HangingNode.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/LinAlg.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/LineAffin.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/LocalProjection.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/NodalFunctional1D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/NodalFunctional2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/QuadAffin.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/QuadIsoparametric.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/QuadBilinear.C") 
#list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/RationalLES.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/RefTrans1D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/RefTrans2D.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/TimeDiscRout.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/TriaAffin.C") 
list(APPEND FE_SOURCES "${PROJECT_SOURCE_DIR}/src/FE/TriaIsoparametric.C") 


# Source files only used in 2D
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/FEVectFunct2D.C")
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/PostProcessing2D.C") 
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/FE/Superconvergence2D.C")

list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/src/FE/FEFunctionInterpolator.C")
list(APPEND FE_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/src/Local_routines/TLinElastic2D_routines.C")

# Source files only used in 3D
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/Aux2D3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/AuxParam3D.C")  
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/BaseFunct3D.C")  
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/DiscreteForm3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FE3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FE3DMapper.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FE3DMapper1Reg.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FEDatabase3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FEDesc3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FEFunction3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FESpace3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/FEVectFunct3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/HexaAffin.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/HexaIsoparametric.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/HexaTrilinear.C")
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/NodalFunctional3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/Output3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/Superconvergence3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/RefTrans3D.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/TetraAffin.C") 
list(APPEND FE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/FE/TetraIsoparametric.C")




list(APPEND PARMOON_SOURCES_2D ${FE_SOURCES} ${FE_SOURCES_2D})
list(APPEND PARMOON_SOURCES_3D ${FE_SOURCES} ${FE_SOURCES_3D})
