# CMakeLists.txt for subdirectory AssembleRoutines of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2017/03/22 Najib Alia: Rework to separate FE and AssembleRoutines files.
#

# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/AssembleRoutines")

# Source files used in 2D and 3D.
list(APPEND ASSEMBLE_SOURCES "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/AlgebraicFluxCorrection.C") 
list(APPEND ASSEMBLE_SOURCES "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/Assemble2D.C")  
list(APPEND ASSEMBLE_SOURCES "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/Assembler4.C")  
list(APPEND ASSEMBLE_SOURCES "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/BoundaryAssembling2D.C")
list(APPEND ASSEMBLE_SOURCES "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/Bulk.C")
list(APPEND ASSEMBLE_SOURCES "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/ConvDiff.C") 


# Source files only used in 2D
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/Assemble2D_edge_Oseen.C") 
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/Brinkman2D_Mixed.C")
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/CD2DErrorEstimator.C") 
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/ConvDiff2D.C") 
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/Drop.C") 
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/DarcyMixed.C")
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/LocalAssembling.C")
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/TCD2D.C")
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/TCD2D_POD_LocalAssemble.C")
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/TimeConvDiff2D.C")
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/TNSE3D_FixPo.C")
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/Upwind.C")

# Source files only used in 3D
list(APPEND ASSEMBLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/Assemble3D.C") 
list(APPEND ASSEMBLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/ConvDiff3D.C")
list(APPEND ASSEMBLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/DarcyMixed.C")
list(APPEND ASSEMBLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/LocalAssembling.C")
list(APPEND ASSEMBLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/TCD3D.C") 
list(APPEND ASSEMBLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/TNSE3D_FixPo.C") 
list(APPEND ASSEMBLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/TNSE3D_Newton.C") 
list(APPEND ASSEMBLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/Upwind3D.C") 


list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/CD_local_assembling_routines.C") 
list(APPEND ASSEMBLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutines/CD_local_assembling_routines.C")

list(APPEND PARMOON_SOURCES_2D ${ASSEMBLE_SOURCES} ${ASSEMBLE_SOURCES_2D})
list(APPEND PARMOON_SOURCES_3D ${ASSEMBLE_SOURCES} ${ASSEMBLE_SOURCES_3D})
