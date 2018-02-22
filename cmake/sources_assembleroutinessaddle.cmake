# CMakeLists.txt for subdirectory AssembleRoutinesSaddle of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2017/03/22 Najib Alia: Rework to separate FE and AssembleRoutines files.
# 2017/03/22 Naveed Ahmed: splitting of the local assembling routnies for the saddle point problems
#

# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/AssembleRoutinesSaddle")

# Source files used in 2D and 3D.
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/TNSE2DGalerkin.C")
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/TNSE2DSmagorinsky.C")
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/NSE2DGalerkin.C")
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/CommonRoutineTNSE2D.C")
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/TNSE2DSUPG.C")
list(APPEND ASSEMBLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/TNSE2DResBasedVMS.C")


list(APPEND ASSEMBLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/TNSE3DSmagorinsky.C")
list(APPEND ASSEMBLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/CommonRoutineTNSE3D.C")
list(APPEND ASSEMBLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/TNSE3DProjBasedVMS.C")
list(APPEND ASSEMBLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/Variational_MultiScale3D.C")

list(APPEND PARMOON_SOURCES_2D ${ASSEMBLE_SOURCES_2D})
list(APPEND PARMOON_SOURCES_3D ${ASSEMBLE_SOURCES_3D})
