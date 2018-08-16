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
list(APPEND ASSEMBLESADDLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/TNSE2DGalerkin.C")
list(APPEND ASSEMBLESADDLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/TNSE2DSmagorinsky.C")
list(APPEND ASSEMBLESADDLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/NSE2DGalerkin.C")
list(APPEND ASSEMBLESADDLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/NSE2DSUPG.C")
list(APPEND ASSEMBLESADDLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/CommonRoutineTNSE2D.C")
list(APPEND ASSEMBLESADDLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/TNSE2DSUPG.C")
list(APPEND ASSEMBLESADDLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/TNSE2DResBasedVMS.C")


list(APPEND ASSEMBLESADDLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/NSE3DGalerkin.C")
list(APPEND ASSEMBLESADDLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/TNSE3DSmagorinsky.C")
list(APPEND ASSEMBLESADDLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/CommonRoutineTNSE3D.C")
list(APPEND ASSEMBLESADDLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/TNSE3DProjBasedVMS.C")
list(APPEND ASSEMBLESADDLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/TNSE3DSUPG.C")
list(APPEND ASSEMBLESADDLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/TNSE3DResBasedVMS.C")

list(APPEND ASSEMBLESADDLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/Variational_MultiScale3D.C")

list(APPEND PARMOON_SOURCES_2D ${ASSEMBLESADDLE_SOURCES_2D})
list(APPEND PARMOON_SOURCES_3D ${ASSEMBLESADDLE_SOURCES_3D})
