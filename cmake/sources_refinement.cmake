# CMakeLists.txt for subdirectory Refinement of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/Refinement")

list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/Refinement.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefLineDesc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefNoRef.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefQuad1Conf0Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefQuad1Conf1Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefQuad1Conf2Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefQuad1Conf3Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefQuad2Conf0Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefQuad2Conf1Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefQuad2Conf2Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefQuad2Conf3Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefQuadBis0Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefQuadBis1Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefQuadRegDesc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefQuadToTri0Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefQuadToTri1Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTriBaryDesc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTriBis01Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTriBis02Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTriBis0Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTriBis10Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTriBis12Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTriBis1Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTriBis20Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTriBis21Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTriBis2Desc.C")
list(APPEND REF_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTriRegDesc.C")

list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefHexaRegDesc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBaryDesc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis01Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis02Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis03Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis04Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis05Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis0Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis10Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis12Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis13Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis14Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis15Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis1Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis20Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis21Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis23Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis24Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis25Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis2Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis30Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis32Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis34Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis35Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis3Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis40Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis41Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis43Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis45Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis4Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis51Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis52Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis53Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis54Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraBis5Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraQuad0Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraQuad1Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraQuad2Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraQuad3Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraReg0Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraReg1Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraReg2Desc.C")
list(APPEND REF_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Refinement/RefTetraRegDesc.C")



list(APPEND PARMOON_SOURCES_2D ${REF_SOURCES_2D})
list(APPEND PARMOON_SOURCES_3D ${REF_SOURCES_2D} ${REF_SOURCES_3D})
