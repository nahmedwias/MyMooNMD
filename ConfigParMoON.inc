#=============================================================================
# configuration file for ParMooN
# modify this file according to your machine type and requirement 
# after that use one of the commands in the terminal to compile ParMooN
#=============================================================================
# make 2D -j       (to compile  2D sequential programs)
# make 3D -j       (to compile  3D sequential programs)
# make 3D_par -j   (to compile  3D parallel programs)

#=============================================================================
# include the main programe file (without .C)
# with full path from the current directory
#=============================================================================
#IN = 3DPrograms/NSE3D_ParMooN
#IN = 2DPrograms/CD2D_ParMooN
IN = 3DPrograms/CD3D_ParMooN
#=============================================================================
# provide the full path from the current directory
# to where the exe file want to be stored
#=============================================================================
OUT=../ParMooN_Output/cylinder/nse3d_mpi.exe
OUT=../ParMooN_Output/cd3d/cd3d_mpi.exe

#=============================================================================
#select the machine type
#=============================================================================
#MACH = LINUX
#MACH = INTEL64 
#MACH = TYRONE64  
MACH = LINUX64
MACH = MAC64
#MACH = CRAY64

#=============================================================================
# select the compilation type
#=============================================================================
FLAG = OPTI
FLAG = FAST
#FLAG = DEBUG

#=============================================================================
# select one of the parallel type: 
# hybrid (takes all, mpi, hybrid, openmp)
# MPI     (takes both mpi and openmp)
# SEQUENTIAL (compilation takes openmp)
# OMPONLY = (takes omp only)
#=============================================================================
PARALLEL_TYPE = SEQUENTIAL
PARALLEL_TYPE = HYBRID
PARALLEL_TYPE = MPI
#PARALLEL_TYPE = OMPONLY

#=============================================================================

