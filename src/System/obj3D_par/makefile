###################################
#
#  Makefile for ParMooN
#  author: Sashikumaar Ganesan
#
###################################

RELPATH = ../../..
SRCPATH = 
LIBPATH = ../../../lib3D/
INC    = -I../../../include/General -I../../../include/FE \
         -I../../../include/Geometry -I../../../include/QuadFormulas \
         -I../../../include/FE2D -I../../../include/FE1D \
         -I../../../include/FE3D -I../../../include/AMG \
         -I../../../include/PBE -I../../../EXT_LIB/UMFPACK/Include \
         -I../../../include/System \
         -I../../../include/Parallel  -I../../../EXT_LIB/MUMPS/include \
         -I../../../EXT_LIB/tecplot/include

include ../../../Make.inc

LIBRARY = $(LIBPATH)lib_System_par.a
   MEMBERS = \
        $(SRCPATH)SystemMatScalar3D.o \
        $(SRCPATH)SystemMatTimeScalar3D.o \
        $(SRCPATH)SystemMatNSE3D.o 
###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) ii_files SunWS_cache

