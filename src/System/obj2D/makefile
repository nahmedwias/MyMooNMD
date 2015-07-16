###################################
#
#  Makefile for ParMooN
#  author: Sashikumaar Ganesan
#
###################################

RELPATH = ../../..
SRCPATH = 
LIBPATH = ../../../lib2D/
INC    = -I../../../include/General -I../../../include/FE \
         -I../../../include/Geometry -I../../../include/QuadFormulas \
         -I../../../include/FE1D -I../../../include/FE2D \
         -I../../../include/AMG \
         -I../../../include/System -I../../../Examples \
         -I../../../include/Parallel \
         -I../../../EXT_LIB/UMFPACK/Include \
         -I../../../EXT_LIB/tecplot/include


include ../../../Make.inc

DEFINEPARAMS = ../../../include/General/DefineParams.h

LIBRARY = $(LIBPATH)lib_System.a
   MEMBERS = \
	$(SRCPATH)CD2D.o \
	$(SRCPATH)SystemMat2D.o \
	$(SRCPATH)SystemDarcy2D.o \
	$(SRCPATH)SystemTCD2D.o \
	$(SRCPATH)SystemTCD2D_ALE.o \
	$(SRCPATH)SystemMatNSE2D.o \
	$(SRCPATH)SystemMatScalar2D.o \
	$(SRCPATH)SystemMatDarcy2D.o \
	$(SRCPATH)NSE2D.o \
	$(SRCPATH)SystemTNSE2D.o \
	$(SRCPATH)SystemMatTimeScalar2D.o \
	$(SRCPATH)SystemMatTNSE2D.o \
###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) ii_files SunWS_cache


