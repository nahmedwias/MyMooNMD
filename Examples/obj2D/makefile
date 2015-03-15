

RELPATH = ../..
SRCPATH = 
LIBPATH = ../../lib2D/

INC    = -I.. \
         -I../../include/General -I../../include/FE \
         -I../../include/Geometry -I../../include/QuadFormulas \
         -I../../include/FE1D -I../../include/FE2D \
         -I../../include/AMG \
         -I../../include/Parallel \
         -I../../EXT_LIB/UMFPACK/Include \
         -I../../EXT_LIB/tecplot/include

include ../../Make.inc

LIBRARY = $(LIBPATH)lib_Examples.a
   MEMBERS = \
	$(SRCPATH)Example2D.o \
	$(SRCPATH)Example_CD2D.o \
	$(SRCPATH)Example_Darcy2D.o \
	$(SRCPATH)Example_NSE2D.o \

###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) 
	