

RELPATH = ../..
SRCPATH = 
LIBPATH = ../../lib3D/

INC    = -I.. \
         -I../../include/General -I../../include/FE \
         -I../../include/Geometry -I../../include/QuadFormulas \
         -I../../include/FE1D -I../../include/FE2D -I../../include/FE3D\
         -I../../include/AMG \
         -I../../include/Parallel \
         -I../../EXT_LIB/UMFPACK/Include \
         -I../../EXT_LIB/tecplot/include

include ../../Make.inc

LIBRARY = $(LIBPATH)lib_Examples3D_par.a
   MEMBERS = \
	$(SRCPATH)Example3D.o \
	$(SRCPATH)Example_CD3D.o \
	$(SRCPATH)Example_NSE3D.o \

###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) 
	
