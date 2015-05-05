#############################
#GEO = 2D
GEO = 3D

IN = main$(GEO)
OUT = main

#MACH = LINUX
#MACH = INTEL64 
#MACH = TYRONE64  
MACH = LINUX64
#MACH = MAC64

#FLAG = OPTI
FLAG = FAST
#FLAG = DEBUG

# profing options
# HP-UX: -G, SUN: -xpg,
# other platform: -pg
#GPROF  = -pg

#EFENCE = YES

#MORTAR = -D__MORTAR__
MORTAR =

#PARALLEL_TYPE = HYBRID
PARALLEL_TYPE = MPI
#PARALLEL_TYPE = OPENMP
#PARALLEL_TYPE = OMPONLY

#hybrid  = takes all, mpi, hybrid, openmp
#MPI     = takes both mpi & openmp
#openmp  = takes only openmp
#seq compilation takes openmp
#omponly = takes omp only
#############################

SRC = src

MAKE = gmake

ARGS = MACH=$(MACH) FLAG=$(FLAG) MORTAR=$(MORTAR) \
        OUT=$(OUT) IN=$(IN) GPROF=$(GPROF) \
       PARALLEL_TYPE=SEQUENTIAL PRG_TYPE=SEQUENTIAL 

ARGS_PAR = MACH=$(MACH) FLAG=$(FLAG) MORTAR=$(MORTAR) \
        OUT=$(OUT) IN=$(IN) GPROF=$(GPROF) \
       PARALLEL_TYPE=$(PARALLEL_TYPE) PRG_TYPE=PARALLEL 

RM = rm -rf

all: $(GEO)
	@echo all

clena: clean
	@echo clena

clean: clean2D clean3D
	@echo clean

clean_2d: clean2D
	@echo clean

clean_3d: clean2D clean3D
	@echo clean

clean_par: clean2D clean3D_PAR
	@echo clean

allclean: clean2D clean2D_PAR clean3D clean3D_PAR
	@echo clean

2D:
	@echo 2D
	@($(RM) $(OUT))
	@(cd $(SRC)/AMG/obj2D; $(MAKE) GEO=2D $(ARGS) )
	@(cd $(SRC)/Refinement/obj2D; $(MAKE) GEO=2D $(ARGS) )
	@(cd $(SRC)/QuadFormulas/obj2D; $(MAKE) GEO=2D $(ARGS) )
	@(cd $(SRC)/Geometry/obj2D; $(MAKE) GEO=2D $(ARGS) )
	@(cd $(SRC)/General/obj2D; $(MAKE) GEO=2D $(ARGS) )
	@(cd $(SRC)/FE/obj2D; $(MAKE) GEO=2D $(ARGS) )
	@(cd Examples/obj2D; $(MAKE) GEO=2D $(ARGS) )
	@(cd $(SRC)/System/obj2D; $(MAKE) GEO=2D $(ARGS) )
	@(cd $(SRC)/PBE/obj2D; $(MAKE) GEO=2D $(ARGS) )
	@$(MAKE) -f Makefile2D GEO=2D EFENCE=$(EFENCE) $(ARGS)


2D_par:
	@echo 2D
	@($(RM) $(OUT))
	@(cd $(SRC)/AMG/obj2D_par; $(MAKE) GEO=2D $(ARGS_PAR) )
	@(cd $(SRC)/Refinement/obj2D_par; $(MAKE) GEO=2D $(ARGS_PAR) )
	@(cd $(SRC)/QuadFormulas/obj2D_par; $(MAKE) GEO=2D $(ARGS_PAR) )
	@(cd $(SRC)/Geometry/obj2D_par; $(MAKE) GEO=2D $(ARGS_PAR) )
	@(cd $(SRC)/General/obj2D_par; $(MAKE) GEO=2D $(ARGS_PAR) )
	@(cd $(SRC)/FE/obj2D_par; $(MAKE) GEO=2D $(ARGS_PAR) )
	@(cd $(SRC)/System/obj2D_par; $(MAKE) GEO=2D $(ARGS_PAR) )
	@(cd Examples/obj2D_par; $(MAKE) GEO=2D $(ARGS) )
	@(cd $(SRC)/PBE/obj2D_par; $(MAKE) GEO=2D $(ARGS_PAR) )
	@(cd $(SRC)/Parallel/obj2D_par; $(MAKE) GEO=2D $(ARGS_PAR) )
	@$(MAKE) -f Makefile2D GEO=2D EFENCE=$(EFENCE) $(ARGS_PAR)

3D:
	@echo 3D
	@($(RM) $(OUT))
	@(cd $(SRC)/AMG/obj3D; $(MAKE) GEO=3D $(ARGS) )
	@(cd $(SRC)/Refinement/obj3D; $(MAKE) GEO=3D $(ARGS) )
	@(cd $(SRC)/QuadFormulas/obj3D; $(MAKE) GEO=3D $(ARGS) )
	@(cd $(SRC)/Geometry/obj3D; $(MAKE) GEO=3D $(ARGS) )
	@(cd $(SRC)/General/obj3D; $(MAKE) GEO=3D $(ARGS) )
	@(cd $(SRC)/FE/obj3D; $(MAKE) GEO=3D $(ARGS) )
	@(cd $(SRC)/System/obj3D; $(MAKE) GEO=3D $(ARGS) )
	@(cd $(SRC)/PBE/obj3D; $(MAKE) GEO=3D $(ARGS) )
	@$(MAKE) -f Makefile3D GEO=3D EFENCE=$(EFENCE) $(ARGS)

3D_par:
	@echo 3D_par
	@($(RM) $(OUT))
	@(cd $(SRC)/AMG/obj3D_par; $(MAKE) GEO=3D $(ARGS_PAR) )
	@(cd $(SRC)/Refinement/obj3D_par; $(MAKE) GEO=3D $(ARGS_PAR) )
	@(cd $(SRC)/QuadFormulas/obj3D_par; $(MAKE) GEO=3D $(ARGS_PAR) )
	@(cd $(SRC)/Geometry/obj3D_par; $(MAKE) GEO=3D $(ARGS_PAR) )
	@(cd $(SRC)/General/obj3D_par; $(MAKE) GEO=3D $(ARGS_PAR) )
	@(cd $(SRC)/FE/obj3D_par; $(MAKE) GEO=3D $(ARGS_PAR) )
	@(cd $(SRC)/PBE/obj3D; $(MAKE) GEO=3D $(ARGS_PAR) )
	@(cd $(SRC)/Parallel/obj3D_par; $(MAKE) GEO=3D $(ARGS_PAR) )
	@(cd $(SRC)/System/obj3D_par; $(MAKE) GEO=3D $(ARGS_PAR) )
	@$(MAKE) -f Makefile3D GEO=3D EFENCE=$(EFENCE) $(ARGS_PAR)

clean2D:
	@echo clean2D
	@(cd $(SRC)/AMG/obj2D; $(MAKE) clean )
	@(cd $(SRC)/Refinement/obj2D; $(MAKE) clean )
	@(cd $(SRC)/QuadFormulas/obj2D; $(MAKE) clean )
	@(cd $(SRC)/Geometry/obj2D; $(MAKE) clean )
	@(cd $(SRC)/General/obj2D; $(MAKE) clean )
	@(cd $(SRC)/FE/obj2D; $(MAKE) clean )
	@(cd $(SRC)/System/obj2D; $(MAKE) clean )
	@(cd $(SRC)/PBE/obj2D; $(MAKE) clean )
	@($(RM) ii_files lib2D/lib_repository lib2D/SunWS_cache )


clean2D_PAR:
	@echo clean2D
	@(cd $(SRC)/AMG/obj2D_par; $(MAKE) clean )
	@(cd $(SRC)/Refinement/obj2D_par; $(MAKE) clean )
	@(cd $(SRC)/QuadFormulas/obj2D_par; $(MAKE) clean )
	@(cd $(SRC)/Geometry/obj2D_par; $(MAKE) clean )
	@(cd $(SRC)/General/obj2D_par; $(MAKE) clean )
	@(cd $(SRC)/FE/obj2D_par; $(MAKE) clean )
	@(cd $(SRC)/System/obj2D_par; $(MAKE) clean )
	@(cd $(SRC)/Parallel/obj2D_par; $(MAKE) clean )
	@($(RM) ii_files lib2D/lib_repository lib2D/SunWS_cache )

clean3D:
	@echo clean3D
	@(cd $(SRC)/AMG/obj3D; $(MAKE) clean )
	@(cd $(SRC)/Refinement/obj3D; $(MAKE) clean )
	@(cd $(SRC)/QuadFormulas/obj3D; $(MAKE) clean )
	@(cd $(SRC)/Geometry/obj3D; $(MAKE) clean )
	@(cd $(SRC)/General/obj3D; $(MAKE) clean )
	@(cd $(SRC)/FE/obj3D; $(MAKE) clean )
	@(cd $(SRC)/System/obj3D; $(MAKE) clean )
	@(cd $(SRC)/PBE/obj3D; $(MAKE) clean )
	@($(RM) ii_files lib3D/lib_repository lib3D/SunWS_cache )

clean3D_PAR:
	@echo clean3D
	@(cd $(SRC)/AMG/obj3D_par; $(MAKE) clean )
	@(cd $(SRC)/Refinement/obj3D_par; $(MAKE) clean )
	@(cd $(SRC)/QuadFormulas/obj3D_par; $(MAKE) clean )
	@(cd $(SRC)/Geometry/obj3D_par; $(MAKE) clean )
	@(cd $(SRC)/General/obj3D_par; $(MAKE) clean )
	@(cd $(SRC)/FE/obj3D_par; $(MAKE) clean )
	@(cd $(SRC)/System/obj3D_par; $(MAKE) clean )
	@(cd $(SRC)/Parallel/obj3D_par; $(MAKE) clean )
	@($(RM) ii_files lib3D/lib_repository lib3D/SunWS_cache )
