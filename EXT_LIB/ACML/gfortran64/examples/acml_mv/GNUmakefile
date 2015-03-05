
# GNUmakefile to build and execute ACML example programs

# Fortran and C compilers and flags (modify if necessary)
#
# Location of ACML installation
ACMLDIR := /homes/sashikum/ParMooN/EXT_LIB/ACML/gfortran64

F77 := gfortran
FLINK := $(F77)
FFLAGS :=  
FLINKFLAGS :=  
FLINKLIBS := -lrt
CC := gcc
CFLAGS := -I$(ACMLDIR)/include 
CLINK := $(CC)
CLINKFLAGS :=  
CLINKLIBS := -lgfortran -lm -lrt
CPLUSPLUS := g++
CPLUSPLUSFLAGS := -I$(ACMLDIR)/include 
CPLUSPLUSLINK := $(CPLUSPLUS)
CPLUSPLUSLINKFLAGS :=  
CPLUSPLUSLINKLIBS := -lgfortran -lrt -lm

obj := o 

.SUFFIXES: .f .c .cpp .$(obj) .exe .res .diff .expected

LIBACML_MV := $(ACMLDIR)/lib/libacml_mv.a



EXPORT_PATH := export LD_LIBRARY_PATH='$(ACMLDIR)/lib:$(LD_LIBRARY_PATH)'

# Source code of Fortran example programs
FSOURCE := testfortranif.f
FEXE := $(FSOURCE:.f=.exe)

# Source code of C example programs
CSOURCE := testCif.c
CEXE := $(CSOURCE:.c=.exe)

# Source code of C++ example programs
CPLUSPLUSSOURCE := testCPPif.cpp
CPLUSPLUSEXE := $(CPLUSPLUSSOURCE:.cpp=.exe)

RESULTS := $(FSOURCE:.f=.res) $(CSOURCE:.c=.res) $(CPLUSPLUSSOURCE:.cpp=.res)

all: $(RESULTS) diffs

# Rule to compile, link and run a Fortran file
.f.res:
	@echo ""
	@echo "Compiling program $*.f:"
	$(F77) -c $(FFLAGS) $*.f -o $*.$(obj)
	@echo "Linking program $*.exe:"
	$(FLINK) $(FLINKFLAGS) $*.$(obj) $(LIBACML_MV) $(FLINKLIBS) -o $*.exe
	@echo "Running program $*.exe:"
	($(EXPORT_PATH); ./$*.exe > $@ 2>&1)
#       Show the results
	@cat $@
#	Compare the results with the expected results
	@diff -w $@ $*.expected 2>&1 | cat > $*.diff
#       Delete the diff file if it is empty
	@if [ -s $*.diff ]; then \
	  echo "!!! $* has diffs !!!"; \
	else \
	  rm -f $*.diff; \
	fi

# Rule to compile, link and run a C file
.c.res:
	@echo ""
	@echo "Compiling program $*.c:"
	$(CC) -c $(CFLAGS) $*.c -o $*.$(obj)
	@echo "Linking program $*.exe:"
	$(CLINK) $(CLINKFLAGS) $*.$(obj) $(LIBACML_MV) $(CLINKLIBS) -o $*.exe
	@echo "Running program $*.exe:"
	($(EXPORT_PATH); ./$*.exe > $@ 2>&1)
#       Show the results
	@cat $@
#	Compare the results with the expected results
	@diff -w $@ $*.expected 2>&1 | cat > $*.diff
#       Delete the diff file if it is empty
	@if [ -s $*.diff ]; then \
	  echo "!!! $* has diffs !!!"; \
	else \
	  rm -f $*.diff; \
	fi

# Rule to compile, link and run a C++ file
.cpp.res:
	@echo ""
	@echo "Compiling program $*.cpp:"
	$(CPLUSPLUS) -c $(CPLUSPLUSFLAGS) $*.cpp -o $*.$(obj)
	@echo "Linking program $*.exe:"
	$(CPLUSPLUSLINK) $(CPLUSPLUSLINKFLAGS) $*.$(obj) $(LIBACML_MV) $(CPLUSPLUSLINKLIBS) -o $*.exe
	@echo "Running program $*.exe:"
	($(EXPORT_PATH); ./$*.exe > $@ 2>&1)
#       Show the results
	@cat $@
#	Compare the results with the expected results
	@diff -w $@ $*.expected 2>&1 | cat > $*.diff
#       Delete the diff file if it is empty
	@if [ -s $*.diff ]; then \
	  echo "!!! $* has diffs !!!"; \
	else \
	  rm -f $*.diff; \
	fi

diffs: $(RESULTS)
	@(ndiff=`ls *.diff 2> /dev/null | wc -w | sed -e 's/ //g'` ; \
	  echo ""; \
	  if [ "$${ndiff}" = "0" ]; then \
	    echo "Testing: no example difference files were generated."; \
	    echo "Test passed OK"; \
	  else \
	    echo "Testing: $${ndiff} example difference files were generated:"; \
	    ls -l *.diff; \
	  fi)

clean:
	rm -f $(FEXE) $(CEXE) $(CPLUSPLUSEXE) \
	      $(FSOURCE:.f=.$(obj)) $(CSOURCE:.c=.$(obj)) \
	      $(CPLUSPLUSSOURCE:.cpp=.$(obj)) \
	      $(RESULTS) *.diff
