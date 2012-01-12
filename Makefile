#include make.inc

F90=gfortran
CC=gcc

RM = \rm -rf
CP = \cp
MV = \mv
LN = \ln -s

SOURCEDIR=./SOURCES
FFTW_INC=-I$(SOURCEDIR)/fftw/mac/include
FFTW_LIB=-L$(SOURCEDIR)/fftw/mac/lib -lfftw3 -lm 
GSL_LIB=-L/opt/local/lib -lgsl -lgslcblas
GSL_INC=-I/opt/local/include
FGSL_INC=-I$(SOURCEDIR)/fgsl/mac/include/$(F90) $(GSL_LIB)
FGSL_LIB=-L$(SOURCEDIR)/fgsl/mac/lib -lfgsl_$(F90)

FPP=-cpp

EXECUTABLE=spectral_analysis

SOURCES_OBJS=ONERROR.o\
			 DATATYPES.o \
			 TIMER.o\
			 INIT.o\
			 IO.o\
			 ALGEBRA.o\
			 STATISTICS.o\
			 PREMIXED.o\
			 SPECTRAL_ANALYSIS.o
#
vpath %.f90 $(SOURCEDIR)
.PHONY: main clean
.INTERMEDIATE: $(SOURCES_OBJS)
#
main : $(SOURCES_OBJS) 
	@echo
	@echo "Linking Object Files   "
	@$(F90) $(OPT) $(FPP) $(FGSL_INC) $(FFTW_INC) -o $@ $(SOURCES_OBJS) $(FFTW_LIB) $(FGSL_LIB) -lm
	@-touch $(EXECUTABLE)
	@chmod a+x $(EXECUTABLE)
	@echo

%.o:%.f90
	@echo
	@echo "****************************************"
	@echo "Compiling     " $<
	@$(F90) $(OPT) $(FPP) -c $(FFTW_INC) $(FGSL_INC) $< $(FFTW_LIB) $(FGSL_LIB) -lm
	@echo "****************************************"
	@echo

#==================================================
#
#  General clean up
#
clean :
	@echo
	@echo "****************************************"
	@echo
	@echo "Removing object files"
	@echo
	@$(RM) $(SOURCES_OBJS)  $(OUTDIR)/$(EXECUTABLE) 
	@echo
	@echo "Removing included files"
	@echo
	@$(RM) *.inc  *.h  *.c *.F *.f *.f90 *.o  *.mod  *.cmn *_LINK*
	@echo "****************************************"
	@echo
