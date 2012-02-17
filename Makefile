OS=$(shell uname)

CC=$(CC_$(OS))
F90=$(F90_$(OS))

CC_Linux=icc
CC_Darwin=gcc
F90_Linux=ifort
F90_Darwin=gfortran

SOURCEDIR=./SOURCES
FFTW_INC=-I$(SOURCEDIR)/fftw/$(OS)/include
FFTW_LIB=-L$(SOURCEDIR)/fftw/$(OS)/lib -lfftw3 -lm 
GSL_LIB=-L$(SOURCEDIR)/gsl/lib -lgsl -lgslcblas
GSL_INC=-I$(SOURCEDIR)/gsl/include
FGSL_INC=-I$(SOURCEDIR)/fgsl/$(OS)/include/$(F90_Linux) $(GSL_LIB)
FGSL_LIB=-L$(SOURCEDIR)/fgsl/$(OS)/lib -lfgsl_$(F90_Linux)

ifeq ($(F90_$(OS)),"gfortran")
	FPP=-cpp
	OPT=-O3 -Wall
else
	FPP=-fpp
	OPT=-O3 -Warn all
endif
RM = \rm -rf

SOURCES=ONERROR NRTYPE TIMER_CLASS INIT IO_CLASS STATISTICS PREMIXED_CLASS  SPECTRAL_ANALYSIS
#ALGEBRA 
SOURCES_OBJS=$(addsuffix .o,$(SOURCES))
#
vpath %.f90 $(SOURCEDIR)
.PHONY: main clean
.INTERMEDIATE: $(SOURCES_OBJS) $(addsuffix .mod,$(SOURCES))
#
main : $(SOURCES_OBJS) 
	$(F90) $(OPT) $(FPP) $(FGSL_INC) $(FFTW_INC) -o $@ $(SOURCES_OBJS) $(FFTW_LIB) $(FGSL_LIB) -lm

%.o %.mod:%.f90
	$(F90) $(OPT) $(FPP) -c $(FFTW_INC) $(FGSL_INC) $< $(FFTW_LIB) $(FGSL_LIB) -lm

clean :
	@$(RM) $(SOURCES_OBJS) main 
	@$(RM) *.inc  *.h  *.c *.F *.f *.f90 *.o  *.mod
