F90=gfortran
CC=gcc

RM = \rm -rf

SOURCEDIR=./SOURCES
FFTW_INC=-I$(SOURCEDIR)/fftw/mac/include
FFTW_LIB=-L$(SOURCEDIR)/fftw/mac/lib -lfftw3 -lm 
GSL_LIB=-L/opt/local/lib -lgsl -lgslcblas
GSL_INC=-I/opt/local/include
FGSL_INC=-I$(SOURCEDIR)/fgsl/mac/include/$(F90) $(GSL_LIB)
FGSL_LIB=-L$(SOURCEDIR)/fgsl/mac/lib -lfgsl_$(F90)

FPP=-cpp
#OPT=-g -O0
OPT=-O3 -Wall

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
