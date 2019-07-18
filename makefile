# Makefile for FMM3D
#
# This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create 
# the file make.inc, which overrides the defaults below (which are 
# for ubunutu linux/gcc system). 

# compiler, and linking from C, fortran

CC=gcc
FC=gfortran

FFLAGS= -fPIC -O3 -march=native   
CFLAGS= -std=c99 
CFLAGS+= $(FFLAGS) 

CLINK = -lgfortran -lm -ldl

LIBS = -lm

# flags for MATLAB MEX compilation..
MFLAGS=-largeArrayDims -lgfortran -DMWF77_UNDERSCORE1 -lm  -ldl
MWFLAGS=-c99complex 

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=../../mwrap-0.33/mwrap


# For your OS, override the above by placing make variables in make.inc
-include make.inc


# objects to compile
#
# Common objects
 
# Helmholtz objects
OBJS = src/hank103.o src/prini.o src/helm_kernels.o  \
	src/formsysmatbac.o 

.PHONY: usage examples matlab debug

default: usage 

all: examples matlab

usage:
	@echo "Makefile for FMM3D. Specify what to make:"
	@echo "  make examples - compile and run fortran examples in examples/"
	@echo "  make matlab - compile matlab interfaces"
	@echo "  make mex - generate matlab interfaces (for expert users only, requires mwrap)"
	@echo "  make clean - also remove lib, MEX, py, and demo executables"


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@


# matlab..
MWRAPFILE = fmm3d
MWRAPFILE2 = fmm3d_legacy
GATEWAY = $(MWRAPFILE)
GATEWAY2 = $(MWRAPFILE2)

matlab:	$(STATICLIB) matlab/$(GATEWAY).c matlab/$(GATEWAY2).c
	$(MEX) matlab/$(GATEWAY).c $(STATICLIB) $(MFLAGS) -output matlab/fmm3d $(MEX_LIBS);
	$(MEX) matlab/$(GATEWAY2).c $(STATICLIB) $(MFLAGS) -output matlab/fmm3d_legacy $(MEXLIBS);


mex:  $(STATICLIB)
	cd matlab; $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAPFILE).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw;\
	$(MEX) $(GATEWAY).c ../$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE) $(MEX_LIBS); \
	$(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY2) -mb $(MWRAPFILE2).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY2) -c $(GATEWAY2).c $(MWRAPFILE2).mw;\
	$(MEX) $(GATEWAY2).c ../$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE2) $(MEX_LIBS);

#
##  examples
#

examples: $(OBJS) examples/ext_dir 
	time -p ./examples/ext_dir_solver

examples/ext_dir:
	$(FC) $(FFLAGS) examples/ext_dir_solver.f $(OBJS) -o examples/ext_dir_solver -llapack -lblas

debug: $(OBJS) examples/hh 
	time -p ./examples/hank

examples/hh:
	$(FC) $(FFLAGS) examples/test_hank103.f $(OBJS) -o examples/hank 



clean: objclean
	rm -f examples/ext_dir

objclean: 
	rm -f $(OBJS) 
	rm -f test/*.o examples/*.o c/*.o
