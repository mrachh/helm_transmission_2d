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

FFLAGS= -fPIC -O3 -march=native -funroll-loops -std=legacy -w 
CFLAGS= -std=c99 
CFLAGS+= $(FFLAGS)
#LDFLAGS="-L/usr/local/opt/openblas/lib"
CPPFLAGS="-I/usr/local/opt/openblas/include" 

CLINK = -lgfortran -lm -ldl

LIBS = -lm

# flags for MATLAB MEX compilation..
MFLAGS=-largeArrayDims -lgfortran -DMWF77_UNDERSCORE1 -lm  -ldl 
MWFLAGS=-c99complex -i8 

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=../../../mwrap/mwrap


# For your OS, override the above by placing make variables in make.inc
-include make.inc


# objects to compile
#
# Common objects
 
# Helmholtz objects
OBJS = src/hank103.o src/prini.o src/helm_kernels.o  \
	src/formsysmatbac.o src/kern_mats.o src/lap_kernels.o \
	src/durdec.o src/corrand4.o src/dumb_conres.o \
	src/legeexps.o src/curve_filtering2.o src/curve_resampler.o \
	src/dfft.o src/cdjseval2d.o src/h2dcommon.o

TOBJS = src/hkrand.o src/dlaran.o

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

LIBNAME = libhelmtrans
STATICLIB = lib-static/$(LIBNAME).a

$(STATICLIB): $(OBJS) 
	ar rcs $(STATICLIB) $(OBJS)

# matlab..
MWRAPFILE = kern_mats
MWRAPFILE2 = helm_kernels
MWRAPFILE3 = lap_kernels
MWRAPFILE4 = curve_resampler
GATEWAY = $(MWRAPFILE)
GATEWAY2 = $(MWRAPFILE2)
GATEWAY3 = $(MWRAPFILE3)
GATEWAY4 = $(MWRAPFILE4)

matlab:	$(STATICLIB) matlab/src/$(GATEWAY).c matlab/src/$(GATEWAY2).c matlab/src/$(GATEWAY3).c matlab/src/$(GATEWAY4).c
	$(MEX) matlab/src/$(GATEWAY).c $(STATICLIB) $(MFLAGS) -output matlab/src/kern_mats $(MEX_LIBS);
	$(MEX) matlab/src/$(GATEWAY2).c $(STATICLIB) $(MFLAGS) -output matlab/src/helm_kernels $(MEXLIBS);
	$(MEX) matlab/src/$(GATEWAY3).c $(STATICLIB) $(MFLAGS) -output matlab/src/lap_kernels $(MEXLIBS);
	$(MEX) matlab/src/$(GATEWAY4).c $(STATICLIB) $(MFLAGS) -output matlab/src/curve_resampler $(MEXLIBS);


mex:  $(STATICLIB)
	cd matlab; cd src;  $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAPFILE).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw;\
	$(MEX) $(GATEWAY).c ../../$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE) $(MEX_LIBS); \
	$(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY2) -mb $(MWRAPFILE2).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY2) -c $(GATEWAY2).c $(MWRAPFILE2).mw;\
	$(MEX) $(GATEWAY2).c ../../$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE2) $(MEX_LIBS);\
	$(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY3) -mb $(MWRAPFILE3).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY3) -c $(GATEWAY3).c $(MWRAPFILE3).mw;\
	$(MEX) $(GATEWAY3).c ../../$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE3) $(MEX_LIBS); \
	$(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY4) -mb $(MWRAPFILE4).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY4) -c $(GATEWAY4).c $(MWRAPFILE4).mw;\
	$(MEX) $(GATEWAY4).c ../../$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE4) $(MEX_LIBS);

#
##  examples
#

#examples: $(OBJS) $(TOBJS) examples/ext_dir examples/trans examples/curve examples/circ
examples: $(OBJS) $(TOBJS) examples/curve
#	time -p ./examples/int2-dir
#	time -p ./examples/int2-trans
	time -p ./examples/int2-curve
#	time -p ./examples/int2-circ

examples/ext_dir:
	$(FC) $(FFLAGS) examples/ext_dir_solver.f $(OBJS) -o examples/int2-dir -lopenblas $(LDFLAGS) 

examples/trans:
	$(FC) $(FFLAGS) examples/trans_solver.f $(OBJS) -o examples/int2-trans -lopenblas $(LDFLAGS)

examples/curve:
	$(FC) $(FFLAGS) examples/test_funcurv_fft.f $(OBJS) $(TOBJS) -o examples/int2-curve 

examples/circ:
	$(FC) $(FFLAGS) examples/circ_test2.f $(OBJS) $(TOBJS) -o examples/int2-circ  -lopenblas $(LDFLAGS)

clean: objclean
	rm -f examples/ext_dir_solver examples/trans_solver

objclean: 
	rm -f $(OBJS) 
	rm -f test/*.o examples/*.o c/*.o
