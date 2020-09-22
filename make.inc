# makefile overrides to use Intel ICC compiler, double prec only.

CC=gcc
FC=gfortran

CFLAGS += -I src 

# OpenMP with gcc on OSX needs the following
OMPFLAGS = -fopenmp
OMPLIBS = -lgomp 

#MATLAB interface:
MFLAGS += -L/usr/local/Cellar/gcc/10.2.0/lib/gcc/10
MEX = $(shell ls -d /Applications/MATLAB_R201*.app)/bin/mex


