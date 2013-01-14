# vim: set ts=8 noexpandtab :

FC := gfortran
FFLAGS := -Wall -g -fbounds-check
FFLAGS := -O2
DEPFLAGS := -cpp -MM
LIBS := -llapack -lblas
