# vim: set ts=8 noexpandtab :

FC := gfortran
FFLAGS := -Wall -g -fbounds-check

tuckerdecomp.o: \
    base.mod linear.mod modeutil.mod

linear.o: \
    base.mod

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

%.mod: %.o
	@true
