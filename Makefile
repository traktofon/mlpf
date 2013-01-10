# vim: set ts=8 noexpandtab :

FC := gfortran
FFLAGS := -Wall -g -fbounds-check

# build rules

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

%.mod: %.o
	@true

# module dependencies

tuckerdecomp.o: \
    base.mod linear.mod modeutil.mod

linear.o: \
    base.mod

tree.o: \
    base.mod dof.mod

dof.o: \
    base.mod

# others

clean:
	rm -f *.o *.mod
