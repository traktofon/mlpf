# vim: set ts=8 noexpandtab :

FC := gfortran
FFLAGS := -Wall -g -fbounds-check
LIBS :=

# program targets

all: test

TEST_OBJS = test.o graphviz.o tree.o dof.o base.o
test: $(TEST_OBJS)
	$(FC) -o $@ $(TEST_OBJS) $(LIBS)

# build rules

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

%.mod: %.o
	@true

# module dependencies

test.o: \
    base.mod dof.mod tree.mod

tuckerdecomp.o: \
    base.mod linear.mod modeutil.mod

linear.o: \
    base.mod

tree.o: \
    base.mod dof.mod

dof.o: \
    base.mod

graphviz.o: \
    dof.mod tree.mod

# others

clean:
	rm -f *.o *.mod
