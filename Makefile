# vim: set ts=8 noexpandtab :

FC := gfortran
FFLAGS := -Wall -g -fbounds-check
DEPFLAGS := -cpp -MM
LIBS :=

SOURCES := \
   base.f90 dof.f90 genpot.f90 graphviz.f90 linear.f90 modeutil.f90 test.f90 tree.f90 tuckerdecomp.f90
MODULES := \
   base.mod dof.mod genpot.mod graphviz.mod linear.mod modeutil.mod tree.mod tuckerdecomp.mod

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

allmod: $(MODULES)

dep: deps.mk

deps.mk: $(MODULES)
	$(FC) $(DEPFLAGS) $(SOURCES) > deps.mk

include deps.mk

# others

clean:
	rm -f *.o *.mod
