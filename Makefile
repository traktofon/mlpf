# vim: set ts=8 noexpandtab :

FC := gfortran
FFLAGS := -Wall -g -fbounds-check
DEPFLAGS := -cpp -MM
LIBS := -llapack -lblas

SOURCES := \
   base.f90 \
   dof.f90 \
   genpot.f90 \
   graphviz.f90 \
   linear.f90 \
   modeutil.f90 \
   test.f90 \
   tree.f90 \
   tuckerdecomp.f90 \
   hiertuck.f90 \
   testfunc.f90
MODULES := \
   base.mod \
   dof.mod \
   genpot.mod \
   testfunc.mod \
   linear.mod \
   modeutil.mod \
   tree.mod \
   graphviz.mod \
   tuckerdecomp.mod \
   hiertuck.mod

# program targets

all: test

TEST_OBJS = \
    test.o \
    hiertuck.o \
    tuckerdecomp.o \
    modeutil.o \
    linear.o \
    graphviz.o \
    tree.o \
    dof.o \
    genpot.o \
    testfunc.o \
    base.o

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

deps.mk: $(MODULES) $(SOURCES)
	$(FC) $(DEPFLAGS) $(SOURCES) > deps.mk

include deps.mk

# others

clean:
	rm -f *.o *.mod
