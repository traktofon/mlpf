# vim: set ts=8 noexpandtab :

include local.mk

SOURCES := \
   base.f90 \
   dof.f90 \
   genpot.f90 \
   graphviz.f90 \
   linear.f90 \
   modeutil.f90 \
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

all: test test_pes3c

LIB_OBJS = \
    hiertuck.o \
    tuckerdecomp.o \
    modeutil.o \
    linear.o \
    graphviz.o \
    tree.o \
    dof.o \
    genpot.o \
    testfunc.o \
    pes3cvpd.o \
    base.o

TEST_OBJS = test.o $(LIB_OBJS)

test: $(TEST_OBJS)
	$(FC) -o $@ $(TEST_OBJS) $(LIBS)

TEST_PES3C_OBJS = test_pes3c.o $(LIB_OBJS)

test_pes3c: $(TEST_PES3C_OBJS)
	$(FC) -o $@ $(TEST_PES3C_OBJS) $(LIBS)

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
