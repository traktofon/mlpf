# vim: set ts=8 noexpandtab :

TGT := $(dd)/test

SOURCES := \
   testfunc_m.f90 \
   test.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

MODS_$(dd) := \
   $(OBJDIR)/testfunc_m.mod

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))

PROGS := $(PROGS) $(TGT)
ALLSOURCES := $(ALLSOURCES) $(SRC_$(dd))

# build rule

$(TGT): $(OBJS_test) $(OBJS_COMMON)
	$(FC) -o $@ $+ $(LIBS)

