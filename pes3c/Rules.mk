# vim: set ts=8 noexpandtab :

TGT := $(dd)/pes3ctest

SOURCES90 := \
   pes3c_m.f90 \
   pes3ctest.f90
SOURCES77 := \
   pes3cvpd.f

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES77)) \
   $(addprefix $(dd)/,$(SOURCES90))

MODS_$(dd) := \
   $(OBJDIR)/pes3c_m.mod

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES77:.f=.o)) \
   $(addprefix $(OBJDIR)/,$(SOURCES90:.f90=.o))

PROGS := $(PROGS) $(TGT)
ALLSOURCES := $(ALLSOURCES) $(SRC_$(dd))

# build rule

$(TGT): $(OBJS_pes3c) $(OBJS_core)
	$(FC) -o $@ $+ $(LIBS)

