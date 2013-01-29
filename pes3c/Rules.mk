# vim: set ts=8 noexpandtab :

TGT := $(dd)/pes3ctest
DEP := $(dd)/Deps.mk

PROGS := $(PROGS) $(TGT)

SOURCES90 := \
   pes3cmod.f90 \
   pes3ctest.f90
SOURCES77 := \
   pes3cvpd.f

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES77)) \
   $(addprefix $(dd)/,$(SOURCES90))

MODS_$(dd) := \
   $(OBJDIR)/pes3cmod.mod

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES77:.f=.o)) \
   $(addprefix $(OBJDIR)/,$(SOURCES90:.f90=.o))

# build rule

$(TGT): $(OBJS_pes3c) $(OBJS_core)
	$(FC) -o $@ $+ $(LIBS)

# module dependencies

$(DEP): $(MODS_core) $(MODS_pes3c) $(SRC_pes3c)
	$(FC) $(DEPFLAGS) -I$(OBJDIR) -J$(OBJDIR) $(SRC_pes3c) | sed -e "s@^\(\S\)@$(OBJDIR)/\1@" > $@

-include $(DEP)

