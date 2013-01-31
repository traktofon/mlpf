# vim: set ts=8 noexpandtab :

TGT := $(dd)/coul4test
DEP := $(dd)/Deps.mk

PROGS := $(PROGS) $(TGT)

SOURCES := \
   coul4mod.f90 \
   coul4test.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

MODS_$(dd) := \
   $(OBJDIR)/coul4mod.mod

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))

# build rule

$(TGT): $(OBJS_coul4) $(OBJS_core)
	$(FC) -o $@ $+ $(LIBS)

# module dependencies

$(DEP): $(MODS_core) $(MODS_coul4) $(SRC_coul4)
	$(FC) $(DEPFLAGS) -I$(OBJDIR) -J$(OBJDIR) $(SRC_test) | sed -e "s@^\(\S\)@$(OBJDIR)/\1@" > $@

-include $(DEP)

