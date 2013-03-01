# vim: set ts=8 noexpandtab :

TGT := $(dd)/coul4test
DEP := $(dd)/Deps.mk

PROGS := $(PROGS) $(TGT)

SOURCES := \
   coul4_m.f90 \
   coul4test.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

MODS_$(dd) := \
   $(OBJDIR)/coul4_m.mod

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))

# build rule

$(TGT): $(OBJS_coul4) $(OBJS_COMMON)
	$(FC) -o $@ $+ $(LIBS)

# module dependencies

$(DEP): $(MODS_COMMON) $(MODS_coul4) $(SRC_coul4)
	$(FC) $(DEPFLAGS) -I$(OBJDIR) -J$(TMPDIR) $(SRC_coul4) | sed -e "s@^\(\S\)@$(OBJDIR)/\1@" -e "s@$(TMPDIR)/@$(OBJDIR)/@" > $@

-include $(DEP)

