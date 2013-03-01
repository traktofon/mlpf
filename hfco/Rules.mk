# vim: set ts=8 noexpandtab :

TGT := $(dd)/hfcotest
DEP := $(dd)/Deps.mk

PROGS := $(PROGS) $(TGT)

SOURCES90 := \
   hfco_m.f90 \
   hfcotest.f90
SOURCES77 := \
   hfco.f

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES77)) \
   $(addprefix $(dd)/,$(SOURCES90))

MODS_$(dd) := \
   $(OBJDIR)/hfco_m.mod

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES77:.f=.o)) \
   $(addprefix $(OBJDIR)/,$(SOURCES90:.f90=.o))

# build rule

$(TGT): $(OBJS_hfco) $(OBJS_COMMON)
	$(FC) -o $@ $+ $(LIBS)

# module dependencies

$(DEP): $(MODS_COMMON) $(MODS_hfco) $(SRC_hfco)
	$(FC) $(DEPFLAGS) -I$(OBJDIR) -J$(TMPDIR) $(SRC_hfco) | sed -e "s@^\(\S\)@$(OBJDIR)/\1@" -e "s@$(TMPDIR)/@$(OBJDIR)/@" > $@

-include $(DEP)

