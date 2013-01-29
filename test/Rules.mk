# vim: set ts=8 noexpandtab :

TGT := $(dd)/test
DEP := $(dd)/Deps.mk

PROGS := $(PROGS) $(TGT)

SOURCES := \
   testfunc.f90 \
   test.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

MODS_$(dd) := \
   $(OBJDIR)/testfunc.mod

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))

# build rule

$(TGT): $(OBJS_test) $(OBJS_core)
	$(FC) -o $@ $+ $(LIBS)

# module dependencies

$(DEP): $(MODS_core) $(MODS_test) $(SRC_test)
	$(FC) $(DEPFLAGS) -I$(OBJDIR) -J$(OBJDIR) $(SRC_test) | sed -e "s@^\(\S\)@$(OBJDIR)/\1@" > $@

-include $(DEP)

