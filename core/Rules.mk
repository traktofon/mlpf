# vim: set ts=8 noexpandtab :

DEP := $(dd)/Deps.mk

SOURCES := \
   base_m.f90 \
   strutil_m.f90 \
   logging_m.f90 \
   tokenize_m.f90 \
   units_m.f90 \
   linear_m.f90 \
   modeutil_m.f90 \
   tree_m.f90 \
   tuckerdecomp_m.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

MODS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.mod))

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))


# module dependencies

$(DEP): $(MODS_core) $(SRC_core)
	$(FC) $(DEPFLAGS) -I$(OBJDIR) -J$(OBJDIR) $(SRC_core) | sed -e "s@^\(\S\)@$(OBJDIR)/\1@" > $@

-include $(DEP)

