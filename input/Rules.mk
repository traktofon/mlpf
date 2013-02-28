# vim: set ts=8 noexpandtab :

DEP := $(dd)/Deps.mk

SOURCES := \
   parse_pbasis_m.f90 \
   parse_tree_m.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

MODS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.mod))

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))


# module dependencies

$(DEP): $(MODS_core) $(MODS_parse) $(MODS_dof) $(MODS_input) $(SRC_input)
	$(FC) $(DEPFLAGS) -I$(OBJDIR) -J$(OBJDIR) $(SRC_input) | sed -e "s@^\(\S\)@$(OBJDIR)/\1@" > $@

-include $(DEP)

