# vim: set ts=8 noexpandtab :

DEP := $(dd)/Deps.mk

SOURCES := \
   tokenize_m.f90 \
   units_m.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

MODS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.mod))

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))


# module dependencies

$(DEP): $(MODS_core) $(MODS_map) $(MODS_parse) $(SRC_parse)
	$(FC) $(DEPFLAGS) -I$(OBJDIR) -J$(TMPDIR) $(SRC_parse) | sed -e "s@^\(\S\)@$(OBJDIR)/\1@" -e "s@$(TMPDIR)/@$(OBJDIR)/@" > $@

-include $(DEP)

