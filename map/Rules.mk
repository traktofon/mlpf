# vim: set ts=8 noexpandtab :

DEP := $(dd)/Deps.mk

SOURCES := \
   map_str2dbl_m.f90 \
   map_str2int_m.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

MODS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.mod))

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))


# module dependencies

$(DEP): $(MODS_core) $(MODS_map) $(SRC_map)
	$(FC) $(DEPFLAGS) -I$(OBJDIR) -Jtmp $(SRC_map) | sed -e "s@^\(\S\)@$(OBJDIR)/\1@" -e "s@tmp/@obj/@" > $@

-include $(DEP)

