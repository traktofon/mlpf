# vim: set ts=8 noexpandtab :

DEP := $(dd)/Deps.mk

SOURCES := \
   genpot_m.f90 \
   graphviz_m.f90 \
   hiertuck_m.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

MODS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.mod))

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))


# module dependencies

$(DEP): $(MODS_core) $(MODS_dof) $(MODS_mlpf) $(SRC_mlpf)
	$(FC) $(DEPFLAGS) -I$(OBJDIR) -J$(OBJDIR) $(SRC_mlpf) | sed -e "s@^\(\S\)@$(OBJDIR)/\1@" > $@

-include $(DEP)

