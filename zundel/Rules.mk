# vim: set ts=8 noexpandtab :

TGT := $(dd)/zundeltest
DEP := $(dd)/Deps.mk

PROGS := $(PROGS) $(TGT)

SOURCES := \
   zundel_m.f90 \
   zundeltest.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

MODS_$(dd) := \
   $(OBJDIR)/zundel_m.mod

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o)) \
   $(OBJDIR)/h5o2.pes.o

# build rules

$(TGT): $(OBJS_zundel) $(OBJS_COMMON)
	$(FC) -o $@ $+ $(LIBS)

$(OBJDIR)/h5o2.pes.o: h5o2.pes.f90
	$(FC) $(FFLAGS) -ffixed-form -c -o $@ $<

# module dependencies

$(DEP): $(MODS_COMMON) $(MODS_zundel) $(SRC_zundel)
	$(FC) $(DEPFLAGS) -I$(OBJDIR) -J$(OBJDIR) $(SRC_zundel) | sed -e "s@^\(\S\)@$(OBJDIR)/\1@" > $@

-include $(DEP)

