# vim: set ts=8 noexpandtab :

all: targets

OBJDIR := obj
TMPDIR := scr
include local.mk

PROGS :=
ALLSOURCES :=

dd := core
include $(dd)/Rules.mk
dd := test
include $(dd)/Rules.mk
dd := pes3c
include $(dd)/Rules.mk
dd := zundel
include $(dd)/Rules.mk
dd := hfco
include $(dd)/Rules.mk
dd := coul4
include $(dd)/Rules.mk

# where to look for source files

vpath %.f90 core
vpath %.f90 test
vpath %.f90 pes3c
vpath %.f   pes3c
vpath %.f90 zundel
vpath %.f90 hfco
vpath %.f   hfco
vpath %.f90 coul4

# general build rules

$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/%.mod: $(OBJDIR)/%.o
	@true

$(OBJDIR)/%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

# other rules

targets: $(PROGS)

clean:
	rm -f $(OBJDIR)/* $(PROGS)

# dependencies

-include deps.mk

dep:
	bin/mkdep $(OBJDIR) $(ALLSOURCES) > deps.mk

