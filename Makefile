# vim: set ts=8 noexpandtab :

all: targets

PROGS :=
OBJDIR := obj
include local.mk

dd := core
include $(dd)/Rules.mk
dd := test
include $(dd)/Rules.mk
dd := pes3c
include $(dd)/Rules.mk

# where to look for source files

vpath %.f90 core
vpath %.f90 test
vpath %.f90 pes3c
vpath %.f   pes3c

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
