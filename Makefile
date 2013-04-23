# vim: set ts=8 noexpandtab :

all: targets

OBJDIR := obj

include local.mk

LINK := $(FC) $(LDFLAGS)

PROGS :=
ALLSOURCES :=

dd := core
include $(dd)/Rules.mk
dd := main
include $(dd)/Rules.mk

# general build rules

$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -I$(OBJDIR) $(MODFLAG)$(OBJDIR) -c -o $@ $<

$(OBJDIR)/%.mod: $(OBJDIR)/%.o
	@true

$(OBJDIR)/%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

# other rules

core : $(OBJS_core)

targets: $(PROGS)

clean:
	rm -f $(OBJDIR)/* $(PROGS) .hgstamp

# dependencies

-include Deps.mk

dep:
	bin/mkdep $(OBJDIR) $(ALLSOURCES) > Deps.mk

# ctags

tags: $(ALLSOURCES)
	ctags --fortran-kinds=+i $(ALLSOURCES)
