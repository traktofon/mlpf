# vim: set ts=8 noexpandtab :

all: targets

OBJDIR := obj
BINDIR := bin

include local.mk

dir_guard = @mkdir -p $(@D)

LINK := $(FC) $(LDFLAGS)

PROGS :=
ALLSOURCES :=

dd := core
include $(dd)/Rules.mk
dd := main
include $(dd)/Rules.mk
dd := vinfo
include $(dd)/Rules.mk
dd := mlpf2npot
include $(dd)/Rules.mk

# general build rules

$(OBJDIR)/%.o: %.f90
	$(dir_guard)
	$(FC) $(FFLAGS) -I $(OBJDIR) $(MODFLAG) $(OBJDIR) -c -o $@ $<

$(OBJDIR)/%.mod: $(OBJDIR)/%.o
	@true

$(OBJDIR)/%.o: %.f
	$(dir_guard)
	$(FC) $(FFLAGS) -c -o $@ $<

$(OBJDIR)/%.o: %.c
	$(dir_guard)
	$(CC) $(CFLAGS) -c -o $@ $<

# other rules

core : $(OBJS_core)

targets: $(PROGS)

clean:
	rm -f $(OBJDIR)/* $(PROGS)

# dependencies

-include Deps.mk

dep: $(ALLSOURCES)
	tools/mkdep $(OBJDIR) $(ALLSOURCES) > Deps.mk

# ctags

tags: $(ALLSOURCES)
	ctags --fortran-kinds=+i $(ALLSOURCES)
