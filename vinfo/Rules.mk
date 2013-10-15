# vim: set ts=8 noexpandtab :

TGT := $(dd)/vinfo

vpath %.f90 $(dd)

SOURCES := \
   vinfo.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))

PROGS += $(TGT)
ALLSOURCES += $(SRC_$(dd))

# build rule

$(TGT): $(OBJS_vinfo) $(OBJS_core)
	$(LINK) -o $@ $+ $(LIBS)

