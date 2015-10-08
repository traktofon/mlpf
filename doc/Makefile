# vim: set ts=8 noexpandtab :

MDFILES = $(wildcard *.md)
HTMLFILES = $(patsubst %.md, %.html, $(MDFILES))

DOTFILES = $(wildcard *.dot)
PNGFILES = $(patsubst %.dot, %.png, $(DOTFILES))

.PHONY: all clean html graphs

all: html graphs

html: $(HTMLFILES)

graphs: $(PNGFILES)

%.html: %.md
	cat $< | sed -e 's/\.md)/.html)/g' | pandoc -f markdown -t html -c common.css -s -o $@

%.png: %.dot
	dot -Tpng -o$@ $<

clean:
	rm -f $(HTMLFILES) $(PNGFILES)
