# vim: set ts=8 noexpandtab :

MDFILES = $(wildcard *.md)
HTMLFILES = $(patsubst %.md, %.html, $(MDFILES))

.PHONY: all clean

all: $(HTMLFILES)

%.html: %.md
	cat $< | sed -e 's/\.md)/.html)/g' | pandoc -f markdown -t html -c common.css -s -o $@

clean:
	rm -f $(HTMLFILES)
