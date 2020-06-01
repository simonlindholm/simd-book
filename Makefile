MAKEFLAGS += --no-builtin-rules
LATEXCMD = pdflatex -output-directory build

.PHONY: fast bib clean

fast: | build
	$(LATEXCMD) book.tex </dev/null
	cp build/book.pdf book.pdf

bib: | build
	bibtex build/book

clean:
	$(RM) -r build/

build:
	mkdir -p build/

