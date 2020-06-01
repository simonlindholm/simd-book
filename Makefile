MAKEFLAGS += --no-builtin-rules
LATEXCMD = pdflatex -output-directory build

.PHONY: fast bib clean

fast: | build
	$(LATEXCMD) book.tex </dev/null
	cp build/book.pdf book.pdf

bib: | build
	bibtex build/book
	cd binomial-coefficients/ && bibtex build/freestanding

bincoef: | binomial-coefficients/build
	cd binomial-coefficients/ && pdflatex -output-directory build freestanding.tex </dev/null
	cp binomial-coefficients/build/freestanding.pdf bincoef.pdf

clean:
	$(RM) -r build/ binomial-coefficients/

build:
	mkdir -p $@

binomial-coefficients/build:
	mkdir -p $@

