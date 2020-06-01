MAKEFLAGS += --no-builtin-rules

.PHONY: fast bib clean

fast:
	pdflatex book.tex

bib:
	bibtex book

clean:
	$(RM) book.{aux,bbl,blg,log,out,pdf,toc}

