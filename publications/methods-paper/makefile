all: heterogeneity-of-telomeres-revealed-by-long-read-sequencing.pdf

%.pdf: %.tex
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	pdflatex $<
	rm -f $(basename $<).aux $(basename $<).bbl $(basename $<).blg
	rm -f $(basename $<)-blx.bib $(basename $<).log $(basename $<).run.xml
	rm -f $(basename $<).toc

%.odt: %.tex
	./tex2office -i references.bib figures tables --svg heterogeneity-of-telomeres-revealed-by-long-read-sequencing.tex
