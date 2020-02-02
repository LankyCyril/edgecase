all: heterogeneity-of-telomeres-revealed-by-long-read-sequencing.pdf heterogeneity-of-telomeres-revealed-by-long-read-sequencing.docx

%.pdf: %.tex
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	pdflatex $<
	rm -f $(basename $<).aux $(basename $<).bbl $(basename $<).blg
	rm -f $(basename $<)-blx.bib $(basename $<).log $(basename $<).run.xml
	rm -f $(basename $<).toc

%.docx: %.tex
	pandoc -s $< --bibliography=references.bib --csl=csl/springer-basic-brackets.csl -o $@
