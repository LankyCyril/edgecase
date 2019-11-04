all: grigorev-poster.pdf

%.aux: %.tex baposter.cls references.bib
	pdflatex $<

%.bbl: %.aux references.bib
	bibtex $(basename $<)

%.pdf: %.tex %.bbl
	rm -f $(basename $<)-blx.bib $(basename $<).blg $(basename $<).run.xml
	rm -f $(basename $<).aux $(basename $<).bcf
	rm -f $(basename $<).toc $(basename $<).out
	pdflatex $<
	pdflatex $<
	rm -f $(basename $<)-blx.bib $(basename $<).blg $(basename $<).run.xml
	rm -f $(basename $<).aux $(basename $<).bcf
	rm -f $(basename $<).toc $(basename $<).out
