name = haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres

pdf: $(name).pdf submission/$(name)-SUPPLEMENTAL.pdf
all: $(name).pdf submission/$(name)-SUPPLEMENTAL.pdf submission/$(name).odt

%.pdf: %.tex
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	pdflatex $<
	rm -f $(basename $<).aux $(basename $<).bbl $(basename $<).blg
	rm -f $(basename $<)-blx.bib $(basename $<).log $(basename $<).run.xml
	rm -f $(basename $<).toc $(basename $<).bcf

submission/%.pdf: %.tex
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	pdflatex $<
	rm -f $(basename $<).aux $(basename $<).bbl $(basename $<).blg
	rm -f $(basename $<)-blx.bib $(basename $<).log $(basename $<).run.xml
	rm -f $(basename $<).toc $(basename $<).bcf
	mkdir -p submission
	mv $(basename $<).pdf $@

submission/%.odt: %.tex
	mkdir -p submission
	./tex2office -i references.bib figures tables --svg $<
	mv $(basename $<).odt $@
