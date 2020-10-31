all: renders/haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.pdf
all: renders/figures/Figure-1.pdf renders/figures/Figure-2.pdf
all: renders/figures/Figure-3.pdf renders/figures/Figure-4.pdf
all: renders/figures/Figure-S1.pdf renders/figures/Figure-S2.pdf
all: renders/figures/Figure-S3.pdf renders/figures/Figure-S4.pdf
all: renders/figures/Figure-S5.pdf renders/figures/Figure-S6.pdf
all: renders/figures/Figure-S7.pdf renders/figures/Figure-SH.pdf

renders/%.pdf: latex/%.tex
	pdflatex $<
	bibtex $(basename $(@F))
	pdflatex $<
	pdflatex $<
	rm -f $(basename $(@F)).aux $(basename $(@F)).bbl $(basename $(@F)).blg
	rm -f $(basename $(@F))-blx.bib $(basename $(@F)).log
	rm -f $(basename $(@F)).toc $(basename $(@F)).run.xml
	mv $(@F) $@

renders/figures/%.pdf: latex/figures/%.tex
	pdflatex $<
	rm -f $(basename $(@F)).aux $(basename $(@F)).log
	mv $(@F) $@

renders/%.odt: latex/%.tex
	./tex2office -i latex/references.bib latex/figures latex/tables --svg $<
	mv latex/$(@F) $<
