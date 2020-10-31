all: renders/haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.pdf
all: renders/figures/Figure-1.pdf renders/figures/Figure-2.pdf
all: renders/figures/Figure-3.pdf renders/figures/Figure-4.pdf
all: renders/figures/Figure-S1.pdf renders/figures/Figure-S2.pdf
all: renders/figures/Figure-S3.pdf renders/figures/Figure-S4.pdf
all: renders/figures/Figure-S5.pdf renders/figures/Figure-S6.pdf
all: renders/figures/Figure-S7.pdf renders/figures/Figure-SH.pdf
all: submission/figures/Supplemental_Fig_S1.pdf
all: submission/figures/Supplemental_Fig_S2.pdf
all: submission/figures/Supplemental_Fig_S3.pdf
all: submission/figures/Supplemental_Fig_S4.pdf
all: submission/figures/Supplemental_Fig_S5.pdf
all: submission/figures/Supplemental_Fig_S6.pdf
all: submission/figures/Supplemental_Fig_S7.pdf
all: submission/figures/Supplemental_Fig_S8.pdf

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

submission/figures/%.pdf: latex/supplemental-figures-GR/%.tex
	pdflatex $<
	rm -f $(basename $(@F)).aux $(basename $(@F)).log
	mv $(@F) $@

renders/%.odt: latex/%.tex
	./tex2office -i latex/references.bib latex/figures latex/tables --svg $<
	mv latex/$(@F) $<
