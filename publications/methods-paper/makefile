all: renders/haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.pdf
all: renders/figures/Figure-1.pdf renders/figures/Figure-2.pdf
all: renders/figures/Figure-3.pdf renders/figures/Figure-4.pdf
all: renders/figures/Figure-S1.pdf renders/figures/Figure-S2.pdf
all: renders/figures/Figure-S3.pdf renders/figures/Figure-S4.pdf
all: renders/figures/Figure-S5.pdf renders/figures/Figure-S6.pdf
all: renders/figures/Figure-S7.pdf renders/figures/Figure-SH.pdf

submission: renders/figures/Figure-1.png renders/figures/Figure-2.png
submission: renders/figures/Figure-3.png renders/figures/Figure-4.png
submission: submission/figures/Figure_1.pdf submission/figures/Figure_2.pdf
submission: submission/figures/Figure_3.pdf submission/figures/Figure_4.pdf
submission: submission/figures/Figure_5.pdf
submission: submission/figures/Supplemental_Fig_S1.pdf
submission: submission/figures/Supplemental_Fig_S2.pdf
submission: submission/figures/Supplemental_Fig_S3.pdf
submission: submission/figures/Supplemental_Fig_S4.pdf
submission: submission/figures/Supplemental_Fig_S5.pdf
submission: submission/figures/Supplemental_Fig_S6.pdf
submission: submission/figures/Supplemental_Fig_S7.pdf
submission: renders/haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.odt

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

renders/figures/%.png: renders/figures/%.pdf
	convert -density 300 $< $@

submission/figures/Figure_%.pdf: renders/figures/Figure-%.pdf
	cp $< $@

submission/figures/Supplemental_%.pdf: latex/supplemental-figures-GR/Supplemental_%.tex
	pdflatex $<
	rm -f $(basename $(@F)).aux $(basename $(@F)).log
	mv $(@F) $@

renders/%.odt: latex/%.tex
	./tex2office -i latex/references.bib latex/figures latex/tables -- $<
	mv latex/$(@F) $@
