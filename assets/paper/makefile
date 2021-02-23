all: manuscript/haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.pdf
submission: figures/Figure_1.png figures/Figure_2.png figures/Figure_3.png
submission: figures/Figure_4.png figures/Figure_5.png
submission: manuscript/haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.odt

manuscript/%.pdf: manuscript/%.tex
	pdflatex $<
	bibtex $(basename $(@F))
	pdflatex $<
	pdflatex $<
	rm -f $(basename $(@F)).aux $(basename $(@F)).bbl $(basename $(@F)).blg
	rm -f $(basename $(@F))-blx.bib $(basename $(@F)).log
	rm -f $(basename $(@F)).toc $(basename $(@F)).run.xml
	rm -f $(basename $(@F)).bcf
	mv $(@F) $@

figures/%.png: figures/%.pdf
	convert -density 300 $< $@

manuscript/%.odt: latex/%.tex
	./tex2office -i manuscript/references.bib figures tables --svg $<
	mv latex/$(@F) $@
