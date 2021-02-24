all: haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.pdf
submission: ../figures/Figure_1.png ../figures/Figure_2.png ../figures/Figure_3.png ../figures/Figure_4.png ../figures/Figure_5.png
submission: haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.odt

%.pdf: %.tex
	pdflatex $<
	bibtex $(basename $(@F))
	pdflatex $<
	pdflatex $<
	rm -f $(basename $(@F)).aux $(basename $(@F)).bbl $(basename $(@F)).blg
	rm -f $(basename $(@F))-blx.bib $(basename $(@F)).log
	rm -f $(basename $(@F)).toc $(basename $(@F)).run.xml
	rm -f $(basename $(@F)).bcf

%.png: %.pdf
	convert -density 300 $< $@

%.odt: %.tex
	cd .. && ./tools/tex2office -i manuscript/references.bib figures tables --svg --rem-dot-dot manuscript/$<

clean:
	rm -f haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.aux
	rm -f haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.bbl
	rm -f haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres-blx.bib
	rm -f haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.blg
	rm -f haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.log
	rm -f haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.odt
	rm -f haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.pdf
	rm -f haplotype-diversity-and-sequence-heterogeneity-of-human-telomeres.run.xml
