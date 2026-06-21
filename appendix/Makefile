# Compile the scoring appendix. Requires lualatex (OldStandard OpenType fonts; pinned in
# .latexmkrc) + bibtex, and rsvg-convert (librsvg) to turn the generated SVG assets into PDF.
# Assets (figures + alignments) are produced by ../bench/gen_appendix_assets.py.
FIGS = figures/region_corr.pdf figures/loo_prauc.pdf aln/GILGFVFTL_msa.pdf aln/NLVPMVATV_msa.pdf

vdjmatch_scoring.pdf: vdjmatch_scoring.tex refs.bib $(FIGS)
	latexmk -lualatex -interaction=nonstopmode -halt-on-error vdjmatch_scoring.tex

%.pdf: %.svg
	rsvg-convert -f pdf -o $@ $<

.PHONY: assets clean
assets:
	python ../bench/gen_appendix_assets.py
clean:
	latexmk -C vdjmatch_scoring.tex
