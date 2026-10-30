# Build the slides (pdflatex, twice for the progress bar and references).
# Numbers and plots come from ../report/data (see ../report/experiments).
MAIN     = slides
PDFLATEX ?= pdflatex
.PHONY: all clean
all: $(MAIN).pdf
$(MAIN).pdf: $(MAIN).tex ../report/data/results.tex ../report/data/*.dat ../report/figures/tikz/*.tex
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MAIN) >/dev/null
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MAIN) | tail -2
clean:
	rm -f $(MAIN).aux $(MAIN).log $(MAIN).nav $(MAIN).out $(MAIN).snm $(MAIN).toc $(MAIN).vrb
