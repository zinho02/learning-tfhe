# Build the report.  "make data" regenerates data/ from data/raw/ first.
MAIN     = main
PDFLATEX ?= pdflatex
.PHONY: all data clean
all: $(MAIN).pdf
$(MAIN).pdf: $(MAIN).tex refs.bib figures/*.tex figures/tikz/*.tex data/results.tex data/*.dat \
             sbc-template.sty sbc.bst
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MAIN) >/dev/null
	bibtex $(MAIN) >/dev/null
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MAIN) >/dev/null
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MAIN) | tail -2
data:
	python3 experiments/summarize.py
clean:
	rm -f $(MAIN).aux $(MAIN).bbl $(MAIN).blg $(MAIN).log $(MAIN).out
