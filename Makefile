.PHONY: all paper experiments check clean

all: paper experiments

paper:
	mkdir -p build
	latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error \
		-outdir=build paper/moulla_2026_stokes_pairing_gauge.tex

experiments:
	mkdir -p results
	python experiments/automatic_stokes_examples.py \
		--output results/automatic_stokes_examples.json

check: experiments
	python -m compileall -q experiments

clean:
	latexmk -C -outdir=build paper/moulla_2026_stokes_pairing_gauge.tex
	rm -rf build experiments/__pycache__

