.PHONY: install-go install-python download-data demo download- plot-recovery plot-linking plot-uniqueness plot-rounding plot-all pgs-selection solving

install-go:
	./scripts/install-go.sh

install-python:
	./sciprts/install-python.sh

download-data:
	./scripts/download_processed_1000genomes.sh

demo:
	go run eval/run.go -e=demo

download-results:
	./scripts/download_results.sh

plot-recovery:
	python3 analysis/plot.py -e=recovery
plot-linking:
	python3 analysis/plot.py -e=linking
plot-uniqueness:
	python3 analysis/plot.py -e=uniqueness
plot-rounding:
	python3 analysis/plot.py -e=rounding
plot-all:
	python3 analysis/plots.py -e=all

pgs-selection:
	mkdir catalog/
	go run info/preprocess.go -e=pgs-selection

solving:
	go run eval/run.go -e=solve