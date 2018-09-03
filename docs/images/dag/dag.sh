#!/bin/sh

module load Graphviz/2.38.0-foss-2016b

python dag.py dag.dot dag.csv | dot -Tpdf -Grankdir="LR" > dag.pdf
python dag.py dag.dot dag.csv | dot -Tpng -Grankdir="LR" > dag.png

