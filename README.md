# Errors upon log-transformation of scRNA-seq data

This repository contains code and manuscript files for the paper **Overcoming systematic errors caused by log-transformation of normalized single-cell RNA sequencing data**
by [Lun (2018)](https://doi.org/10.1101/404962).

The `report` subdirectory contains LaTeX files for the report.
The `scripts` subdirectory contains R code to reproduce the analyses described in the report:

- `spurious_sim.R`, to examine the non-zero log-fold change induced by the log-transformation.
- `substructure.R`, to examine the introduction of structural artifacts due the log-transformation.
- `alternatives.R`, to examine the performance of other transformations.
- `check_ercc.R`, to examine the behaviour of the log-transformation on real data.
- `power_sim.R`, to examine the consequences on detection power with an increased pseudo-count. 

