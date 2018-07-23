# Pseudo-count choice for log-transformed scRNA-seq data

This repository contains code to explore some of the consequences of pseudo-count choice in log-transformation of scRNA-seq data.
The `report` subdirectory contains LaTeX files for the report.
The `scripts` subdirectory contains R code to reproduce the analyses described in the report:

- `spurious_sim.R`, to examine the non-zero log-fold change induced by the log-transformation.
- `substructure.R`, to examine the introduction of structural artifacts due the log-transformation.
- `alternatives.R`, to examine the performance of other transformations.
- `check_ercc.R`, to examine the behaviour of the log-transformation on real data.
- `power_sim.R`, to examine the consequences on detection power with an increased pseudo-count. 

