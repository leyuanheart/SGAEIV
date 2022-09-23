# Subgroup analysis in linear model with measurement error

This is the official implementation of "Subgroup analysis in linear model with measurement error" in Canadian Journal of Statistics.

## `sgaeiv` package overview

- `iteration_procedures.R`: the proposed algorithm in this paper.

- `iteration_ma.R`: the algorithm in ''Exploration of Heterogeneous Treatment Effects via Concave Fusion (Ma et al., 2020)" .

- `UpdateDelta.R`: functions used in the ADMM algorithm to update the parameter $\boldsymbol{\delta}$, which can be applied to both *SCAD* and *MCP* penalty functions.

## Software and main packages

- R: 3.6.0

- MASS

- snow

- fossil

- rlecuyer

- parallel

- dummy

- car

## Reproducibility

### Simulations

Set the working directory of **R** to the directory of `.R` files in `simulation` folder and run `test_onedim_gau_final.R`, `test_intcpt_3groups`, and `acc_nogroups` to get `example1.RData`, `example2.RData`, and `example3.RData` respectively in `simulation/results` folder, which are used for figures and tables.

### Real data

For the consideration of privacy, we can not provide the original data of the Lifestyle Education for Activity and Nutrition (LEAN) study.

We have done some pre-processing on the LEAN data. If you have the access to the LEAN data, you can run `data_generation.R` in `real_data` folder and get the data version used in this paper, `nine_dat.RData`. Then you can run `real_data.R` to get `real_data.RData` for figures and tables.

Note that you still need to set the working directory to the directory of these `.R` files.

### Figures and Tables

Set the working directory of **R** to the directory of `figures_tables.R` and run this file to get all figures and tables in this paper and Appendix in the supplementary material.
