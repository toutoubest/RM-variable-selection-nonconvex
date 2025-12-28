# RM-variable-selection-nonconvex

This repository contains code for robust and missing-data-aware variable selection with nonconvex penalties in high-dimensional regression, as described in our paper.

## File overview

- **helper_funcs and synthetic .R**  
  Contains helper functions and the synthetic dataset implementation code.

- **Real datasets example code.R**  
  Example code for running our method on two real datasets.

## Getting started

- Requires R (>= 4.0) and the following packages: `ncvreg`, `glmnet`, `robustbase`, `data.table`, `ggplot2`, etc.
- See comments in the scripts for how to run synthetic and real data examples.

## Real datasets:
- Riboflavin dataset is in the R package hdi
- Breast cancer dataset is in our github breast_expression_original.csv 

This repository contains the code associated with the paper:

**Zhang, Linging.**  
*Introducing HYBRID and ENSEMBLE: novel nonconvex penalization strategies for robust variable selection under missing data*,  
**Computational Statistics**, Springer Nature, 2025/12/26.

The paper has been published online and can be accessed via the Springer Nature SharedIt link:  
https://rdcu.be/eWph2
DOI: https://doi.org/10.1007/s00180-025-01707-1
