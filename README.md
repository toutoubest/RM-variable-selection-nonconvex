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
- Breast cancer dataset is in out github Breast_expression_original.csv 
