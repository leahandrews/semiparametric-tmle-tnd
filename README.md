# semiparametric-tmle-tnd
Repository for "Targeted maximum likelihood estimation of vaccine effectiveness and immune correlates in test-negative design studies with missing data" by Leah I. B. Andrews, Lars van der Laan, and Peter B. Gilbert.

This code was developed using R version 4.4.0.

#### Contents


`README.md`: Describes files  
\


**Simulation Folder**: Contains files necessary to conduct the two-phase TND immune correlates simulation study

- `params.R`: Generates toy RCT dataset to obtain covariate distribution from (i.e., plasmode simulations) and describes additional parameters for the data-generating mechanism

- `methods.R`: Includes functions for the data-generating mechanism, TMLE under semiparametric logistic regression method, pseudo-likelihood logistic regression, and ordinary logistic regression

- `run_sim.R`: Runs simulation using high performance computing cluster and calls `params.R` and `methods.R`

- `sim_results.R`: Generates figures and tables of simulation results from the manuscript and supplement 
    
  
    \
    
**Data Application Folder**: Contains files necessary to assess COVID-19 vaccine effectiveness and antibody marker correlates of COVID-19 from TND study cohorts derived from RCT data

- `COVID-19 Vaccine Effectiveness Application.Rmd`: Implements the COVID-19 vaccine effectiveness data application described in the manuscript 

- `COVID-19-Vaccine-Effectiveness-Application.pdf`: Knitted PDF of `COVID-19 Vaccine Effectiveness Application.Rmd` file

- `COVID VE TMLE Super Learner Model 041326.rds` (not uploaded): RDS file of the TMLE estimator that is described in `COVID-19 Vaccine Effectiveness Application.Rmd` (the .Rmd file loads this RDS file rather than running it to minimize compilation time). This file is too large to be included on GitHub but can be run using the code provided in `COVID-19 Vaccine Effectiveness Application.Rmd`.

- `Exposure-Proximal Correlates of COVID-19 Application.Rmd`: Implements the COVID-19 immune marker data application described in the manuscript 

- `Exposure-Proximal-Correlates-of-COVID-19-Application.pdf`: Knitted PDF of `Exposure-Proximal Correlates of COVID-19 Application.Rmd` file
    
    