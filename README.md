# semiparametric-tmle-tnd
Repository for "Targeted maximum likelihood estimation of vaccine effectiveness and immune correlates in test-negative design studies with missing data" by Leah I. B. Andrews, Lars van der Laan, and Peter B. Gilbert.

## Contents

  1. `README.md`: Describes files  

  2. Simulation Folder: Contains files necessary to conduct the two-phase TND immune correlates simulation study
    a). `params.R`: Generates toy RCT dataset to obtain covariate distribution from (i.e., plasmode simulations) and describes additional parameters for the data-generating mechanism 
    
    b). `methods.R`: Includes functions for the data-generating mechanism, TMLE under semiparametric logistic regression method, pseudo-likelihood logistic regression, and ordinary logistic regression
    
  c). `methods.R`: Runs simulation using high performance computing cluster and calls `params.R` and `methods.R`
    
  d). `sim_results.R`: Generates figures and tables of simulation results from the manuscript and supplement 
    
  3. Data Application Folder: Contains files necessary to assess COVID-19 vaccine effectiveness and antibody marker correlates of COVID-19 from TND study cohorts derived from RCT data
  a). `COVID-19 Vaccine Effectiveness Application.Rmd`: Implements the COVID-19 vaccine effectiveness data application described in the manuscript 
  b). `COVID-19 Vaccine Effectiveness Application.pdf`: Knitted PDF of .Rmd file
  c). `COVID VE TMLE Super Learner Model 041326.rds`: RDS file of the TMLE estimator that is described in `COVID-19 Vaccine Effectiveness Application.Rmd` (the .Rmd file loads this RDS file rather than running it to minimize compilation time)
    
    