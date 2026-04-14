### Data-Generating Function and Statistical Methods for: "Targeted maximum likelihood 
# estimation of vaccine effectiveness and immune correlates 
# in test-negative design studies with missing data" 
# Code Created By: Leah I. B. Andrews
# Date: 04/10/26


#### ph2_biased_data_func: Creates Simulated Two-Phase TND Study Dataset ###########

### Description
# A healthcare-seeking population of pop_size is generated using the
# covariate relationships for A, Y, W, and C specified in param.R.
# sample_size participants with C = 1 are sampled without replacement
# from the population to form the phase one TND study cohort.
# A biased two-phase sampling design is applied to the phase one 
# TND participants. All phase one TND cases are included in phase two.
# (samp_ratio * # of cases) are sampled without replacement from the
# TND study cohort as noncases using four sampling probabilities for
# the covariate strata defined by sex and comorbidity status:
#   40% females without comorbidities
#   10% males without comorbidities
#   10% females with comorbidities
#   40% males with comorbidities
# For example, if samp_ratio = 3 and phase one TND study cohort had 100 cases,
# then 100 cases and 300 noncases would be included in phase two TND study cohort.
# There would be 120 female noncases without comorbidities, 30 male noncases 
# without comorbidities, 30 female noncases with comorbidities, and 120 male 
# noncases with comorbidities. 
# Note: The sampling probabilities and covariate strata are fixed. Would need to 
# change the function itself to edit the sampling probabilities or covariate strata.

### Function Inputs
# a_betas: vector of all beta (log OR) coefficients involved in generating
#   exposure variable A in population. 
#   Vector must have glm model coefficient names.
#   (Be careful with categorical variable names, e.g., "RISKGR1At Risk")
# a_formula: string of glm formula to generate exposure variable A
#   (formula must include all variables mentioned in a_betas)
#    e.g., "~ FEMALE + RISKGR1 + FEMALE:RISKGR1 + CALTIME + CALTIMEs2 + CALTIMEs3"
# a_y_beta: integer of beta coefficient for A to generate SARS-CoV-2 infection Y in 
#    population (i.e., log OR of SARS-CoV-2 infection Y comparing exposure 
#    status A = 1 vs. A = 0). Integer can be named or unnamed.
# not_a_y_betas: vector of all beta (log OR) coefficients involved in generating
#    SARS-CoV-2 infection Y in population, except beta coefficient for A (a_y_beta). 
#    Vector must have glm model coefficient names.   
# y_formula: string of glm formula to generate SARS-CoV-2 infection Y
#    (formula must include all variables mentioned in y_betas)
#    e.g., "~ A + FEMALE + RISKGR1 + FEMALE:RISKGR1 + CALTIME + CALTIMEs2 + CALTIMEs3"
# a_w_beta: integer of beta coefficient for A to generate other infection variable W
#    in population (i.e., log OR of other infection W comparing exposure 
#    status A = 1 vs. A = 0). Integer can be named or unnamed.
# not_a_w_betas: vector of all beta (log OR) coefficients involved in generating
#    other infection W in population, except beta coefficient for A (a_w_beta). 
#    Vector must have glm model coefficient names.   
# w_formula: string of glm formula to generate other infection W
#    (formula must include all variables mentioned in w_betas)
#    e.g., "~ A + FEMALE + RISKGR1+ CALTIME + CALTIMEs2 + CALTIMEs4"
# a_c_beta: integer of beta coefficient for A to generate meeting the symptom 
#    definition C in population (i.e., log OR of meeting symptom definition C 
#    comparing exposure status A = 1 vs. A = 0). Integer can be named or unnamed.
# not_a_c_betas: vector of all beta (log OR) coefficients involved in generating
#    meeting symptom definition C in population, except beta coefficient for  
#    A (a_c_beta). Vector must have glm model coefficient names.   
# c_formula: string of glm formula to generate meeting symptom definition C
#   (formula must include all variables mentioned c_betas)
#     e.g., "~Y + W + A + Y*FEMALE + Y*RISKGR1 + W*FEMALE + W*RISKGR1"
# samp_ratio: desired ratio of total noncases per total case in phase two sample, 
#    must be a numeric. All cases are included in phase two sample. If want all
#    phase one noncases included (no missing data in A), set samp_ratio = 100 or 
#    do not specify a samp_ratio (default is include all noncases).
#    If not enough noncases are available in phase one cohort to achieve desired 
#    samp_ratio, then all noncases are selected for phase two and a warning is outputted. 
#     e.g, if want 1:1 case:noncase, then samp_ratio = 1
#     e.g., if want 1:3 case:noncase, then samp_ratio = 3
#     e.g., if want all cases and noncases, then samp_ratio = NA (default)
# pop_size: integer of desired population sample size
# sample_size: integer of desired phase one TND sample size
# rct.df: RCT data set (without any missing data) to obtain covariate distribution
#     from. rct.df must include all covariates mentioned in the formulas. 

### ph2_biased_data_func outputs a list with 3 items:
# Population: data.frame of the simulated population (retains original rct.df covariates)
# Sample: data.frame of the phase one TND study cohort 
#     (subset on D = 1 to obtain phase two TND study cohort)
# DataParams: data.frame of all the beta coefficients (log OR) used to generate 
#     the population (a_beta, a_y_beta, not_a_y_beta, a_w_beta, not_a_w_beta, 
#     a_c_beta, not_a_c_beta). Each beta coefficient represents the log OR of the  
#     column variable (variable being generated) comparing one unit difference in 
#     the row variable. noncase_d_probs column describes the noncases in phase two:
#       (Intercept) = Number of male noncases without comorbidities
#       FEMALE = Number of female noncases without comorbidities
#       FEMALE:RISKGR1At Risk = Number female noncases with comorbidities 
#       RISKGR1At Risk = Number male noncases with comorbidities
#       Y = desired case:noncase ratio (samp_ratio)
#       Y:RISKGR1At Risk = actual case:noncase ratio 

ph2_biased_data_func <- function(a_betas,a_formula,
                                 a_y_beta, not_a_y_betas, y_formula,
                                 a_w_beta, not_a_w_betas,w_formula,
                                 a_c_beta, not_a_c_betas,c_formula,
                                 samp_ratio = 100,
                                 pop_size,sample_size,rct.df){
  
  # Resample full RCT data with replacement to get new population
  pop_ids <- sample(1:nrow(rct.df), pop_size, replace = TRUE)
  simpop.df <- rct.df[pop_ids, ]
  
  # All A's are observed in the population
  simpop.df$D <- 1
  
  ### Create exposure status A
  
  # Design matrix for A
  x_a_matrix <- model.matrix(eval(parse(text = a_formula)), simpop.df)
  
  # Order a_betas to match design matrix
  reorder_a_betas <- a_betas[match(colnames(x_a_matrix), names(a_betas))]
  

  # Logistic regression to generate A
  simpop.df$A <- rbinom(pop_size , 1, plogis(x_a_matrix %*% reorder_a_betas))
 
  # Saving coefficients that generated A for later
  a_betas.df <- data.frame(variable = names(reorder_a_betas),
                          a_beta = reorder_a_betas)
  
  ### Create SARS-CoV-2 Infection Y
  
  # Design matrix for Y
  x_y_matrix <- model.matrix(eval(parse(text = y_formula)), simpop.df)
  
  # Order y_betas to match design matrix
  names(a_y_beta) <- "A"
  y_betas <- c(a_y_beta, not_a_y_betas)
  reorder_y_betas <- y_betas[match(colnames(x_y_matrix), names(y_betas))]
  
  # Logistic regression to generate Y
  simpop.df$Y <- rbinom(pop_size, 1, plogis(x_y_matrix %*% reorder_y_betas))

  # Saving coefficients that generated Y for later
  y_betas.df <- data.frame(variable = names(reorder_y_betas),
                          y_beta = reorder_y_betas)
  
  ### Create Other Infection W
  
  # Design matrix for Non-SARS Infection W
  x_w_matrix <- model.matrix(eval(parse(text = w_formula)), simpop.df) 
  
  # Order w_betas to match design matrix
  names(a_w_beta) <- "A" 
  w_betas <- c(a_w_beta, not_a_w_betas)
  reorder_w_betas <- w_betas[match(colnames(x_w_matrix), names(w_betas))]
  
  # Logistic Regression to generate W
  simpop.df$W <- rbinom(pop_size, 1, plogis(x_w_matrix %*% reorder_w_betas))
  
  # Saving coefficients that generated W for later
  w_betas.df <- data.frame(variable = names(reorder_w_betas),
                          w_beta = reorder_w_betas)
  
  ### Create meeting symptom definition C (labelled D in manuscript)
  
  # Design matrix for C
  x_c_matrix <- model.matrix(eval(parse(text = c_formula)), simpop.df)
  
  # Order c_betas to match design matrix
  names(a_c_beta) <- "A"
  c_betas <- c(a_c_beta, not_a_c_betas)
  reorder_c_betas <- c_betas[match(colnames(x_c_matrix), names(c_betas))]
  
  # Logistic regression to generate C
  simpop.df$C <- rbinom(pop_size, 1, plogis(x_c_matrix %*% reorder_c_betas))
  
  # Saving coefficients that generated C for later
  c_betas.df <- data.frame(variable = names(reorder_c_betas),
                          c_beta = reorder_c_betas) 
  
  ### Create phase one TND cohort (S = 1)
  simpop.df$S <- 0
  
  if(sum(simpop.df$C == 1) < sample_size){
    sample_size <- sum(simpop.df$C == 1)
    warning(paste0("Fewer symptomatic individuals than desired sample size. ",
                   sample_size," individuals sampled instead.") )}
  
  samp_ids <- sample(which(simpop.df$C == 1), sample_size, replace=FALSE) 
  simpop.df[samp_ids,"S"] <- 1
  sample.df <- simpop.df[samp_ids, ]
  
  
  ### Create phase two TND cohort (who gets exposure status A measured, 
  # Delta = 1). In this code, D is shorthand for Delta and C is symptoms
  
  # Define four covariate strata for two-phase sampling design
  sample.df$STRATA <- interaction(sample.df$FEMALE, sample.df$RISKGR1)
  
  # First initialize number of cases and noncases in phase one cohort
  nsampcase <- sum(sample.df$Y)
  n_0nocom_noncase <- sum(sample.df$Y == 0 & sample.df$STRATA == "0.Not At Risk")
  n_1nocom_noncase <- sum(sample.df$Y == 0 & sample.df$STRATA == "1.Not At Risk")
  n_0com_noncase <- sum(sample.df$Y == 0 & sample.df$STRATA == "0.At Risk")
  n_1com_noncase <- sum(sample.df$Y == 0 & sample.df$STRATA == "1.At Risk")
  
  # If user provides samp_ratio, then measure A only on proportion of noncases.
  # If user does not provide samp_ratio or chooses samp_ratio = 100, 
  # then A measured on all phase one participants.
  if(samp_ratio != 100){
    
    # Initially, no participants in phase two cohort
    sample.df$D <- 0
    
    # All cases are in phase two TND cohort
    sample.df[sample.df$Y == 1 , "D"]  <- 1
    
    # Number of noncases that can be obtained from each stratum when 
    # following the biased noncase sampling ratio (0.4, 0.1, 0.1, 0.4)
    # These numbers may differ than the expected numbers according to the sampling ratio
    n_0com_noncase <- min( floor(nsampcase * samp_ratio * .4), 
                           sum(sample.df$Y==0 & sample.df$STRATA == "0.At Risk")) 
    n_1com_noncase <- min( ceiling(nsampcase * samp_ratio * .5) - n_0com_noncase , 
                           sum(sample.df$Y==0 & sample.df$STRATA == "1.At Risk")) 
    n_0nocom_noncase <- min( floor(nsampcase * samp_ratio * .1) , 
                             sum(sample.df$Y==0 & sample.df$STRATA == "0.Not At Risk")) 
    n_1nocom_noncase <- min( ceiling(nsampcase * samp_ratio)- 
                               n_0com_noncase-n_1com_noncase-n_0nocom_noncase, 
                             sum(sample.df$Y==0 & sample.df$STRATA == "1.Not At Risk")) 
    
    # Throw warning when actual noncase strata sampling ratio isn't met 
    # (i.e., not enough participants in at least one phase one stratum)
    if( n_0com_noncase <  floor(nsampcase * samp_ratio * .4) |  
        n_0nocom_noncase <  floor(nsampcase * samp_ratio * .1)|
        n_1com_noncase <  ceiling(nsampcase * samp_ratio * .1) ){
      
      warning(paste0("Noncase Sampling Ratio Not Met: Actual (Expected) Male Comorbid Noncases = ", n_0com_noncase, " (",  
             floor(nsampcase * samp_ratio * .4), "), Female Comorbid Noncases = ", n_1com_noncase, " (",  
             ceiling(nsampcase * samp_ratio * .1), "), Male NonComorbid Noncases = ", n_0nocom_noncase, " (",  
             floor(nsampcase * samp_ratio * .1), "), Female NonComorbid Noncases = ", n_1nocom_noncase, " (",  
             ceiling(nsampcase * samp_ratio * .4),")" ))
    }
    
    # Throw warning when case:noncase sampling ratio not obtained (too few noncases)
    if(n_0com_noncase + n_0nocom_noncase + n_1com_noncase + n_1nocom_noncase < nsampcase * samp_ratio){
      warning(paste0("Case:Noncase Sampling Ratio Not Obtained (Actual Ratio = 1:",
             round((n_0com_noncase + n_0nocom_noncase +n_1com_noncase + n_1nocom_noncase)/nsampcase,2),
             " with ", n_0com_noncase + n_0nocom_noncase +n_1com_noncase + n_1nocom_noncase," Noncases",
             ", Desired Ratio = 1:",samp_ratio, " with ",nsampcase * samp_ratio," Noncases."))
    }
    
    # Sample without replacement from phase one strata
    p2_0com_noncase_ids <- sample( which(sample.df$Y==0 & sample.df$STRATA == "0.At Risk"),
                                   n_0com_noncase, replace=FALSE) 
    p2_1com_noncase_ids <- sample( which(sample.df$Y==0 & sample.df$STRATA == "1.At Risk"),
                                   n_1com_noncase, replace=FALSE) 
    p2_0nocom_noncase_ids <- sample( which(sample.df$Y==0 & sample.df$STRATA == "0.Not At Risk"),
                                     n_0nocom_noncase, replace=FALSE) 
    p2_1nocom_noncase_ids <- sample( which(sample.df$Y==0 & sample.df$STRATA == "1.Not At Risk"),
                                     n_1nocom_noncase, replace=FALSE) 
    
    sample.df[c(p2_0com_noncase_ids, p2_0nocom_noncase_ids,
                p2_1com_noncase_ids, p2_1nocom_noncase_ids), "D"] <- 1
  }
  
  # Saving information on noncases for later
  noncase_d_probs.df <- data.frame(variable =  c("(Intercept)","RISKGR1At Risk","FEMALE",
                                                 "FEMALE:RISKGR1At Risk","Y","Y:RISKGR1At Risk"),
                                   noncase_d_probs = c(
                                     
                                     # male noncomorbid noncases
                                     n_0nocom_noncase,
                                     
                                     # male comorbid noncases
                                     n_0com_noncase,
                                     
                                     # female noncomorbid noncases
                                     n_1nocom_noncase,
                                     
                                     # female comorbid noncases
                                     n_1com_noncase,
                                     
                                     # expected case:noncase ratio
                                     samp_ratio,
                                     
                                     # actual case:noncase ratio
                                     (n_0com_noncase + n_0nocom_noncase + 
                                        n_1com_noncase + n_1nocom_noncase)/nsampcase))
  
  
  # Collect all data-generating parameters and information
  
  param.df <- merge(a_betas.df, y_betas.df,  all.x = TRUE, all.y = TRUE)
  param.df <- merge(param.df, w_betas.df,  all.x = TRUE, all.y = TRUE)
  param.df <- merge(param.df, c_betas.df,  all.x = TRUE, all.y = TRUE)
  param.df <- merge(param.df, noncase_d_probs.df,  all.x = TRUE, all.y = TRUE)
  
  ### Output list of simulated population, simulated phase one TND study 
  # cohort, and data-generating information
  return(list(Population = simpop.df,
              Sample = sample.df,
              DataParams = param.df))
}

#### ordlogit_func: Runs Ordinary Logistic Regression ########

### Function inputs:
# poi_var: string of exposure of interest variable name
#  e.g., "A"
# outcome_var: string of outcome variable name
#  e.g., "Y"
# model_formula: string of desired logistic regression glm formula
#  e.g., "Y ~ A + FEMALE + RISKGR1 + CALTIME"
# missing_var: string of indicator of missing the exposure variable name
#  e.g., "D"
# oracle_formula: optional string of true data-generating glm formula for 
#  outcome variable (more relevant when no missing data). Included so user 
#  can choose (using oracle_flag) if a logistic regression should be fit
#  with the model_formula or oracle_formula. If no input is given, 
#  oracle_formula = model_formula
# oracle_flag: optional T/F flag indicating if the oracle_formula should 
#  be used in the ordinary logistic regression instead of the model_formula.
#  oracle_flag = F is default (i.e., model_formula is used)
# learner: unused argument included to give identical argument structure
#  for all statistical methods evaluated in the simulation study.
# df: dataset with complete information on all variables in model_formula
#  except poi_var
#  e.g., sim.df

### ordlogit_func outputs:
# list that contains information about the fitted logistic regression
# and about df (collected as sanity checks in the simulations)
ordlogit_func <- function(poi_var, outcome_var, model_formula,
                          missing_var = "D",
                          oracle_formula = model_formula,
                          oracle_flag = F,learner = "glm", df){
  
  # Keep individuals with observed exposure variable 
  # (D = 1; D is shorthand for Delta)
  complete.df <- df %>% filter(get(missing_var) == 1)
  
  # Run model using the oracle formula
  if(oracle_flag==T){
    model_formula <- oracle_formula
  }
  
  # Run logistic regression model
  mod <- glm(eval(parse(text = model_formula)), 
             data = complete.df, family = "binomial")
  glm.mod <- summary(mod)
  
  # Calculating relevant statistics
  coef <- glm.mod$coefficients[poi_var, "Estimate"]
  se <- glm.mod$coefficients[poi_var, "Std. Error"]
  ve <- (1 - exp(coef))*100
  z <- qnorm(0.975)
  veCIlow <- (1-exp(coef + z*se))*100
  veCIhigh <- (1-exp(coef - z*se))*100
  pval = glm.mod$coefficients[poi_var, "Pr(>|z|)"]
  
  # Large list of coefficients and information about df (that was
  # collected as a sanity check in the simulations)
  return(list(
        StatMethod = "Ordinary", 
        Beta = coef,
        SE = se,
        VE = ve,
        VELL = veCIlow,
        VEUL = veCIhigh,
        PValue = pval,
        Reject = as.integer(pval < 0.05),
        nsamp = nrow(df),
        ncomplete = nrow(complete.df),
        nVac = with(complete.df, sum(get(outcome_var) == 1)),
        nCase = with(complete.df, sum(get(poi_var) == 1)),
        nW = with(complete.df, sum(W == 1)),
        nC = with(complete.df, sum(C == 1)),
        nVacCase = with(complete.df, sum(get(outcome_var) == 1 & get(poi_var) == 1)),
        nUnvacCase = with(complete.df, sum(get(outcome_var) == 0 & get(poi_var) == 1)),
        nVacCaseCom = with(complete.df, sum(get(outcome_var) == 1 & 
                                            get(poi_var) == 1 & RISKGR1 == "At Risk")),
        nUnvacCaseCom = with(complete.df, sum(get(outcome_var) == 0 & 
                                              get(poi_var) == 1 & RISKGR1 == "At Risk")),
        nVacCaseNoCom = with(complete.df, sum(get(outcome_var) == 1 & 
                                              get(poi_var) == 1 & 
                                                RISKGR1 == "Not At Risk")),
        nUnvacCaseNoCom = with(complete.df, sum(get(outcome_var) == 0 & 
                                                get(poi_var) == 1 & 
                                                  RISKGR1 == "Not At Risk")),
        nVacNoncaseCom = with(complete.df, sum(get(outcome_var) == 1 & 
                                               get(poi_var) == 0 & 
                                                 RISKGR1 == "At Risk")),
        nUnvacNoncaseCom = with(complete.df, sum(get(outcome_var) == 0 & 
                                                 get(poi_var) == 0 & 
                                                   RISKGR1 == "At Risk")),
        nVacNoncaseNoCom = with(complete.df, sum(get(outcome_var) == 1 & 
                                                 get(poi_var) == 0 & 
                                                   RISKGR1 == "Not At Risk")),
        nUnvacNoncaseNoCom = with(complete.df, sum(get(outcome_var) == 0 & 
                                                   get(poi_var) == 0 & 
                                                     RISKGR1 == "Not At Risk")),
        nFemCom = with(complete.df, sum(FEMALE == 1 & RISKGR1 == "At Risk")),
        nMaleCom = with(complete.df, sum(FEMALE == 0 & RISKGR1 == "At Risk")),
        nFemNoCom = with(complete.df, sum(FEMALE == 1 & RISKGR1 == "Not At Risk")),
        nMaleNoCom = with(complete.df, sum(FEMALE == 0 & RISKGR1 == "Not At Risk")),
        VacObs = with(complete.df, sum(get(outcome_var) == 1))/with(df, sum(get(outcome_var) == 1)),
        UnvacObs = with(complete.df, sum(get(outcome_var)==0))/with(df, sum(get(outcome_var) == 0)),
        ComObs = with(complete.df, sum(RISKGR1 == "At Risk"))/with(df, sum(RISKGR1 == "At Risk")),
        NoComObs = with(complete.df, sum(RISKGR1 == "Not At Risk"))/with(df, sum(RISKGR1 == "Not At Risk")),
        FemObs = with(complete.df, sum(FEMALE == 1))/with(df, sum(FEMALE == 1)),
        MaleObs = with(complete.df, sum(FEMALE == 0))/with(df, sum(FEMALE == 0)),
        missing_var = missing_var, 
        Formula = model_formula,
        Model = mod)) }




#### twoph_e_logit_func: Two-Phase Logistic Regression Approach w/ Empirical Variance ####

### Function inputs:
# poi_var: string of exposure of interest variable name
#  e.g., "A"
# outcome_var: string of outcome variable name
#  e.g., "Y"
# model_formula: string of desired logistic regression formula
#  e.g., "Y ~ A + FEMALE + RISKGR1 + CALTIME"
# missing_var: string of indicator of missing the exposure variable name
#  e.g., "D"
# oracle_formula: optional string of true data-generating formula for 
#  outcome variable (more relevant when no missing data). Included so user 
#  can choose (using oracle_flag) if a logistic regression should be fit
#  with the model_formula or oracle_formula. If no input is given, 
#  oracle_formula = model_formula
# oracle_flag: optional T/F flag indicating if the oracle_formula should 
#  be used in the two-phase logistic regression instead of the model_formula.
#  oracle_flag = F is default (i.e., model_formula is used)
# learner: string of method of estimation allowed from tps function 
#   e.g., pseudo-likelihood ("PL"), weighted likelihood ("WL"), or 
#   maximum likelihood ("ML")
# df: dataset with complete information on all variables in model_formula
#  except poi_var and a variable named STRATA that represents the covariate 
#  strata involved in two-phase sampling design. The covariate strata will 
#  be accounted for using an offset.
#  e.g., sim.df

### twoph_e_logit_func outputs:
# list that contains information about the fitted pseudo-likelihood 
# logistic regression and about df (collected as sanity checks 
# in the simulations)

twoph_e_logit_func <- function(poi_var = "A", outcome_var = "Y",
                                  model_formula,missing_var = "D",
                           oracle_formula = model_formula, oracle_flag = F,
                           learner = "PL", df){
  
  # Table of STRATA and cases
  strata_tbl <- table(df$STRATA, df$Y)
  
  
  # Keep individuals with observed exposure 
  # (D = 1; D is shorthand for Delta)
  complete.df <- df %>% filter(get(missing_var) == 1)
  
  # Run model using the oracle formula
  if(oracle_flag == T){
    model_formula <- oracle_formula
  }
  
  # Run logistic regression model
  mod <- tps(eval(parse(text = model_formula)), data = complete.df,
      nn0 = strata_tbl[ , "0"], nn1 = strata_tbl[ ,"1"],
      group = complete.df$STRATA, method = learner, cohort = FALSE)

  
  # Calculating relevant statistics
  coef <-  mod$coef[which(names(mod$coef) == poi_var)]
  se <- sqrt( mod$cove[poi_var, poi_var] ) # Empirical Variance
  ve <- (1 - exp(coef))*100
  z <-qnorm(0.975)
  veCIlow <- (1 - exp(coef + z*se))*100
  veCIhigh <- (1 - exp(coef - z*se))*100
  pval = 2*(1 - pnorm(abs(coef/se)))
  
  # Large list of coefficients and information about df (that was
  # collected as a sanity check in the simulations)
  return(list(StatMethod = paste0(learner, " Logit EmpSE"), 
        Beta = coef,
        SE = se,
        VE = ve,
        VELL = veCIlow,
        VEUL = veCIhigh,
        PValue = pval,
        Reject = as.integer(pval<0.05),
        nsamp = nrow(df),
        ncomplete = nrow(complete.df),
        nVac = with(complete.df, sum(get(outcome_var) == 1)),
        nCase = with(complete.df, sum(get(poi_var) == 1)),
        nW = with(complete.df, sum(W == 1)),
        nC = with(complete.df, sum(C == 1)),
        nVacCase = with(complete.df, sum(get(outcome_var) == 1 & get(poi_var) == 1)),
        nUnvacCase = with(complete.df, sum(get(outcome_var) == 0 & get(poi_var) == 1)),
        nVacCaseCom = with(complete.df, sum(get(outcome_var) == 1 &
                                            get(poi_var) == 1&RISKGR1 == "At Risk")),
        nUnvacCaseCom = with(complete.df, sum(get(outcome_var) == 0 &
                                              get(poi_var) == 1 & RISKGR1 == "At Risk")),
        nVacCaseNoCom = with(complete.df, sum(get(outcome_var) == 1 &
                                              get(poi_var) == 1&RISKGR1 == "Not At Risk")),
        nUnvacCaseNoCom = with(complete.df, sum(get(outcome_var) == 0 &
                                                get(poi_var) == 1 &
                                                  RISKGR1 == "Not At Risk")),
        nVacNoncaseCom = with(complete.df, sum(get(outcome_var) == 1 &
                                               get(poi_var) == 0 &
                                                 RISKGR1 == "At Risk")),
        nUnvacNoncaseCom = with(complete.df, sum(get(outcome_var) == 0 &
                                                 get(poi_var) == 0 &
                                                   RISKGR1 == "At Risk")),
        nVacNoncaseNoCom = with(complete.df, sum(get(outcome_var) == 1 &
                                                 get(poi_var) == 0 &
                                                   RISKGR1 == "Not At Risk")),
        nUnvacNoncaseNoCom = with(complete.df, sum(get(outcome_var) == 0 &
                                                   get(poi_var) == 0 &
                                                     RISKGR1 == "Not At Risk")),
        nFemCom = with(complete.df, sum(FEMALE == 1&RISKGR1 == "At Risk")),
        nMaleCom = with(complete.df, sum(FEMALE == 0&RISKGR1 == "At Risk")),
        nFemNoCom = with(complete.df, sum(FEMALE == 1&RISKGR1 == "Not At Risk")),
        nMaleNoCom = with(complete.df, sum(FEMALE == 0&RISKGR1 == "Not At Risk")),
        VacObs = with(complete.df, sum(get(outcome_var) == 1))/with(df, sum(get(outcome_var) == 1)),
        UnvacObs = with(complete.df, sum(get(outcome_var) == 0))/with(df, sum(get(outcome_var) == 0)),
        ComObs = with(complete.df, sum(RISKGR1 == "At Risk"))/with(df, sum(RISKGR1 == "At Risk")),
        NoComObs = with(complete.df, sum(RISKGR1 == "Not At Risk"))/with(df, sum(RISKGR1 == "Not At Risk")),
        FemObs = with(complete.df, sum(FEMALE == 1))/with(df, sum(FEMALE == 1)),
        MaleObs = with(complete.df, sum(FEMALE == 0))/with(df, sum(FEMALE == 0)),
        missing_var = missing_var, 
        Formula = model_formula,
        Model = mod)) }

#### twoph_m_logit_func: Two-Phase Logistic Regression Approach w/ Model Variance ####

### Function inputs:
# poi_var: string of exposure of interest variable name
#  e.g., "A"
# outcome_var: string of outcome variable name
#  e.g., "Y"
# model_formula: string of desired logistic regression formula
#  e.g., "Y ~ A + FEMALE + RISKGR1 + CALTIME"
# missing_var: string of indicator of missing the exposure variable name
#  e.g., "D"
# oracle_formula: optional string of true data-generating formula for 
#  outcome variable (more relevant when no missing data). Included so user 
#  can choose (using oracle_flag) if a logistic regression should be fit
#  with the model_formula or oracle_formula. If no input is given, 
#  oracle_formula = model_formula
# oracle_flag: optional T/F flag indicating if the oracle_formula should 
#  be used in the two-phase logistic regression instead of the model_formula.
#  oracle_flag = F is default (i.e., model_formula is used)
# learner: string of method of estimation allowed from tps function 
#   e.g., pseudo-likelihood ("PL"), weighted likelihood ("WL"), or 
#   maximum likelihood ("ML")
# df: dataset with complete information on all variables in model_formula
#  except poi_var and a variable named STRATA that represents the covariate 
#  strata involved in two-phase sampling design. The covariate strata will 
#  be accounted for using an offset.
#  e.g., sim.df

### twoph_e_logit_func outputs:
# list that contains information about the fitted pseudo-likelihood 
# logistic regression and about df (collected as sanity checks 
# in the simulations)

twoph_m_logit_func <- function(poi_var = "A", outcome_var = "Y",
                                  model_formula, missing_var = "D",
                                  oracle_formula = model_formula, oracle_flag = F,
                                  learner = "PL", df){
  
  # Table of STRATA and cases
  strata_tbl <- table(df$STRATA, df$Y)
  
  
  # Keep individuals with observed exposure 
  # (D = 1; D is shorthand for Delta)
  complete.df <- df %>% filter(get(missing_var) == 1)
  
  # Run model using the oracle formula
  if(oracle_flag == T){
    model_formula <- oracle_formula
  }
  
  ## Run logistic regression model
  mod <- tps(eval(parse(text = model_formula)), data = complete.df,
             nn0 = strata_tbl[ , "0"], nn1 = strata_tbl[ ,"1"],
             group = complete.df$STRATA, method = learner, cohort = FALSE)
  
  
  ## Calculating relevant statistics
  coef <-  mod$coef[which(names(mod$coef) == poi_var)]
  se <- sqrt( mod$covm[poi_var, poi_var] ) # Model Variance
  ve <- (1 - exp(coef))*100
  z <-qnorm(0.975)
  veCIlow <- (1 - exp(coef + z*se))*100
  veCIhigh <- (1 - exp(coef - z*se))*100
  pval = 2*(1 - pnorm(abs(coef/se)))
  
  # Large list of coefficients and information about df (that was
  # collected as a sanity check in the simulations)
  return(list(StatMethod = paste0(learner, " Logit ModSE"), 
        Beta = coef,
        SE = se,
        VE = ve,
        VELL = veCIlow,
        VEUL = veCIhigh,
        PValue = pval,
        Reject = as.integer(pval < 0.05),
        nsamp = nrow(df),
        ncomplete = nrow(complete.df),
        nVac = with(complete.df, sum(get(outcome_var) == 1)),
        nCase = with(complete.df, sum(get(poi_var) == 1)),
        nW = with(complete.df, sum(W == 1)),
        nC = with(complete.df, sum(C == 1)),
        nVacCase = with(complete.df, sum(get(outcome_var) == 1 & get(poi_var) == 1)),
        nUnvacCase = with(complete.df, sum(get(outcome_var) == 0 & get(poi_var) == 1)),
        nVacCaseCom = with(complete.df, sum(get(outcome_var) == 1 &
                                              get(poi_var) == 1 & RISKGR1 == "At Risk")),
        nUnvacCaseCom = with(complete.df, sum(get(outcome_var) == 0 &
                                                get(poi_var) == 1 & 
                                                RISKGR1 == "At Risk")),
        nVacCaseNoCom = with(complete.df, sum(get(outcome_var) == 1 & 
                                                get(poi_var) == 1 &
                                                RISKGR1 == "Not At Risk")),
        nUnvacCaseNoCom = with(complete.df, sum(get(outcome_var) == 0 &
                                                  get(poi_var) == 1 &
                                                  RISKGR1 == "Not At Risk")),
        nVacNoncaseCom = with(complete.df, sum(get(outcome_var) == 1 &
                                                 get(poi_var) == 0 &
                                                 RISKGR1 == "At Risk")),
        nUnvacNoncaseCom = with(complete.df, sum(get(outcome_var) == 0 &
                                                   get(poi_var) == 0 &
                                                   RISKGR1 == "At Risk")),
        nVacNoncaseNoCom = with(complete.df, sum(get(outcome_var) == 1 &
                                                   get(poi_var) == 0 &
                                                   RISKGR1 == "Not At Risk")),
        nUnvacNoncaseNoCom = with(complete.df, sum(get(outcome_var) == 0 &
                                                     get(poi_var) == 0 &
                                                     RISKGR1 == "Not At Risk")),
        nFemCom = with(complete.df, sum(FEMALE == 1 & RISKGR1 == "At Risk")),
        nMaleCom = with(complete.df, sum(FEMALE == 0 & RISKGR1 == "At Risk")),
        nFemNoCom = with(complete.df, sum(FEMALE == 1 & RISKGR1 == "Not At Risk")),
        nMaleNoCom = with(complete.df, sum(FEMALE == 0 & RISKGR1 == "Not At Risk")),
        VacObs = with(complete.df, sum(get(outcome_var) == 1))/with(df, sum(get(outcome_var) == 1)),
        UnvacObs = with(complete.df, sum(get(outcome_var) == 0))/with(df, sum(get(outcome_var) == 0)),
        ComObs = with(complete.df, sum(RISKGR1 == "At Risk"))/with(df, sum(RISKGR1 == "At Risk")),
        NoComObs = with(complete.df, sum(RISKGR1 == "Not At Risk"))/with(df, sum(RISKGR1 == "Not At Risk")),
        FemObs = with(complete.df, sum(FEMALE == 1))/with(df, sum(FEMALE == 1)),
        MaleObs = with(complete.df, sum(FEMALE == 0))/with(df, sum(FEMALE == 0)),
        missing_var = missing_var, 
        Formula = model_formula,
        Model = mod)) }


#### semilogit_func: TMLE Approach under Semiparametric Logistic Regression Model ####
# with estimation using partially linear first-order smooth 
# highly adaptive lasso that allows for two-way covariate interactions.
# TMLE approach constructed from spglm function in causalglm package

### Function inputs:
# poi_var: string of exposure of interest variable name
#  e.g., "A"
# outcome_var: string of outcome variable name
#  e.g., "Y"
# model_formula: string of desired logistic regression formula
#  the function will modify the formula to match spglm formula requirements
#  e.g., "Y ~ A + FEMALE + RISKGR1 + CALTIME"
# missing_var: string of indicator of missing the exposure variable name
#  e.g., "D"
# oracle_formula: optional string of true data-generating formula for 
#  outcome variable (more relevant when no missing data). Included so user 
#  can choose (using oracle_flag) if a logistic regression should be fit
#  with the model_formula or oracle_formula. If no input is given, 
#  oracle_formula = model_formula
# oracle_flag: optional T/F flag indicating if the oracle_formula should 
#  be used in the logistic regression instead of the model_formula.
#  oracle_flag = F is default (i.e., model_formula is used)
# learner: unused argument included to give identical argument structure
#  for all statistical methods evaluated in the simulation study.  
#  Regardless of input, a partially linear first-order smooth HAL that
#  allows for two-way covariate interactions will be used for estimation.
# df: dataset with complete information on all variables in model_formula
#  except poi_var.
#  e.g., sim.df

### semilogit_func outputs:
# list that contains information about the fitted pseudo-likelihood 
# logistic regression and about df (collected as sanity checks 
# in the simulations)

semilogit_func <- function(poi_var, outcome_var, model_formula,
                           missing_var = "D",
                           oracle_formula = model_formula, oracle_flag = F,
                           learner = "HAL", df){
  
  # Keep individuals with observed exposure variable 
  # (D = 1; D is shorthand for Delta)
  complete.df <- df %>% filter(get(missing_var) == 1)
  
  # Run model using the oracle formula
  if(oracle_flag == T){
    model_formula <- oracle_formula
  }
  
  # Use model.matrix to get data in right format (i.e., dataset with 
  # only numeric and binary variables, including outcome)
  model_formula <- gsub(paste0(outcome_var, " *~"),
                                paste0("~", outcome_var, "+"),
                                model_formula)
  
  model.df <- as.data.frame(model.matrix(eval(parse(text = model_formula)),
                                         data = complete.df))
  
  # Run TMLE under partially linear logistic regression model using first-order 
  # smooth highly adaptive lasso that allows for two-way covariate interactions
  mod <- spglm( ~1, model.df,
                W = setdiff(names(model.df), c("(Intercept)", outcome_var, poi_var)), 
                A = poi_var, Y = outcome_var,
                learning_method = learner,
                estimand = "OR",
                HAL_args_Y0W = list(smoothness_orders = 1, 
                                    max_degree = 2, 
                                    num_knots = c(10, 10)))
  
  # Large list of coefficients and information about df (that was
  # collected as a sanity check in the simulations)
  return(list(StatMethod = "Semiparametric", 
        Beta = coef(mod)$tmle_est,
        SE = coef(mod)$se,
        VE = (1-coef(mod)$psi_exp)*100,
        VELL = (1-coef(mod)$upper_exp)*100,
        VEUL = (1-coef(mod)$lower_exp)*100,
        PValue = coef(mod)$p_value,
        Reject = as.integer(coef(mod)$p_value<0.05),
        nsamp = nrow(df),
        ncomplete = nrow(complete.df),
        nVac = with(complete.df, sum(get(outcome_var) == 1)),
        nCase = with(complete.df, sum(get(poi_var) == 1)),
        nW = with(complete.df, sum(W == 1)),
        nC = with(complete.df, sum(C == 1)),
        nVacCase = with(complete.df, sum(get(outcome_var) == 1 & get(poi_var) == 1)),
        nUnvacCase = with(complete.df, sum(get(outcome_var) == 0 & get(poi_var) == 1)),
        nVacCaseCom = with(complete.df, sum(get(outcome_var) == 1 &
                                              get(poi_var) == 1 & RISKGR1 == "At Risk")),
        nUnvacCaseCom = with(complete.df, sum(get(outcome_var) == 0 &
                                                get(poi_var) == 1 &
                                                RISKGR1 == "At Risk")),
        nVacCaseNoCom = with(complete.df, sum(get(outcome_var) == 1 &
                                                get(poi_var) == 1 &
                                                RISKGR1 == "Not At Risk")),
        nUnvacCaseNoCom = with(complete.df, sum(get(outcome_var) == 0 &
                                                  get(poi_var) == 1 &
                                                  RISKGR1 == "Not At Risk")),
        nVacNoncaseCom = with(complete.df, sum(get(outcome_var) == 1 &
                                                 get(poi_var) == 0 &
                                                 RISKGR1 == "At Risk")),
        nUnvacNoncaseCom = with(complete.df, sum(get(outcome_var) == 0 &
                                                   get(poi_var) == 0 &
                                                   RISKGR1 == "At Risk")),
        nVacNoncaseNoCom = with(complete.df, sum(get(outcome_var) == 1 &
                                                   get(poi_var) == 0 &
                                                   RISKGR1 == "Not At Risk")),
        nUnvacNoncaseNoCom = with(complete.df, sum(get(outcome_var) == 0 &
                                                     get(poi_var) == 0 &
                                                     RISKGR1 == "Not At Risk")),
        nFemCom = with(complete.df, sum(FEMALE == 1 & RISKGR1 == "At Risk")),
        nMaleCom = with(complete.df, sum(FEMALE == 0 & RISKGR1 == "At Risk")),
        nFemNoCom = with(complete.df, sum(FEMALE == 1 & RISKGR1 == "Not At Risk")),
        nMaleNoCom = with(complete.df, sum(FEMALE == 0 & RISKGR1 == "Not At Risk")),
        VacObs = with(complete.df, sum(get(outcome_var) == 1))/with(df, sum(get(outcome_var) == 1)),
        UnvacObs = with(complete.df, sum(get(outcome_var) == 0))/with(df, sum(get(outcome_var) == 0)),
        ComObs = with(complete.df, sum(RISKGR1 == "At Risk"))/with(df, sum(RISKGR1 == "At Risk")),
        NoComObs = with(complete.df, sum(RISKGR1 == "Not At Risk"))/with(df, sum(RISKGR1 == "Not At Risk")),
        FemObs = with(complete.df, sum(FEMALE == 1))/with(df, sum(FEMALE == 1)),
        MaleObs = with(complete.df, sum(FEMALE == 0))/with(df, sum(FEMALE == 0)),
        missing_var = missing_var, 
        Formula = model_formula,
        Model = mod)) }


