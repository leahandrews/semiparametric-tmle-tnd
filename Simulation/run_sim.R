## Running Simulation Using High Performance Computing Cluster for: 
# "Targeted maximum likelihood estimation of vaccine effectiveness 
# and immune correlates in test-negative design studies with missing data" 
# Code Created By: Leah I. B. Andrews
# Last Upated: 04/21/26


# Set library path
.libPaths(c("/home/users/landrew2/R_lib"))

library(SimEngine)
library(sl3)


print(Sys.time())

(t1 <- Sys.time())

run_on_cluster(
  first = {
    #file_path <- "/home/users/landrew2/Dissertation/Project1/"
    file_path <- "~/Desktop/UW Stuff/IS TND/Code/TND TMLE/semiparametric-tmle-tnd local/Simulation/"
    source(paste0(file_path,"params.R"))
    source(paste0(file_path,"methods.R"))
    nsims <- 1000
    ncores <- 70
    nseed <-  2024 
    stop_error_flag <- FALSE 
    
    # Population size
    npop <- 50000
    
    sim <- new_sim()
    
    sim %<>% set_levels(
      
      # Comparing seven logistic regression estimators defined in methods.R
      Estimator = list(
        
        # nMLE: Naïve ordinary logistic regression (main effects)
        "OrdLogitRev" = list(func = "ordlogit_func", oracle_flag = F,
                       poi_var = "A", outcome_var = "Y",
                       learner = "glm",
                       missing_var = "D",
                       model_formula = "Y ~ A + FEMALE + RISKGR1 + CALTIME" ),
        
        # MLEx: Ordinary logistic regression(main effects and interaction)
        "OrdLogitRevx" = list(func = "ordlogit_func", oracle_flag = F,
                      poi_var = "A", outcome_var = "Y",
                      learner = "glm",
                      missing_var = "D",
                      model_formula = "Y ~ A + FEMALE + RISKGR1 + FEMALE:RISKGR1 + CALTIME"),
        
        # nPLE: Naïve pseudo-likelihood logistic regression with empirical 
        # variance estimates (main effects)
        "PLLogitE0" = list(func = "twoph_e_logit_func", oracle_flag = F,
                     poi_var = "A", outcome_var = "Y",
                     learner = "PL",
                     missing_var = "D",
                     model_formula = "Y ~ A + FEMALE + RISKGR1 + CALTIME"),
        
        # nPLM: Naïve pseudo-likelihood logistic regression with model 
        # variance estimates (main effects)
        "PLLogitM0" = list(func = "twoph_m_logit_func", oracle_flag = F,
                     poi_var = "A", outcome_var = "Y",
                     learner = "PL",
                     missing_var = "D",
                     model_formula = "Y ~ A + FEMALE + RISKGR1 + CALTIME" ),
        
        # PLEx: Pseudo-likelihood logistic regression with empirical 
        # variance estimates (main effects and interaction)
        "PLLogitEx" = list(func = "twoph_e_logit_func", oracle_flag = F,
                     poi_var = "A", outcome_var = "Y",
                     learner="PL",
                     missing_var = "D",
                     model_formula = "Y ~ A + FEMALE + RISKGR1 + FEMALE:RISKGR1 + CALTIME" ),
        
        # PLMx: Pseudo-likelihood logistic regression with model 
        # variance estimates (main effects and interaction)
        "PLLogitMx" = list(func = "twoph_m_logit_func", oracle_flag = F,
                     poi_var = "A", outcome_var = "Y",
                     learner="PL",
                     missing_var = "D",
                     model_formula = "Y ~ A + FEMALE + RISKGR1 + FEMALE:RISKGR1 + CALTIME" ),
        
        # TMLE: Targeted maximum likelihood estimation under semiparametric 
        # logistic regression. Estimated using highly adaptive lasso.
        "SemiLogit0" = list(func = "semilogit_func", oracle_flag=F,
                    poi_var = "Y", outcome_var = "A",
                    learner = "HAL",
                    missing_var="D",
                    model_formula = "A ~ Y + FEMALE + RISKGR1 + CALTIME")),
      
      # Three parameters for inference on A (will be converted to 
      # beta = log OR coefficients from VE = (1 - logOR_A)*100 )
      trueVE = c(0, 30, 80),
      
      nseed = c(nseed), 
      
      # Three confounding settings:
      distribution = list(
        
        # Scenario 3/Splines Setting: Main Effects + Sex-Comorbidity Interaction + 
        # Calendar Date Splines
        "Scenario 3" = list(
            a_betas = a_betas3,
            y_betas = y_betas3,
            w_betas = w_betas,
            c_betas = c_betas,
            oracle_formula = "A ~ Y + FEMALE + RISKGR1 + FEMALE:RISKGR1 + CALTIME + CALTIMEs2 + CALTIMEs3"),
        
        # Scenario 2/Interaction Setting: Main Effects + Sex-Comorbidity Interaction
        "Scenario 2" = list(a_betas = a_betas2,
            y_betas = y_betas2,
            w_betas = w_betas,
            c_betas = c_betas,
            oracle_formula = "A ~ Y + FEMALE + RISKGR1 + FEMALE:RISKGR1 + CALTIME"),
        
         # Scenario 1/Main Effects Setting: Sex, Comorbidity, Calendar Date Main Effect
         "Scenario 1" = list(a_betas = a_betas1,
            y_betas = y_betas1,
            w_betas = w_betas,
            c_betas = c_betas,
            oracle_formula = "A ~ Y + FEMALE + RISKGR1 + CALTIME")),
     
      # Two-Phase Sampling Design for A 
      # (1:1 case:noncase , 1:3 case:noncase, all participants)
      samprate = c(1, 3, 100), 
      
      # Phase One TND Sample Size
      nsamp = c(2000, 500, 1000, 3000))
    
    sim %<>% set_config(
      num_sim = nsims,
      n_cores = ncores,
      batch_levels = c("trueVE", "distribution", "samprate", "nsamp"),
      return_batch_id = T,
      seed = nseed, 
      stop_at_error = stop_error_flag,
      packages = c("causalglm", "osDesign", "dplyr", "readr", "sl3"))
    
    # Set_script is like a for loop
    sim %<>% set_script(function() {
      
      # Apply all statistical methods to same simulated dataset
      batch({  
        
        # Simulate data
        data <- ph2_biased_data_func(a_betas = L$distribution$a_betas, 
                                     a_formula = a_formula,
                                     a_y_beta = log(1 - L$trueVE/100),
                                     not_a_y_betas = L$distribution$y_betas,
                                     y_formula = y_formula,
                                     a_w_beta = imm_w_beta, 
                                     not_a_w_betas = L$distribution$w_betas,
                                     w_formula = w_formula,
                                     a_c_beta = imm_c_beta, 
                                     not_a_c_betas = L$distribution$c_betas,
                                     c_formula = c_formula,
                                     samp_ratio = L$samprate,
                                     pop_size = npop, 
                                     sample_size = L$nsamp, 
                                     rct.df = simtest)
      })
      # Apply each statistical method to same simulated dataset 
      # (this is why we need identical argument structure)
      estimates <- use_method(L$Estimator$func, 
                              list(poi_var = L$Estimator$poi_var,
                                   outcome_var = L$Estimator$outcome_var,
                                   model_formula = L$Estimator$model_formula, 
                                   missing_var = L$Estimator$missing_var,
                                   oracle_formula = L$distribution$oracle_formula,
                                   oracle_flag = L$Estimator$oracle_flag,
                                   learner = L$Estimator$learner,
                                   df = data$Sample))
      
      # Output large list of coefficients and information about simulated data
      return(list(
        "truelogOR" = log(1 - L$trueVE/100),
        "logOR" = estimates$Beta,
        "logORSE" = estimates$SE,
        "logORVar" = (estimates$SE)^2,
        "VE" = estimates$VE,
        "VELL" = estimates$VELL,
        "VEUL" = estimates$VEUL,
        "PValue"= estimates$PValue,
        "Reject"= estimates$Reject,
        "PopMargA" = mean(data$Population$A),
        "PopMargY" = mean(data$Population$Y),
        "PopMargW" = mean(data$Population$W),
        "PopMargC" = mean(data$Population$C),
        "PopMargFem" = mean(data$Population$FEMALE),
        "PopMargCom" = prop.table(table(data$Population$RISKGR1))[2],
        "Phase1MargA" = mean(data$Sample$A),
        "Phase1MargY" = mean(data$Sample$Y),
        "Phase1MargW" = mean(data$Sample$W),
        "Phase1MargC" = mean(data$Sample$C),
        "Phase1MargFem" = mean(data$Sample$FEMALE),
        "Phase1MargCom" = prop.table(table(data$Sample$RISKGR1))[2],
        "Phase2MargA" = estimates$nVac/estimates$ncomplete,
        "Phase2MargY" = estimates$nCase/estimates$ncomplete,
        "Phase2MargW" = estimates$nW/estimates$ncomplete,
        "Phase2MargC" = estimates$nC/estimates$ncomplete,
        "Phase2MargnA" = estimates$nVac,
        "Phase2MargnY" = estimates$nCase,
        "SampMargnW" = estimates$nW,
        "SampMargnC" = estimates$nC,
        "SampVacCase" = estimates$nVacCase, 
        "SampUnvacCase" = estimates$nUnvacCase,
        
        # Total cases with comorbidities in phase one
        "SampCaseComTot" = with(data$Sample, sum( RISKGR1 == "At Risk" & Y == 1 )), 
        
        # Total cases without comorbidities in phase one
        "SampCaseNoComTot" = with(data$Sample,sum( RISKGR1 == "Not At Risk" & Y == 1 )), 
        
        # Total noncases with comorbidities in phase one
        "SampNoncaseComTot" = with(data$Sample, sum( RISKGR1 == "At Risk" & Y == 0 )), 
        
        # Total noncases without comorbidities in phase one
        "SampNoncaseNoComTot" = with(data$Sample,sum( RISKGR1 == "Not At Risk" & Y == 0 )), 
        
        # Total noncases with comorbidities in phase two
        "SampNoncaseComObs" = with(data$Sample, sum( D == 1 & RISKGR1 == "At Risk" & Y == 0 )), 
        
        # Total noncases without comorbidities in phase two
        "SampNoncaseNoComObs" = with(data$Sample, sum( D == 1 & RISKGR1 == "Not At Risk" & Y == 0 )), 
        
        "ObsSampRatio" = tail(data$DataParams$noncase_d_probs, n=1),
        
        # 1 - missing rate
        "SampObsProp"= estimates$ncomplete/estimates$nsamp, 
        
        # Proportion of phase one participants with A = 1 in phase two
        "SampVacObs" = estimates$VacObs,  
        
        # Proportion of phase one participants with A = 0 that are in phase two
        "SampUnvacObs" = estimates$UnvacObs, 
        
        # Proportion of phase one participants with comorbidities that are in phase two
        "SampComObs" = estimates$ComObs,  
        
        # Proportion of phase one participants without comorbidities that are in phase two
        "SampNoComObs" = estimates$NoComObs, 
        
        # Proportion of phase one females that are in phase two
        "SampFemObs" = estimates$FemObs, 
        
        # Proportion of phase one males that are in phase two
        "SampMaleObs" = estimates$MaleObs, 
        
        "SampnFemCom" = estimates$nFemCom,
        "SampnMaleCom" = estimates$nMaleCom,
        "SampnFemNoCom" = estimates$nFemNoCom,
        "SampnMaleNoCom" = estimates$nMaleNoCom,
        "ModelSamp" = estimates$ncomplete,
        "MissingVar"= estimates$missing_var,
        "check" = data$Sample$SUBJID[1],
        "model_formula" = estimates$Formula,
        "SampVacCaseCom" = estimates$nVacCaseCom,
        "SampUnvacCaseCom" = estimates$nUnvacCaseCom,
        "SampVacCaseNoCom" = estimates$nVacCaseNoCom,
        "SampUnvacCaseNoCom" = estimates$nUnvacCaseNoCom,
        "SampVacNoncaseCom"= estimates$nVacNoncaseCom,
        "SampUnvacNoncaseCom" = estimates$nUnvacNoncaseCom,
        "SampVacNoncaseNoCom" = estimates$nVacNoncaseNoCom,
        "SampUnvacNoncaseNoCom"= estimates$nUnvacNoncaseNoCom,
        
        # List that can store larger, complex objects from each simulation
        ".complex" = list(
          "Table" = with(data$Sample, table( Y,A,FEMALE, RISKGR1, D)))
         # "Model" = estimates$Model,
         #"DataParams"= data$DataParams,
         # "Population" = data$Population,
         #  "Sample" = data$Sample
      ))
    })
  },
  main = {
    sim %<>% run()
  },
  last = {
    runtime<- sim %>% SimEngine::vars("total_runtime")
    print(paste0("Total run time is ", round(runtime/60), 
                 " minutes ", "(",round(runtime/60/60,2), " hours)"))
  },
  cluster_config = list(js = "slurm")
)

saveRDS(sim$results,
        file = paste0(file_path, "Output/", "Sim_Raw_", 
                      npop/1000, "K_", nsims, "_sims_", Sys.Date(), ".rds"))
saveRDS(sim$levels,
        file = paste0(file_path, "Output/", "Sim_Betas_", 
                      npop/1000, "K_", nsims, "_sims_", Sys.Date(), ".rds"))
saveRDS(sim$warnings,
        file = paste0(file_path, "Output/", "Sim_Warnings_", 
                      npop/1000, "K_", nsims, "_sims_", Sys.Date(), ".rds"))
saveRDS(sim$results_complex,
        file = paste0(file_path, "Output/","Sim_Complex_",  
                      npop/1000, "K_", nsims, "_sims_", Sys.Date(), ".rds"))

#sim$warnings
sim$errors

runtime <- sim %>% SimEngine::vars("total_runtime")
print(Sys.time())
print(paste0("This core started at ", t1," and ran for ",
             round(runtime/60), " minutes ", "(",round(runtime/60/60,2), " hours)"))
