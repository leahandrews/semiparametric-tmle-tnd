### Simulation Parameters for: "Targeted maximum likelihood 
# estimation of vaccine effectiveness and immune correlates 
# in test-negative design studies with missing data" 
# Code Created By: Leah I. B. Andrews
# Date: 04/10/26

library(readr)
library(dplyr)
library(sl3)




#### Generate Moderna COVE Toy Dataset for Plasmode Simulations ######

# simtest is dataset of 27,976 Moderna COVE per protocol individuals who 
# had race or ethnicity date recorded
# contains SUBJID, AGE, SEX, CALTIME (number of days until primary 
# COVID-19 endpoint or censoring from Sept 1, 2020)

# In our paper, we used data from blinded phase of Moderna COVE to generate 
# realistic covariate distributions for our simulations (i.e., plasmode simulations)
# Since we are unable to share Moderna COVE data, 
# we will generate a toy dataset of simtest in place of Moderna COVE data
# to demonstrate how we constructed our simulation study

### Toy Moderna COVE Dataset
set.seed(2024)
nsubj <- 27976
demo.df0 <- data.frame(SUBJID = paste0("ID", 1:nsubj),
                       AGE = round(rnorm(nsubj, mean = 50, sd = 10)),
                       SEX = sample(c("F", "M"), nsubj, replace = TRUE),

                       # Presence of Comorbidities
                       RISKGR1 = sample(c("At Risk", "Not At Risk"), nsubj,
                                       replace = TRUE, prob = c(1/4, 3/4)),

                       # Date of COVID-19 or Censoring in Blinded Phase
                       ADT = sample(seq(as.Date('2020/09/01'), as.Date('2021/03/31'), by="day"),
                                        nsubj, replace = TRUE))

### Adding Indicator Variables and Linear Splines to Toy COVE Dataset
simtest<- demo.df0 %>%
  mutate(
    # Creating Linear Calendar Splines
    CALTIME  = as.numeric(ADT - as.Date("2020-09-01")),
    CALTIMEs2 = I(CALTIME >= 90)*(CALTIME - 90),
    CALTIMEs3 = I(CALTIME >= 135)*(CALTIME - 135),
    CALTIMEs4 = I(CALTIME >= 180)*(CALTIME - 180),

    # Creating Indicators
    FEMALE = ifelse(SEX == "F", 1, 0),
    RISKGR1 = relevel(factor(RISKGR1), ref = "Not At Risk"),
    RISKGR1BIN = ifelse(RISKGR1 == "At Risk", 1, 0))
    #CALTIMEc = c(scale(CALTIME)))

### A, Y, W, and C are Generated According to the Following Logistic Regression Models #####

## Some Parameters Have Different Values Depending on the Confounding Setting:
# - Main Effects Setting: Sex, Comorbidity, Calendar Date Main Effects
# - Interaction Setting: Main Effects + Sex-Comorbidity Interaction
# - Splines Setting: Main Effects+ Sex-Comorbidity Interaction + Calendar Date Splines

####### High Immune Marker Level (A) Parameters #########

# Logistic Regression Formula to Generate A
a_formula <- "~ FEMALE + RISKGR1 + FEMALE:RISKGR1 + CALTIME + CALTIMEs2 + CALTIMEs3"

### Sex
sex_a_or <- 3 
sex_a_beta <- log(sex_a_or)
names(sex_a_beta) <- "FEMALE"

### Comorbidities
comorb_a_or <- 1/4 
comorb_a_beta <- log(comorb_a_or)
names(comorb_a_beta) <- "RISKGR1At Risk"

### Sex:Comorbidity Interaction

# Beta = 0 in Main Effects Setting
sexxcomorb_a_beta0 <- 0
names(sexxcomorb_a_beta0) <- "FEMALE:RISKGR1At Risk"

# Nonzero Beta in Interaction and Splines Setting
sexxcomorb_a_ror <-  4 
sexxcomorb_a_beta <- log(sexxcomorb_a_ror) # Log of Ratio of Two ORs  
names(sexxcomorb_a_beta) <- "FEMALE:RISKGR1At Risk"

### Calendar Time

# At Time = 0 (September 1, 2020)
sept_a_prob <- 0.3
sept_a_odds <- sept_a_prob/(1 - sept_a_prob)

# At Time = 90 (Beginning of December)
dec_a_prob <- 0.6 
dec_a_odds <- dec_a_prob/(1 - dec_a_prob)

# Main Effect Parameter
septdec_a_or <- dec_a_odds/sept_a_odds
septdec_a_beta <- log(septdec_a_or)/90

# At Time = 135 (Mid January)
jan_a_prob <- 0.7 
jan_a_odds <- jan_a_prob/(1 - jan_a_prob)
decjan_a_or <- jan_a_odds/dec_a_odds

# Beta = 0 in Main Effects and Interaction Setting
decjan_a_beta0 <- 0

# Nonzero Beta in Splines Setting
decjan_a_beta<- log(decjan_a_or)/45 - septdec_a_beta

# At Time 210 (End of March)
mar_a_prob <- 0.4 
mar_a_odds <- mar_a_prob/(1 - mar_a_prob)
janmar_a_or <- mar_a_odds/jan_a_odds

# Beta = 0 in Main Effects and Interaction Setting
janmar_a_beta0 <- 0

# Nonzero Beta in Splines Setting
janmar_a_beta<- log(janmar_a_or)/75 - decjan_a_beta-septdec_a_beta

# Calendar Time Betas in Main Effects and Interaction Settings
time_a_betas0<- c(septdec_a_beta, decjan_a_beta0, janmar_a_beta0)
names(time_a_betas0) <- c("CALTIME","CALTIMEs2","CALTIMEs3")

# Calendar Time Betas in Splines Setting
time_a_betas<- c(septdec_a_beta, decjan_a_beta,janmar_a_beta)
names(time_a_betas) <- c("CALTIME","CALTIMEs2","CALTIMEs3")

### Intercept (Odds of High Immune Marker Level 
# for Male w/o Comorbidities in September)
ref_a_odds <- .25/(1 - .25) 
names(ref_a_odds) <-"(Intercept)" 

### Main Effects Setting Parameters: Interaction and Splines Parameters = 0
# Note: parameters can be in any order, but must be named
a_betas1 <- c(log(ref_a_odds), time_a_betas0,
               comorb_a_beta, sex_a_beta, sexxcomorb_a_beta0)

### Interaction Setting Parameters: Splines Parameters = 0  
# Note: parameters can be in any order, but must be named
a_betas2 <- c(log(ref_a_odds), time_a_betas0,
               comorb_a_beta, sex_a_beta, sexxcomorb_a_beta)

### Splines Setting Parameters
# Note: parameters can be in any order, but must be named
a_betas3 <- c(log(ref_a_odds), time_a_betas,
               comorb_a_beta, sex_a_beta, sexxcomorb_a_beta)

######## SARS-CoV-2 Infection (Y) Parameters ########

# Logistic Regression Formula to Generate Y
y_formula <- "~ A + FEMALE + RISKGR1 + FEMALE:RISKGR1 + CALTIME + CALTIMEs2 + CALTIMEs3"

### Sex 
sex_y_or <- 3 
names(sex_y_or) <- "FEMALE"
sex_y_beta<- log(sex_y_or)

### Comorbidities
comorb_y_or <- 4 
names(comorb_y_or) <-"RISKGR1At Risk"
comorb_y_beta <- log(comorb_y_or)

### Sex:Comorbidity Interaction

# Beta = 0 in Main Effects Setting
sexxcomorb_y_beta0 <- 0
names(sexxcomorb_y_beta0) <-"FEMALE:RISKGR1At Risk"

# Nonzero Beta in Interaction and Splines Setting
sexxcomorb_y_ror <- 0.25 
sexxcomorb_y_beta <-  log(sexxcomorb_y_ror)  # Log of Ratio of Two ORs  
names(sexxcomorb_y_beta) <-"FEMALE:RISKGR1At Risk"

### Calendar Time

# At Time = 0 (September 1, 2020)
sept_y_prob <- 0.008 
sept_y_odds <- sept_y_prob/(1 - sept_y_prob)

# At Time = 90 (Beginning of December)
dec_y_prob <- 0.003 
dec_y_odds <- dec_y_prob/(1 - dec_y_prob)

# Main Effect Parameter
septdec_y_or <- dec_y_odds/sept_y_odds
septdec_y_beta <- log(septdec_y_or)/90

# At Time = 135 (Mid January)
jan_y_prob <- 0.0015 
jan_y_odds <- jan_y_prob/(1 - jan_y_prob)
decjan_y_or <- jan_y_odds/dec_y_odds

# Beta = 0 in Main Effects and Interaction Setting
decjan_y_beta0 <- 0

# Nonzero Beta in Splines Setting
decjan_y_beta<- log(decjan_y_or)/45 - septdec_y_beta

# At Time 210 (End of March)
mar_y_prob <- 0.006 
mar_y_odds <- mar_y_prob/(1 - mar_y_prob) 
janmar_y_or <- mar_y_odds/jan_y_odds

# Beta = 0 in Main Effects and Interaction Setting
janmar_y_beta0 <- 0

# Nonzero Beta in Splines Setting
janmar_y_beta <- log(janmar_y_or)/75 - decjan_y_beta - septdec_y_beta

# Calendar Time Betas in Main Effects and Interaction Settings
time_y_betas0 <- c(septdec_y_beta, decjan_y_beta0, janmar_y_beta0)
names(time_y_betas0) <- c("CALTIME","CALTIMEs2","CALTIMEs3")

# Calendar Time Betas in Splines Setting
time_y_betas <- c(septdec_y_beta, decjan_y_beta, janmar_y_beta)
names(time_y_betas) <- c("CALTIME","CALTIMEs2","CALTIMEs3")

### Intercept (Odds of SARS-CoV-2 Infection 
# for Males w/o Comorbidities in September)
ref_y_prob <- 0.13 
ref_y_odds <- ref_y_prob/(1 - ref_y_prob)  
names(ref_y_odds) <-"(Intercept)" 

### Main Effects Setting Parameters: Interaction and Splines Parameters = 0
# Note: parameters can be in any order, but must be named
# Note: log OR for A is omitted in y_betas1, but will be included as
# a separate input in data-generating function in methods.R
y_betas1 <- c(log(ref_y_odds), time_y_betas0, 
               comorb_y_beta, sex_y_beta, sexxcomorb_y_beta0)

### Interaction Setting Parameters: Spline Parameters = 0
# Note: parameters can be in any order, but must be named
# Note: log OR for A is omitted in y_beta2, but will be included as
# a separate input in data-generating function in methods.R
y_betas2 <- c(log(ref_y_odds), time_y_betas0, 
               comorb_y_beta, sex_y_beta, sexxcomorb_y_beta)

### Splines Setting Parameters 
# Note: parameters can be in any order, but must be named
# Note: log OR for A is omitted in y_betas3, but will be included as
# a separate input in data-generating function in methods.R
y_betas3 <- c(log(ref_y_odds), time_y_betas, 
               comorb_y_beta, sex_y_beta, sexxcomorb_y_beta)

##### Infections Other than SARS-CoV-2 that Could Induce Symptoms (W) Parameters ########

# Logistic Regression Formula to Generate W
w_formula <- "~ A + FEMALE + RISKGR1+ CALTIME + CALTIMEs2 + CALTIMEs4"

### Sex
sex_w_or <- 2 
names(sex_w_or) <- "FEMALE"
sex_w_beta <- log(sex_w_or)

### Comorbidities 
comorb_w_or <- 2
names(comorb_w_or) <- "RISKGR1At Risk"
comorb_w_beta <- log(comorb_w_or)

### Calendar Time

# At Time = 0 (September 1, 2020)
sept_w_prob <- 0.07/0.3 
sept_w_odds <- sept_w_prob/(1 - sept_w_prob)

# At Time = 90 (Beginning of December)
dec_w_prob <- 0.14/0.3 
dec_w_odds <- dec_w_prob/(1 - dec_w_prob)
septdec_w_or <- dec_w_odds/sept_w_odds
septdec_w_beta<- log(septdec_w_or)/90

# At Time = 180 (Beginning of March, Different Cutpoint than for Y)
mar_w_prob <- 0.17/0.3
mar_w_odds <- mar_w_prob/(1 - mar_w_prob)
decmar_w_or <- mar_w_odds/dec_w_odds
decmar_w_beta<- log(decmar_w_or)/90 - septdec_w_beta

# At Time = 300 (End of June)
jun_w_prob <- 0.08/0.3
jun_w_odds <- jun_w_prob/(1 - jun_w_prob)
marjun_w_or <- jun_w_odds/mar_w_odds
marjun_w_beta <- log(marjun_w_or)/120 - decmar_w_beta - septdec_w_beta

# Calendar Time Betas for All Settings
time_w_betas<- c(septdec_w_beta, decmar_w_beta, marjun_w_beta)
names(time_w_betas) <- c("CALTIME","CALTIMEs2","CALTIMEs4")

### Immune Marker Level A
imm_w_or <- 1 # Assume Immune Marker Level Doesn't Affect W
imm_w_beta <- log(imm_w_or)
names(imm_w_beta) <- "A" 

### Odds of Non-SARS Infection for Reference Group
ref_w_odds <- 0.1 
names(ref_w_odds) <-"(Intercept)" 

# Parameters for W 
# Note: parameters can be in any order, must be named
# Note: imm_w_beta is omitted in w_betas, but will be included as
# a separate input in data-generating function in methods.R
w_betas <- c(log(ref_w_odds), time_w_betas, 
             comorb_w_beta, sex_w_beta)

#### Meeting the Symptom Definition C #####
# Note: In our manuscript, D represents meeting the symptom definition
# In our code, C represents meeting the symptom definition 
# because we use D as shorthand for Delta (observing A)

# Logistic Regression Formula to Generate C
c_formula <- "~Y + W + A + Y*FEMALE + Y*RISKGR1 + W*FEMALE + W*RISKGR1"

### Intercept (Odds of Having Symptoms when Male, No Comorbidities, 
# Low Immune Marker Level, No Infections)
ref_c_odds <- 0.1 
ref_c_beta <- log(ref_c_odds)
names(ref_c_beta) <- "(Intercept)"

### Immune Marker Level
imm_c_or <- 1 # Assume Immune Marker Level Doesn't Affect Symptoms
imm_c_beta <- log(imm_c_or)
names(imm_c_beta) <- "A" 

### No Infection 

# Sex Main Effect
noinf_sex_c_beta <- 0
names(noinf_sex_c_beta) <- "FEMALE"

# Comorbidities Main Effect
noinf_comorb_c_or <- 2 
names(noinf_comorb_c_or) <- "RISKGR1At Risk"
noinf_comorb_c_beta <- log(noinf_comorb_c_or)

### Other Infections W

# W Main Effect (log OR of C Comparing W Status in Males without Comorbidities)
other_c_or <- 4 
other_c_beta <- log(other_c_or)
names(other_c_beta) <- "W"

# W:Sex Interaction
other_fem_c_prob <- 0.96
other_fem_c_odds <- other_fem_c_prob/(1 - other_fem_c_prob)
other_male_c_prob <- 0.80
other_male_c_odds <- other_male_c_prob/(1 - other_male_c_prob)
other_sex_c_or <-other_fem_c_odds/other_male_c_odds
other_sex_c_beta<- log(other_sex_c_or) - noinf_sex_c_beta
names(other_sex_c_beta) <- "W:FEMALE"

# W:Comorbidities Interaction
other_comorb_c_or <- 1.06 
names(other_comorb_c_or) <- "W:RISKGR1At Risk"
other_comorb_c_beta<- log(other_comorb_c_or) - noinf_comorb_c_beta

### SARS-CoV-2 Infection 

# Y Main Effects Term (log OR of C Comparing Y Status in Males without Comorbidities)
sars_c_prob <- 0.6
sars_c_odds <- sars_c_prob/(1 - sars_c_prob)
noinf_c_prob <- 0.1
noinf_c_odds <- noinf_c_prob/(1-noinf_c_prob)
sars_c_or <- sars_c_odds/noinf_c_odds
sars_c_beta <- log(sars_c_or)
names(sars_c_beta) <- "Y"

# Y:Sex Interaction
sars_sex_c_or <- 1/0.93
names(sars_sex_c_or) <- "Y:FEMALE"
sars_sex_c_beta <- log(sars_sex_c_or) - noinf_sex_c_beta

# Y:Comorbidities Interaction
sars_comorb_c_or <- 1.06
names(sars_comorb_c_or) <- "Y:RISKGR1At Risk"
sars_comorb_c_beta <- log(sars_comorb_c_or) - noinf_comorb_c_beta

# Parameters for C 
# Note: parameters can be in any order, must be named
# Note: imm_c_beta is omitted in w_betas, but will be included as a 
# separate input in data-generating function in methods.R
c_betas <- c(ref_c_beta, 
             noinf_sex_c_beta, noinf_comorb_c_beta, 
             other_c_beta,  other_sex_c_beta,other_comorb_c_beta, 
             sars_c_beta, sars_sex_c_beta, sars_comorb_c_beta)