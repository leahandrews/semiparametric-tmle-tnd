## Simulation figures and tables for: 
# "Targeted maximum likelihood estimation of vaccine effectiveness 
# and immune correlates in test-negative design studies with missing data"
# Note: exact results will differ from manuscript because manuscript uses
# Moderna COVE data for covariate distributions, not the toy dataset 
# simulated in params.R
# Code Created By: Leah I. B. Andrews
# Date: 04/10/26


library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(flextable)

### Load saved simulation data
raw.df0<- readRDS("~/Desktop/UW Stuff/IS TND/Code/TND TMLE/Output/073125/Sim_Raw_50K_1000_sims_2025-07-31.rds")
warn.df0<- readRDS("~/Desktop/UW Stuff/IS TND/Code/TND TMLE/Output/073125/Sim_Warnings_50K_1000_sims_2025-07-31.rds")
betas.df<- readRDS("~/Desktop/UW Stuff/IS TND/Code/TND TMLE/Output/073125/Sim_Betas_50K_1000_sims_2025-07-31.rds")


# Modify variable names for figures and tables
raw.df <- raw.df0 %>% 
  filter(Estimator %in% c("OrdLogitRev","OrdLogitRevx",
                          "PLLogitEx","PLLogitMx",
                          "PLLogitE0","PLLogitM0",
                          "SemiLogit0"))%>%
  mutate(
    Estimator0 = Estimator,
    Estimator = factor(case_when(Estimator0 == "OrdLogitRev" ~ "nMLE",
                                 Estimator0 == "OrdLogitRevx" ~ "MLEx",
                                 Estimator0 == "PLLogitE0" ~ "nPLE",
                                 Estimator0 == "PLLogitM0" ~ "nPLM",
                                 Estimator0 == "PLLogitEx" ~ "PLEx",
                                 Estimator0 == "PLLogitMx" ~ "PLMx",
                                 Estimator0 == "SemiLogit0" ~ "TMLE"),
                       levels = c("TMLE","PLMx","PLEx","MLEx",
                                  "nPLM","nPLE", "nMLE")),
    samprate0 = samprate,
    samprate = factor(case_when( 
      samprate0 == "100" ~ "All",  #"Complete Data"
      samprate0 == "3" ~ "1:3",
      samprate0 == "1" ~ "1:1")),
    sampratemath = case_match(samprate, "All" ~ "'All'",
                              "1:3" ~ "'1:3'",
                              "1:1" ~ "'1:1'",
                              .ptype = factor(levels = c("'1:1'", "'1:3'", 
                                                         "'All'"))),
    truelogORc = factor(paste0("log(",round(exp(truelogOR),2),")"), 
                        levels=c("log(1)","log(0.7)","log(0.2)")),
    truelogORmath = case_match(truelogORc, "log(1)" ~ "beta(P[F]) == log(1)", 
                               "log(0.7)" ~ "beta(P[F]) == log(0.7)",
                               "log(0.2)" ~ "beta(P[F]) == log(0.2)", 
                               .ptype = factor(levels = c("beta(P[F]) == log(1)",
                                                          "beta(P[F]) == log(0.7)",
                                                          "beta(P[F]) == log(0.2)"))),
    distributionc = factor(case_when(distribution == "Scenario 1" ~ "Main Effects",
                                     distribution == "Scenario 2" ~ "Interaction",
                                     distribution == "Scenario 3" ~ "Splines"),
                           levels = c("Main Effects", "Interaction", "Splines")),
    distributionmath = case_match(distributionc, "Main Effects" ~ "'Main Effects'",
                                  "Interaction" ~ "'Interaction'",
                                  "Splines" ~ "'Splines'", 
                                  .ptype = factor(levels = c("'Main Effects'",
                                                             "'Interaction'",
                                                             "'Splines'")))
  )



# Summarize simulation results by simulation setting and estimator
summ.df <- raw.df %>% 
  dplyr:: group_by(trueVE, distribution, nsamp, samprate, Estimator) %>%  
  dplyr::summarise(
    nsamp = first(nsamp),
    nsims = n(),
    truelogOR = first(truelogOR),
    truelogORc = first(truelogORc),
    truelogORmath = first(truelogORmath),
    distributionc = first(distributionc),
    distributionmath = first(distributionmath),
    sampratemath = first(sampratemath),
    ModelSamp = mean(ModelSamp),
    SampObsMean = mean(SampObsProp, na.rm=TRUE),
    logORMean = mean(logOR, na.rm = TRUE), 
    logORBias = mean(logOR-truelogOR, na.rm = TRUE), 
    VEMean = mean(VE, na.rm = TRUE), 
    VEBias = mean(VE - trueVE, na.rm = TRUE), 
    logORMCSD = sd(logOR, na.rm = TRUE), 
    logORMCVar = var(logOR, na.rm = TRUE), 
    logORSEMean = mean(logORSE, na.rm = TRUE), 
    logORVarMean = mean(logORVar, na.rm = TRUE),  
    VECoverage = mean(trueVE > VELL & trueVE < VEUL , na.rm = TRUE), 
    Power = mean(Reject, na.rm = TRUE), 
    logORMSE = mean((logOR-truelogOR)^ 2, na.rm = TRUE),
    PopMargAMean = mean(PopMargA, na.rm = TRUE),
    PopMargYMean = mean(PopMargY, na.rm = TRUE),
    PopMargWMean = mean(PopMargW, na.rm = TRUE),
    PopMargCMean = mean(PopMargC, na.rm = TRUE),
    Phase1MargAMean = mean(Phase1MargA, na.rm = TRUE),
    Phase1MargYMean = mean(Phase1MargY, na.rm = TRUE),
    Phase1MargWMean = mean(Phase1MargW, na.rm = TRUE),
    Phase1MargCMean = mean(Phase1MargC, na.rm = TRUE),
    Phase1MargFemMean = mean(Phase1MargFem, na.rm = TRUE),
    Phase1MargComMean = mean(Phase1MargCom, na.rm = TRUE),
    
    Phase2MargAMean = mean(Phase2MargA, na.rm = TRUE),
    Phase2MargYMean = mean(Phase2MargY, na.rm = TRUE),
    Phase2MargnAMean = mean(Phase2MargnA, na.rm = TRUE),
    Phase2MargnYMean = mean(Phase2MargnY, na.rm = TRUE),
    Phase2MargWMean = mean(Phase2MargW, na.rm = TRUE),
    Phase2MargCMean = mean(Phase2MargC, na.rm = TRUE),
    Phase2MargFemMean = mean(SampFemObs, na.rm = TRUE),
    Phase2MargComMean = mean(SampComObs, na.rm = TRUE),
    SampVacObsMean = mean(SampVacObs, na.rm = TRUE),
    SampUnvacObsMean = mean(SampUnvacObs, na.rm = TRUE),
    ObsSampRatioMean = mean(ObsSampRatio, na.rm = TRUE),
    SampObsPropMean = mean(SampObsProp, na.rm = TRUE)) 

group.colors <- c(MLEx ="#D55E00", 
                  nMLE =  "#E69F00", 
                  PLEx = "#0072B2",
                  PLMx = "#882255",
                  nPLE = "#56B4E9",
                  nPLM = "#CC79A7", 
                  TMLE =  "#009E73") 

group.linetype <- c(MLEx = "dotdash" , 
                    nMLE= "longdash", 
                    PLEx =  "dotted", 
                    PLMx =  "twodash", 
                    nPLE = "12345678",
                    nPLM = "4C88C488",
                    TMLE =  "solid")


############# Data-Generating Mechanism Simulation Plots ###################

### Generating Immune Marker Level A Under Main Effects Setting
(a1 <- data.frame(FEM = c(rep(0,211), rep(0,211), rep(1,211), rep(1,211)),
                  COM = c(rep(0,211), rep(1,211), rep(0,211), rep(1,211)),
                  TIME = rep(c(0:210), 4)) %>% mutate(
                  y=c(plogis(
                     betas.df$distribution$`Scenario 1`$a_betas["(Intercept)"] +
                     betas.df$distribution$`Scenario 1`$a_betas["FEMALE"]*FEM +
                     betas.df$distribution$`Scenario 1`$a_betas["RISKGR1At Risk"]*COM +
                     betas.df$distribution$`Scenario 1`$a_betas["FEMALE:RISKGR1At Risk"]*FEM*COM +
                     betas.df$distribution$`Scenario 1`$a_betas["CALTIME"]*TIME +
                     betas.df$distribution$`Scenario 1`$a_betas["CALTIMEs2"]*I(TIME >= 90)*(TIME - 90) +
                     betas.df$distribution$`Scenario 1`$a_betas["CALTIMEs3"]*I(TIME >= 135)*(TIME - 135))),
                  Sex = ifelse(FEM == 1, "Female", "Male"),
                  Comorbidity = ifelse(COM == 1, "Comorbidity", "No Comorbidity")) %>%
    ggplot(aes(x = TIME,y = y, col = Sex, linetype = Comorbidity)) + ylim(0, 1) +
    scale_color_manual(values = c(Male = "#1E88E5", Female = "#D81B60"))+
    geom_line() + theme_bw() + scale_linetype_manual( values = c(`No Comorbidity` = "longdash",
                                                             Comorbidity = "dotted"))+
    labs(x = "Time", y = "High Immune Marker\nLevel Probability", 
         title = "Main Effects: Immune Marker") +
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "bottom") )

### Generating Immune Marker Level A Under Interaction Setting
(a2 <- data.frame(FEM=c(rep(0,211), rep(0,211), rep(1,211), rep(1,211)),
                  COM =c(rep(0,211), rep(1,211), rep(0,211), rep(1,211)),
                  TIME = rep(c(0:210), 4)) %>% mutate(
                  y= c(plogis(
                      betas.df$distribution$`Scenario 2`$a_betas["(Intercept)"] +
                      betas.df$distribution$`Scenario 2`$a_betas["FEMALE"]*FEM +
                      betas.df$distribution$`Scenario 2`$a_betas["RISKGR1At Risk"]*COM +
                      betas.df$distribution$`Scenario 2`$a_betas["FEMALE:RISKGR1At Risk"]*FEM*COM +
                      betas.df$distribution$`Scenario 2`$a_betas["CALTIME"]*TIME +
                      betas.df$distribution$`Scenario 2`$a_betas["CALTIMEs2"]*I(TIME >= 90)*(TIME - 90) +
                      betas.df$distribution$`Scenario 2`$a_betas["CALTIMEs3"]*I(TIME >= 135)*(TIME - 135))),
                  Sex = ifelse(FEM == 1, "Female", "Male"),
                  Comorbidity = ifelse(COM == 1, "Comorbidity", "No Comorbidity")) %>%
    ggplot(aes(x = TIME,y = y, col = Sex, linetype = Comorbidity)) + ylim(0, 1) +
    scale_color_manual(values = c(Male = "#1E88E5", Female = "#D81B60")) +
    geom_line() + theme_bw() + 
    scale_linetype_manual( values =  c(`No Comorbidity` = "longdash", Comorbidity = "dotted"))+
    labs(x="Time", y="", title = "Interaction: Immune Marker")+
    theme(axis.title.x = element_text( size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "bottom") )

### Generating Immune Marker Level A Under Splines Setting
(a3 <- data.frame(FEM = c(rep(0, 211), rep(0, 211), rep(1, 211), rep(1, 211)),
                  COM = c(rep(0, 211), rep(1, 211), rep(0, 211), rep(1, 211)),
                  TIME = rep(c(0:210), 4)) %>% mutate(
                  y=c(plogis(betas.df$distribution$`Scenario 3`$a_betas["(Intercept)"] +
                     betas.df$distribution$`Scenario 3`$a_betas["FEMALE"]*FEM +
                     betas.df$distribution$`Scenario 3`$a_betas["RISKGR1At Risk"]*COM +
                     betas.df$distribution$`Scenario 3`$a_betas["FEMALE:RISKGR1At Risk"]*FEM*COM +
                     betas.df$distribution$`Scenario 3`$a_betas["CALTIME"]*TIME +
                     betas.df$distribution$`Scenario 3`$a_betas["CALTIMEs2"]*I(TIME >= 90)*(TIME - 90) +
                     betas.df$distribution$`Scenario 3`$a_betas["CALTIMEs3"]*I(TIME >= 135)*(TIME - 135))),
                  Sex = ifelse(FEM == 1, "Female","Male"),
                  Comorbidity = ifelse(COM == 1, "Comorbidity", "No Comorbidity")) %>%
    ggplot(aes(x = TIME, y = y, col = Sex, linetype = Comorbidity)) +  ylim(0, 1) +
    scale_color_manual(values = c(Male = "#1E88E5", Female = "#D81B60")) +
    geom_line() + theme_bw() + scale_linetype_manual( values = c(`No Comorbidity` = "longdash",
                                                             Comorbidity = "dotted"))+
    labs(x = "Time", y = "", title = "Splines: Immune Marker")+
    theme(axis.title.x = element_text( size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "bottom"))

### Generating SARS-CoV-2 Infection Y Under Main Effects Setting
(y1 <- data.frame(FEM = c(rep(0, 211), rep(0, 211), rep(1, 211), rep(1, 211)),
                  COM = c(rep(0, 211), rep(1, 211), rep(0, 211), rep(1, 211)),
                  TIME = rep(c(0:210), 4)) %>% mutate(
                  y= c(plogis(betas.df$distribution$`Scenario 1`$y_betas["(Intercept)"] +
                      betas.df$distribution$`Scenario 1`$y_betas["FEMALE"]*FEM +
                      betas.df$distribution$`Scenario 1`$y_betas["RISKGR1At Risk"]*COM +
                      betas.df$distribution$`Scenario 1`$y_betas["FEMALE:RISKGR1At Risk"]*FEM*COM +
                      betas.df$distribution$`Scenario 1`$y_betas["CALTIME"]*TIME +
                      betas.df$distribution$`Scenario 1`$y_betas["CALTIMEs2"]*I(TIME >= 90)*(TIME - 90) +
                      betas.df$distribution$`Scenario 1`$y_betas["CALTIMEs3"]*I(TIME >= 135)*(TIME - 135))),
                  Sex = ifelse(FEM == 1, "Female","Male"),
                  Comorbidity = ifelse(COM == 1, "Comorbidity", "No Comorbidity")) %>%
    ggplot(aes(x = TIME,y = y, col = Sex, linetype = Comorbidity))+  ylim(0,1) +
    scale_color_manual(values = c(Male = "#1E88E5", Female = "#D81B60")) +
    geom_line() + theme_bw() + 
    scale_linetype_manual( values = c(`No Comorbidity` = "longdash", Comorbidity = "dotted"))+
    labs(x = "Time", y = "SARS-CoV-2\nInfection Probability", 
         title = "Main Effects: SARS-CoV-2") +
    theme(axis.title.x = element_text( size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "bottom"))

### Generating SARS-CoV-2 Infection Y Under Interaction Setting
(y2 <- data.frame(FEM=c(rep(0, 211), rep(0, 211), rep(1, 211), rep(1, 211)),
                  COM =c(rep(0, 211), rep(1, 211), rep(0, 211), rep(1, 211)),
                  TIME = rep(c(0:210), 4)) %>% mutate(
                  y= c(plogis(betas.df$distribution$`Scenario 2`$y_betas["(Intercept)"] +
                      betas.df$distribution$`Scenario 2`$y_betas["FEMALE"]*FEM +
                      betas.df$distribution$`Scenario 2`$y_betas["RISKGR1At Risk"]*COM +
                      betas.df$distribution$`Scenario 2`$y_betas["FEMALE:RISKGR1At Risk"]*FEM*COM +
                      betas.df$distribution$`Scenario 2`$y_betas["CALTIME"]*TIME +
                      betas.df$distribution$`Scenario 2`$y_betas["CALTIMEs2"]*I(TIME >= 90)*(TIME - 90)+
                      betas.df$distribution$`Scenario 2`$y_betas["CALTIMEs3"]*I(TIME >= 135)*(TIME - 135))),
                  Sex = ifelse(FEM == 1, "Female","Male"),
                  Comorbidity = ifelse(COM == 1, "Comorbidity", "No Comorbidity")) %>%
    ggplot(aes(x = TIME,y = y, col = Sex, linetype=Comorbidity)) +  ylim(0,1) +
    scale_color_manual(values = c(Male = "#1E88E5", Female = "#D81B60"))+
    geom_line() + theme_bw() + 
    scale_linetype_manual( values = c(`No Comorbidity` = "longdash", Comorbidity = "dotted"))+
    labs(x = "Time", y = "",  title = "Interaction: SARS-CoV-2") +
    theme(axis.title.x = element_text( size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "bottom") )

### Generating SARS-CoV-2 Infection Y Under Splines Setting
(y3 <- data.frame(FEM = c(rep(0, 211), rep(0, 211), rep(1, 211), rep(1, 211)),
                  COM = c(rep(0, 211), rep(1, 211), rep(0, 211), rep(1, 211)),
                  TIME = rep(c(0:210), 4)) %>% mutate(
                  y = c(plogis(betas.df$distribution$`Scenario 3`$y_betas["(Intercept)"] +
                      betas.df$distribution$`Scenario 3`$y_betas["FEMALE"]*FEM +
                      betas.df$distribution$`Scenario 3`$y_betas["RISKGR1At Risk"]*COM +
                      betas.df$distribution$`Scenario 3`$y_betas["FEMALE:RISKGR1At Risk"]*FEM*COM +
                      betas.df$distribution$`Scenario 3`$y_betas["CALTIME"]*TIME +
                      betas.df$distribution$`Scenario 3`$y_betas["CALTIMEs2"]*I(TIME >= 90)*(TIME - 90)+
                      betas.df$distribution$`Scenario 3`$y_betas["CALTIMEs3"]*I(TIME >= 135)*(TIME - 135))),
                  Sex = ifelse(FEM == 1, "Female", "Male"),
                  Comorbidity = ifelse(COM == 1, "Comorbidity", "No Comorbidity")) %>%
    ggplot(aes(x = TIME, y = y, col = Sex, linetype = Comorbidity)) + ylim(0,1) +
    scale_color_manual(values = c(Male = "#1E88E5", Female = "#D81B60"))+
    geom_line() + theme_bw() + 
    scale_linetype_manual(values = c(`No Comorbidity` = "longdash", Comorbidity = "dotted"))+
    labs(x = "Time", y = "", title = "Splines: SARS-CoV-2")+
    theme(axis.title.x = element_text( size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "bottom"))

### Figure S1
pdf(file = "~/Desktop/UW Stuff/IS TND/Code/TND TMLE/semiparametric-tmle-tnd/Simulation/Output/AY Data Gen 040926.pdf",
    height = 5, width = 10.5)
ggarrange( a1, a2, a3, y1, y2, y3,  nrow = 2, ncol = 3,
           common.legend = TRUE, legend = "bottom" , 
          labels = c("A1", "A2", "A3", "B1", "B2","B3") )
dev.off()

########## Simulation Result: Bias Figures #################

### Figure S2
pdf(file = "~/Desktop/UW Stuff/IS TND/Code/TND TMLE/semiparametric-tmle-tnd/Simulation/Output/Bias Main Effects.pdf",
    height = 6.5, width = 10)
summ.df %>% filter(distribution == "Scenario 1") %>% 
  mutate(trueVEdec = trueVE/100) %>%
  ggplot( aes(x = nsamp, y = logORBias, col = Estimator, linetype = Estimator)) +
  geom_line(size = 0.7) +
  geom_hline(yintercept = 0, color = "black", lty = 2, size = .5) + theme_bw() +
  facet_grid( truelogORmath ~ distributionmath + sampratemath, 
              label = "label_parsed")+
  labs(y = expression(paste("Log OR(",P[F],") Bias")), x = "Sample Size") +
  theme(legend.position = "bottom")+ guides(color = guide_legend(nrow = 1)) +
  scale_color_manual(values = group.colors) +
  scale_linetype_manual(values = group.linetype) +
  coord_cartesian(
    xlim = NULL,
    ylim = c(-.41, .15),
    expand = TRUE,
    default = FALSE,
    clip = "on"
  ) + theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          strip.text = element_text( face = "bold", size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.position = "bottom") 
dev.off()

### Figure 2
pdf(file = "~/Desktop/UW Stuff/IS TND/Code/TND TMLE/semiparametric-tmle-tnd/Simulation/Output/Bias Inx Splines.pdf",
    height = 6.5, width = 10)
summ.df %>% filter(distribution != "Scenario 1") %>% 
  mutate(trueVEdec = trueVE/100) %>%   
  ggplot( aes(x = nsamp, y = logORBias, col = Estimator, linetype = Estimator)) +
  geom_line(size = 0.7) +
  geom_hline(yintercept = 0, color = "black", lty = 2, size = .5) + theme_bw() +
  facet_grid( truelogORmath ~ distributionmath + sampratemath, 
              label = "label_parsed" ) +
 labs(y = expression(paste("Log OR(",P[F],") Bias")), x = "Sample Size") +
  theme(legend.position = "bottom") + guides(color = guide_legend(nrow = 1)) +
  scale_color_manual(values = group.colors) +
  scale_linetype_manual(values = group.linetype)+ 
  coord_cartesian(
    xlim = NULL,
    ylim = c(-.41, .15), 
    expand = TRUE,
    default = FALSE,
    clip = "on"
  )+theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          strip.text = element_text( face = "bold", size = 14), 
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.position = "bottom") 
dev.off()


########## Simulation Result: Coverage Figures ##############

### Figure S3
pdf(file = "~/Desktop/UW Stuff/IS TND/Code/TND TMLE/semiparametric-tmle-tnd/Simulation/Output/Coverage Main Effects 040926.pdf",
    height=6.5, width=10) 
summ.df%>% filter(distributionc == "Main Effects")%>%
  mutate(trueVEdec = trueVE/100) %>%  
  ggplot( aes(x = nsamp, y = VECoverage, col = Estimator, 
              linetype=Estimator)) +
  geom_line( size=.7)+
  geom_hline(yintercept = 0.95, color = "black", lty = 2)+ theme_bw() +
  facet_grid( truelogORmath~distributionmath+sampratemath, 
              label= "label_parsed" )+
  labs(y="95% Confidence Interval Coverage",x= "Sample Size")+ 
  theme(legend.position = "bottom")+ guides(color = guide_legend(nrow = 1))+
  scale_color_manual(values=group.colors)+
  scale_linetype_manual(values= group.linetype)+
  coord_cartesian(
    xlim = NULL,
    ylim = c(.63,.97), 
    expand = TRUE,
    default = FALSE,
    clip = "on"
  )+theme(axis.title.x = element_text( size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size = 14),
          strip.text = element_text( face = "bold", size=14), 
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.position = "bottom") 
dev.off()

### Figure 3
pdf(file = "~/Desktop/UW Stuff/IS TND/Code/TND TMLE/semiparametric-tmle-tnd/Simulation/Output/Coverage Inx and Splines 040926.pdf",
    height = 6.5, width = 10) 
summ.df %>% filter(distributionc != "Main Effects") %>%
  mutate(trueVEdec = trueVE/100) %>%  
  ggplot( aes(x = nsamp, y = VECoverage, col = Estimator, 
              linetype = Estimator)) +
  geom_line( size = .7) +
  geom_hline(yintercept = 0.95, color = "black", lty = 2) + theme_bw() +
  facet_grid(truelogORmath ~ distributionmath + sampratemath, 
              label = "label_parsed" )+
  labs(y = "95% Confidence Interval Coverage", x = "Sample Size")+ 
  theme(legend.position = "bottom")+
  scale_color_manual(values = group.colors) + guides(color = guide_legend(nrow = 1))+
  scale_linetype_manual(values = group.linetype) +
  coord_cartesian(
    xlim = NULL,
    ylim = c(.63, .97), 
    expand = TRUE,
    default = FALSE,
    clip = "on"
  )+theme(axis.title.x = element_text( size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          strip.text = element_text( face = "bold", size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.position = "bottom") 
dev.off()

########## Simulation Result: Type 1 Error and Power Figures ##############

### Figure S4
pdf(file = "~/Desktop/UW Stuff/IS TND/Code/TND TMLE/semiparametric-tmle-tnd/Simulation/Output/Type 1 Error.pdf",
    height = 6, width = 10)
summ.df %>% filter(trueVE == 0 ) %>% 
  mutate(trueVEdec = trueVE/100) %>%
  ggplot( aes(x = nsamp, y = Power, col = Estimator, linetype = Estimator)) +
  geom_line() +
  geom_hline(yintercept = 0.05, color = "black", lty = 2) + theme_bw() +
  facet_grid(distributionc ~ samprate) +
  labs(y = "Type 1 Error", x = "Sample Size") + 
  theme(legend.position = "bottom") + guides(color = guide_legend(nrow = 1))+
  scale_color_manual(values = group.colors) +
  scale_linetype_manual( values = group.linetype)+
  coord_cartesian(
    xlim = NULL,
    ylim = c(0, .4), 
    expand = TRUE,
    default = FALSE,
    clip = "on")+ 
  theme(axis.title.x = element_text( size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text( face = "bold", size = 14), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = "bottom") 
dev.off()

### Figure S5
pdf(file = "~/Desktop/UW Stuff/IS TND/Code/TND TMLE/semiparametric-tmle-tnd/Simulation/Output/Power.pdf",
    height = 10.5, width = 10) 
summ.df%>% filter(trueVE != 0 ) %>% 
  mutate(trueVEdec = trueVE/100) %>%  
  ggplot( aes(x = nsamp, y = Power, col = Estimator, linetype = Estimator)) +
  geom_line()+
  theme_bw() +
  facet_grid(distributionmath + truelogORmath ~ sampratemath, 
             label = "label_parsed")+
  labs(y = "Power",x = "Sample Size") + 
  theme(legend.position = "bottom")+ guides(color = guide_legend(nrow = 1))+
  scale_color_manual(values = group.colors)+
  scale_linetype_manual(values = group.linetype)+
  coord_cartesian(
    xlim = NULL,
    expand = TRUE,
    default = FALSE,
    clip = "on"
  ) + theme(axis.title.x = element_text( size = 14),
           axis.title.y = element_text(size = 14),
           axis.text.x = element_text(size = 14),
           axis.text.y = element_text(size = 14),
           strip.text = element_text( face = "bold", size=14),
           legend.title = element_text(size = 14),
           legend.text = element_text(size = 14),
           legend.position = "bottom") 
dev.off()

########## Simulation Result: Variability ##############

### Figure S6
pdf(file = "~/Desktop/UW Stuff/IS TND/Code/TND TMLE/semiparametric-tmle-tnd/Simulation/Output/MCSD and Mean SE Complete Data.pdf",
    height = 6.6, width = 10)
summ.df %>% filter(samprate == "All") %>%
  ggplot(aes(x = logORSEMean,y = logORMCSD, color = Estimator))+
  theme_bw()+ 
  geom_point(aes(shape = factor(nsamp)), size = 3)+ 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  facet_grid( truelogORmath ~ distributionmath + sampratemath, 
              label = "label_parsed" ) +
  labs(x = "Mean Estimated Standard Error", y = "Monte Carlo Standard Deviation",
       shape = "Sample Size")+ 
  scale_color_manual(values = group.colors) +
  scale_shape_manual(values = c(0, 2, 1, 5)) + 
  guides(color = guide_legend(nrow = 1, order = 1),
         shape = guide_legend(nrow = 1, order = 2))+
  coord_cartesian(
    expand = TRUE,
    default = FALSE,
    clip = "on"
  ) + theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          strip.text = element_text( face = "bold", size = 14), 
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.margin = margin()) 
dev.off()

### Table S2
summ.wide.df <- summ.df %>% 
  pivot_wider(names_from = Estimator, 
              values_from = c(nsims, ModelSamp:SampObsPropMean))

var.tbl0 <- summ.wide.df %>% select(
  trueVE ,truelogORc, distribution, distributionc, nsamp, samprate, 
  logORMCSD_TMLE, logORSEMean_TMLE,
  logORMCSD_PLMx, logORSEMean_PLMx,
  logORSEMean_PLEx, # excluded logORMCSD_PLEx because equals logORMCSD_PLMx
  logORMCSD_MLEx, logORSEMean_MLEx,
  logORMCSD_nPLM, logORSEMean_nPLM, 
   logORSEMean_nPLE, # excluded logORMCSD_nPLE because equals logORMCSD_nPLM
  logORMCSD_nMLE, logORSEMean_nMLE,
  nsims_TMLE, nsims_PLMx) %>% arrange(distribution) 

var.tbl <- var.tbl0 %>% mutate(
  across(where(is.numeric), \(x) format(round(x, digits = 2), nsmall = 2)))

flextable::save_as_docx(flextable(var.tbl),
    path = "~/Desktop/UW Stuff/IS TND/Code/TND TMLE/semiparametric-tmle-tnd/Simulation/Output/Variance Table TMLE.docx")

write.csv(var.tbl, 
  "~/Desktop/UW Stuff/IS TND/Code/TND TMLE/semiparametric-tmle-tnd/Simulation/Output/Variance TMLE.csv")
