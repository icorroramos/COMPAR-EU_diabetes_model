######################################################################################################
########## Simulation model for self-management interventions (SMI) in type II diabetes ##############
######################################################################################################

# Original work Copyright (C) 2020 Isaac Corro Ramos & ... 
# Institute for Medical Technology Assessment (iMTA) of the Erasmus University Rotterdam.

# This is the R code used to generate baseline patient characteristics that will be used in the health economic
# simulation model for SMI interventions in type II diabetes developed as part of the COMPAR-EU project (https://self-management.eu/).


# Before the simulation code starts, make sure that all the packages below are installed in your computer.
# Then load the packages.

rm(list = ls(all.names = TRUE))

# library(lattice)
# library(MASS)
# library(survival)
# library(plyr)
#library(tidyverse)

# When you are reading external files and exporting results you may set a working directory.
# This can be for example the folder where you saved some previous results.
# Do not forget to change this path into your personal working directory.

# wd <- setwd("\\\\campus.eur.nl/shared/groups/IMTA-600-H2020COMPAR-EU/600 Projectuitvoer/Diabetes/R code")

# The code below is based on the data generation process described in Appendix B of the paper by Brown et al. 2000 - The global
# diabetes model: user frinedly version 3.0. That appendix describes how to generate data from the KPNW "prevalent" reference population.
# In the same appendix, it is mentioned that since UKPDS data were not available to the authors, the UKPDS reference population was built
# from prediction equations from the KPNW "incident" adjusted to published mean values for the UKPDS baseline population.
# Note then that this is a two-step process.

# QUESTION: in that paper they present equations for the KPNW "prevalent" population. Is it ok to adjust this "prevalent" population to
# match UKPDS or do we need to get data from the "incident" population? Does it matter at all? Or maybe what we have to do is to match
# SMI trial data? Is is important then whether the initial data come from "incident" or "prevalent"?


sample_size <- 10

baseline_characteristics <- c("FEMALE", "MALE", "AFRO", "INDIAN", "ASIAN", "OTHER", "WHITE", "HISPANIC", "CURR.AGE", "YEAR", "AGE.DIAG",
                              "HbA1c", "Tryglicerides", "HDL", "LDL", "SBP", "SMOKER", "MI.HIST", "STROKE.HIST", "CHF.HIST")


# Note I have added MALE just to make the HbA1c claculation below easier

simulation_baseline_patients <- data.frame(matrix(vector(), sample_size, length(baseline_characteristics), 
                                                  dimnames=list(c(), baseline_characteristics)),
                                           stringsAsFactors=F)


### Step 1: dichotomous variables

# Gender
simulation_baseline_patients$FEMALE <- rbinom(sample_size, 1, 1 - 0.511) # 0.511 is the proportion of males in KPNW data
simulation_baseline_patients$MALE <- 1-simulation_baseline_patients$FEMALE

# Ethnicity: in UKPDS equations we have only two ethnicity variables AFRO and INDIAN (the HR referent is caucasian).
# In the KPNW there are more: Black (AFRO), Native American (INDIAN), Asian, Other race, White (Caucasian) 
ETHNIC <- rmultinom(sample_size, size = 1, prob = c(0.027,0.049-0.027,0.075-0.049,0.083-0.075, 1-0.083))

simulation_baseline_patients$AFRO   <- ETHNIC[1,]
simulation_baseline_patients$INDIAN <- ETHNIC[2,]
simulation_baseline_patients$ASIAN  <- ETHNIC[3,]
simulation_baseline_patients$OTHER  <- ETHNIC[4,]
simulation_baseline_patients$WHITE  <- ETHNIC[5,]

# Patients of any race/ethnicity can be Hispanic
simulation_baseline_patients$HISPANIC <- rbinom(sample_size, 1, 0.019) # 0.019 is the proportion of males in KPNW data

# Age and duration of diabetes are sampled a normal distribution.
# In the UKPDS we have AGE.DIAG = age diagnosis, YEAR = diabetes duration and CURR.AGE = current age
# Here we are sampling values of YEAR and CURR.AGE at baseline, i.e. these values are later updated in the simulation 
# AGE.DIAG is defined as CURR.AGE - YEAR

simulation_baseline_patients$CURR.AGE <- round(rnorm(sample_size, 64.671, 12.505),2) # mean and sd in KPNW data
simulation_baseline_patients$YEAR     <- pmax(0, round(rnorm(sample_size, 6.374, 5.126),2)) # mean and sd in KPNW data

# QUESTION: it seems counterintuitive to me that these two variables are not correlated. Moreover, sd in YEAR is relatively large and 
# results in negative values which is impossible.Also, for age the sd seems large... wide range of ages considered. Do no change this 
# for the moment

simulation_baseline_patients$AGE.DIAG <- pmax(0,simulation_baseline_patients$CURR.AGE - simulation_baseline_patients$YEAR)

# QUESTION: although very unlikely, it might be possible to get AGE.DIAG < 0

# Remaining characteristics are obtained from regression equations (OLS) estimated from KPNW data 
# Coefficients below for HbA1c, Tryglicerides, HDL, LDL and SBP were obtained from Table 4 in Appendix B of the paper by Brown et al. 2000  

# HbA1c: predictors are CURR.AGE, MALE, YEAR, AFRO, INDIAN and HISPANIC (plus the regression intercept and the RMSE)
HbA1c_predictors <- c("CURR.AGE", "MALE", "YEAR", "AFRO", "INDIAN", "HISPANIC")
simulation_baseline_patients$HbA1c <- (8.817 + as.matrix(simulation_baseline_patients[HbA1c_predictors]) %*% c(-0.020, 0.266, 0.032, -0.016, -0.130, -0.061)) + rnorm(sample_size, 0, 1.672) 

# Tryglicerides: predictors are the same as the previous ones plus HbA1c (plus the regression intercept and the RMSE)

# QUESTION: I think this is not in UKPDS but it is needed here nevertheless to predict HDL, LDL and SBP

Tryglicerides_predictors <- c(HbA1c_predictors, "HbA1c")
simulation_baseline_patients$Tryglicerides <- (309.622 + as.matrix(simulation_baseline_patients[Tryglicerides_predictors]) %*% c(-1.899, -7.863, -1.513, -71.97, -6.497, -23.543, 8.212)) + rnorm(sample_size, 0, 257.140) 

# QUESTION: Tryglicerides results in negative values

# HDL: predictors are the same as the previous ones plus Tryglicerides (plus the regression intercept and the RMSE)
HDL_predictors <- c(Tryglicerides_predictors, "Tryglicerides")
simulation_baseline_patients$HDL <- (43.688 + as.matrix(simulation_baseline_patients[HDL_predictors]) %*% c(0.048, -4.078, 0.084, 1.992, -1.079, -0.349, 0.109, -0.006)) + rnorm(sample_size, 0, 12.555) 


# LDL: predictors are the same as the previous ones plus HDL (plus the regression intercept and the RMSE)
LDL_predictors <- c(HDL_predictors, "HDL")
simulation_baseline_patients$LDL <- (121.497 + as.matrix(simulation_baseline_patients[LDL_predictors]) %*% c(-0.040, -1.678, -0.038, -1.371, 5.005, 2.368, 0.268, 0.000, 0.059)) + rnorm(sample_size, 0, 35.103) 

# SBP: predictors are the same as the previous ones plus LDL (plus the regression intercept and the RMSE)
SBP_predictors <- c(LDL_predictors, "LDL")
simulation_baseline_patients$SBP <- (125.727 + as.matrix(simulation_baseline_patients[SBP_predictors]) %*% c(0.209, -3.277, 0.135, 2.555, -0.851, -4.746, -0.095, 0.002, 0.007, 0.000)) + rnorm(sample_size, 0, 19.685) 

# Note in the global diabetes model, it is possible to enter target mean values for each varible predicted by OLS.

# QUESTION: come back to this later but we wont do it right now.

# Finally, there are also dichotomous variables generated from logistic regression equations. 
# The set of predictors for these variables is fixed and includes all the previous ones (including SBP)
all_predictors <- c(SBP_predictors, "SBP")

# The remaining predicted variables are then as follows:

# SMOKER
pred_SMOKER <- 1.385 + as.matrix(simulation_baseline_patients[all_predictors]) %*% c(-0.053, 0.134, 0.007, -0.379, 0.567, -0.221, -0.073, 0.000, -0.004, 0.002, -0.001)
prob_SMOKER <- exp(pred_SMOKER)/(1+exp(pred_SMOKER))
simulation_baseline_patients$SMOKER <- sapply(prob_SMOKER, function(x) rbinom(1,1,x))

# QUESTION: there are events in the KPNW data that are not included in the UKPS equations. These "other" events are not used 
# to predict anything else either. For the moment, I will not include them in the model.
# Events in KPNW data but not in UKPDS are CVD [IS THIS THE SAME AS PVD???], CHD, angina, intermittent claudication, PAD, 
# neuropathy, nephropathy, BDR (what's this?), PDR (what's this?), Macular edema and retinal laser treatment.

# QUESTION: for events like MI or stroke, which are included in UKPDS, it is possible to predict them at baseline. So these will be
# for example the variables MI.HIST. We need to make sur ein the simulation whether it is allowed to have history of an event at
# baseline or not (i.e. history can only be applied after a simulated event). We also need to decide what would make more sense.

# MI.HIST
pred_MI <- 0.575 + as.matrix(simulation_baseline_patients[all_predictors]) %*% c(-0.003, 0.169, 0.006, 0.08, -0.520, 0.073, -0.002, -0.00, -0.001, -0.002, -0.005)
prob_MI <- exp(pred_MI)/(1+exp(pred_MI))
simulation_baseline_patients$MI.HIST <- sapply(prob_MI, function(x) rbinom(1,1,x))

# STROKE.HIST
pred_STROKE <- -3.263 + as.matrix(simulation_baseline_patients[all_predictors]) %*% c(0.022, -0.150, 0.01, 0.049, -0.401, -1.165, -0.027, 0.00, 0.002, 0.00, 0.007)
prob_STROKE <- exp(pred_STROKE)/(1+exp(pred_STROKE))
simulation_baseline_patients$STROKE.HIST <- sapply(prob_STROKE, function(x) rbinom(1,1,x))

# CHF.HIST
pred_CHF <- 0.04 + as.matrix(simulation_baseline_patients[all_predictors]) %*% c(0.019, -0.124, 0.024, -0.250, 0.251, 0.317, 0.012, 0.00, -0.005, -0.001, -0.008)
prob_CHF <- exp(pred_CHF)/(1+exp(pred_CHF))
simulation_baseline_patients$CHF.HIST <- sapply(prob_CHF, function(x) rbinom(1,1,x))

# Continue with paper Table 6 in Appendix B with Amputation (LEA), MMALB (microalbuminuria), ESRD (IS THIS RENAL FAILURE?), Severe vision loss (assume to be blindness in UKPDS)

# When finish with above: Check -->> Are all characteristics needed for UKPDS included here???
