######################################################################################################
########## Simulation model for self-management interventions (SMI) in type II diabetes ##############
######################################################################################################

# Original work Copyright (C) 2019 Isaac Corro Ramos & ... 
# Institute for Medical Technology Assessment (iMTA) of the Erasmus University Rotterdam.

# This is the R code of the health economic simulation model for SMI interventions in type II diabetes.
# This model is developed as part of the COMPAR-EU project (https://self-management.eu/)

# Before the simulation code starts, make sure that all the packages below are installed in your computer.
# Then load the packages.

#library(lattice)
#library(MASS)
#library(survival)
#library(plyr)

# When you are reading external files and exporting results you may set a working directory.
# This can be for example the folder where you saved some previous results.
# Do not forget to change this path into your personal working directory.

wd <- setwd("\\\\campus.eur.nl/shared/groups/IMTA-600-H2020COMPAR-EU/600 Projectuitvoer/Diabetes/R code")

# This model makes use of previously estimated regression equations to predict the risk of experiencing diabetes-related 
# complications. These regression equations were sourced from UKPDS 82 and can be found in the paper by Hayes et al. 2013

# UKPDS outcomes model 2: a new version of a model to simulate lifetime health outcomes of patients with type 2 diabetes mellitus
# using data from the 30 year United Kingdom Prospective Diabetes Study: UKPDS 82.
# Hayes AJ, Leal J, Gray AM, Holman RR, Clarke PM. Diabetologia. 2013 Sep;56(9):1925-33. Epub 2013 Jun 22.
# DOI: 10.1007/s00125-013-2940-y

# The simulation model reads these coefficients from previously saved .csv files. The functions below can be used to predict 
# the risk of experiencing diabetes-related complications. --> not from csv, changed my mind

# Notation: as a general rule, names for function input parameters end with "_input".

# Nine complications can occur in the simulation: congestive heart failure (CHF), death, MI, stroke, IHD, amputation, blindness, 
# renal failure and ulcer. For MI, stroke, and amputation, a second ocurrence is possible. 
# The model calculations will rely on the calculations of annual probabilities for those complications.
# The unconditional probability of CHF in the interval t to t+1 is calculated as a function of the cumulative hazard as  
# P = 1 - exp{H(t|x_j) - H(t+1|x_j)}, where H is the cumulative hazard function, t is the duration of diabetes (years) and x_j are 
# the covariates in the equation (also called risk factors). 

# For CHF, IHD, 1st MI for females, 1st and 2nd stroke and 1st amputaiotn with no prior ulcer, 
# H(t|x_j) is assumed to follow a (proportional hazards) Weibull distribution where t is the duration of diabetes.
# The functional form of H(t|x_j) is thus the same, only the regression coefficients and the risk factors will be 
# different for each complication.

# CHECK: if ro = 1 then weibull == exponential and the same function can be used.

cum_haz_weibull <- function(regression_coefficents_input, risk_factors_input, duration_diabetes_input){
  
  risk_factors_input <- as.numeric(risk_factors_input) # delete if not needed
  
  # Note: mind the format of the regression coefficients. It should be a list where the first two coefficients are
  # lambda and ro (with the notation from UKPDS paper - hence the tail(,-2) command) and the rest are the coefficients 
  # for the risk factors associated to each complication.
  linear_predictor <- sum(tail(regression_coefficents_input,n = -2)*c(risk_factors_input))
  
  # Then H returns the value of the cumulative hazard function 
  H <- exp(regression_coefficents_input[1] + linear_predictor)*duration_diabetes_input^regression_coefficents_input[2]
  
  return(list(H = H))
}



### CONTINUE HERE
### NEED TO IMPLEMENT cum_haz_logistic and cum_haz_gompertz



# Regression coefficients from UKPDS paper ESM Table 4: Macrovascular complications 

risk_factors_macrovascular <- c("AFRO", "AGE DIAG", "FEMALE", "INDIAN", "ATFIB", "BMI", "eGFR", "eGFR<60", "HbA1c", "HDL",
                                "LDL", "LDL>35", "MMALB", "PVD", "SBP", "SMOKER", "WBC", "AMP HIST", "CHF HIST", "IHD HIST", 
                                "STROKE HIST", "ULCER HIST")

parameters_macrovascular <- c("lambda", "ro", risk_factors_macrovascular)

CHF <- c(-12.332, 1.514, 0, 0.068,  0, 0, 1.562, 0.072, 0, -0.220, 0, 0,     
         0.012, 0, 0.771, 0.479, 0, 0, 0, 0.658, 0, 0, 0, 0.654)

IHD <- c(-6.709, 1.276, 0, 0.016, -0.532, 0, 0, 0, -0.053, 0, 0, -0.065, 
         0.023, 0, 0, 0.486, 0.058, 0, 0, 0.526, 0.824, 0, 0, 0)

FMIMALE <- c(-8.791, 1, -0.83, 0.045, 0 , 0.279, 0, 0, 0, 0, 0.108, -0.049,
             0.023, 0, 0.203, 0.340, 0.046, 0.277, 0.026, 0.743, 0.814, 0.846, 0.448, 0) # exponential distirbution ro = 1 (2nd element in vector)

FMIFEMALE <- c(-8.708, 1.376, -1.684, 0.041, 0, 0, 0, 0, 0, -0.28, 0.078, 0, 
             0, 0.035, 0.277, 0.469, 0.056, 0.344, 0.07, 0, 0.853, 0.876, 0, 0)

SMI <- c(-4.179, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
         0.021, 0, 0.344, 0, 0, 0, 0, 0, 0, 0, 0, 0) # exponential distirbution ro = 1 (2nd element in vector)

FSTROKE <- c(-13.053, 1.466, 0, 0.066, -0.42, 0, 1.476, 0, 0, -0.190, 0.092, 0, 
             0.016, 0, 0.420, 0, 0.170, 0.331, 0.040, 1.090, 0, 0.481, 0, 0)

SSTROKE <- c(-9.431, 1.956, 0, 0.046, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0.537, 0, 0, 0.656, 0, 0, 0, 0, 0, 0)


macrovascular_risk_equations <- data.frame(CHF, IHD, FMIMALE, FMIFEMALE, SMI, FSTROKE, SSTROKE, row.names = parameters_macrovascular)

macrovascular_risk_equations

# cum_haz_weibull(macrovascular_risk_equations$CHF,c(0,62,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),8)
# cum_haz_weibull(macrovascular_risk_equations$IHD,c(0,62,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),8)
# cum_haz_weibull(macrovascular_risk_equations$FMIMALE,c(0,62,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),8)
# cum_haz_weibull(macrovascular_risk_equations$FMIFEMALE,c(0,62,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),8)


# Regression coefficients from UKPDS paper ESM Table 5: Microvascular complications

risk_factors_microvascular <- c("AFRO", "AGE DIAG", "FEMALE", "ATFIB", "BMI", "eGFR<60", "eGFR>60", "HAEM", "HbA1c", "HDL",
                                "HEART R", "LDL", "MMALB", "PVD", "SBP", "WBC", "AMP HIST", "BLIND HIST", "CHF HIST", 
                                "IHD HIST", "STROKE HIST")

parameters_microvascular <- c("lambda", "ro", risk_factors_microvascular)


BLIND <- c(-11.607, 1, 0, 0.047, 0, 0, 0, 0, 0, 0, 0.171, 0, 0.08, 0, 0, 0, 0.068, 0.052, 0, 0, 0.841, 0.610, 0) #ro =1  --> exponential

ULCER <- c(-11.295, 0, 0, 0.043, 0, 0, 0.053, 0, 0, 0, 0.160, 0, 0, 0, 0, 0.968, 0, 0, 0, 0, 0, 0, 0) # ro is 2nd element. ULCER is logistic so dpont know yet if ro = 1  or 0

FAMPNOULCER <- c(-14.844, 2.067, 0, 0.023, -0.445, 1.088, 0, 0, 0, 0, 0.248, -0.059, 0.098, 0, 0.602, 1.010, 0.086, 0.040, 0, 0, 0, 0, 1.299)

FAMPULCER <- c(-0.811, 1, 0, -0.065, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.769, 0, 0, 0, 0, 0, 0, 0)

SAMP <- c(-3.455, 1, rep(0,8), 0.127, rep(0,12))

RENALF <- c(3.549, 1, 0.686, -0.029, -0.869, 0, -0.054, -1.031, -0.487, -0.268, 0, 0, 0, 0.027, 1.373, 0, 0.085, 0.029, 1.108, 0.732, 0, 0, 0)



microvascular_risk_equations <- data.frame(BLIND, ULCER, FAMPNOULCER, FAMPULCER, SAMP, RENALF, row.names = parameters_microvascular)

microvascular_risk_equations

cum_haz_weibull(microvascular_risk_equations$BLIND,c(0,62,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),8)



# Regression coefficients from UKPDS paper ESM Table 6: Mortality in current year

risk_factors_mortality <- c("FEMALE", "INDIAN", "YEAR", "ATFIB", "BMI1", "BMI3", "CURR AGE", "HDL", "HEART R", 
                            "MMALB", "PVD", "SMOKER", "WBC", "AMP EVENT", "AMP HIST", "AMP2 EVENT", "CHF HIST", 
                            "IHD EVENT", "IHD HIST", "MI EVENT", "MI HIST", "RENAL EVENT", "RENAL HIST", "STROKE EVENT", "STROKE HIST")

parameters_mortality <- c("lambda", "phi", risk_factors_mortality)


# Death in years with no history of events: Gompertz
DEATHNOHIST <- c(-10.908, 0.098, -0.229, rep(0,10), 0.379, rep(0,13))

# Death in first year of events: Logistic
DEATH1YEVENT <- c(-6.916, 0, 0, -0.540, 0.042, rep(0,3), 0.058, 0, 0.124, 0, 0.367, 0.444, 0, 0.321, rep(0,3), 0.423,0, 1.309, 0, 0.584,0, 0.547, 0) #2nd element in logistic not sure if 0 or 1

# Death in years with history but no events: Gompertz
DEATHHISTNOEVENT <- c(-9.207, 0.073, rep(0,4), 1.083, -0.293, rep(0,3), 0.348, 0, 0.374, 0.048, 0, 0.539, 0, 0.632, rep(0, 5), 1.150, 0, 0.473)

# Death in subsequent years of events: Logistic
DEATHYSEVENT <- c(-4.868, 0, rep(0,3), 1.081, 0, 0, 0.050, 0.068, 0, 0, 0.352, 0, 0.089, -1.267, 0.753, -1.727, 0, 0.583, -0.507, 0.982, 0.440, 0, 0.961, -0.619, 0) #2nd element in logistic not sure if 0 or 1


mortality_risk_equations <- data.frame(DEATHNOHIST, DEATH1YEVENT, DEATHHISTNOEVENT, DEATHYSEVENT, row.names = parameters_mortality)

mortality_risk_equations

cum_haz_weibull(mortality_risk_equations$DEATHNOHIST,c(0,62,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),8)
