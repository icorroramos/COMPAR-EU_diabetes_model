######################################################################################################
########## Simulation model for self-management interventions (SMI) in type II diabetes ##############
######################################################################################################

# Original work Copyright (C) 2020 Isaac Corro Ramos & ... 
# Institute for Medical Technology Assessment (iMTA) of the Erasmus University Rotterdam.

# This is the R code of the health economic simulation model for SMI interventions in type II diabetes.
# This model is developed as part of the COMPAR-EU project (https://self-management.eu/)

# Before the simulation code starts, make sure that all the packages below are installed in your computer.
# Then load the packages.

rm(list = ls(all.names = TRUE))

# library(lattice)
# library(MASS)
# library(survival)
# library(plyr)
# library(tidyverse)

# When you are reading external files and exporting results you may set a working directory.
# This can be for example the folder where you saved some previous results.
# Do not forget to change this path into your personal working directory.

wd <- setwd("\\\\campus.eur.nl/shared/groups/IMTA-600-H2020COMPAR-EU/600 Projectuitvoer/Diabetes/R code")


#######################################################
########## PART I - UKPDS RISK EQUATIONS ##############
#######################################################

# This model makes use of previously estimated regression equations to predict the risk of experiencing diabetes-related 
# complications. These regression equations were sourced from UKPDS 82 and can be found in the paper by Hayes et al. 2013:

# UKPDS outcomes model 2: a new version of a model to simulate lifetime health outcomes of patients with type 2 diabetes mellitus
# using data from the 30 year United Kingdom Prospective Diabetes Study: UKPDS 82.
# Hayes AJ, Leal J, Gray AM, Holman RR, Clarke PM. Diabetologia. 2013 Sep;56(9):1925-33. Epub 2013 Jun 22.
# DOI: 10.1007/s00125-013-2940-y

# The simulation model reads these coefficients and uses the functions below to predict the annual risk of experiencing the following 
# nine diabetes-related complications: 

# 1. Congestive heart failure (CHF) [Assumption: It can happen just once as event. After that, it is considered as history.]
# 2. Myocardial Infarction (MI) [Assumption: There is a distinction between first and second. Third and above is not possible.]
# 3. Stroke [Assumption: There is a distinction between first and second. Third and above is not possible.]
# 4. Ischaemic Heart Disease (IHD) [Assumption: It can happen just once as event. After that, it is considered as history.]
# 5. Amputation [Assumption: There is a distinction between first and second. Third and above is not possible.]
# 6. Blindness [Assumption: It can happen just once as event, even though it is defined as blindness in one eye (UKPDS OM1 paper). After that, it is considered as history.]
# 7. Renal failure [Assumption: It can happen just once as event. After that, it is considered as history.]
# 8. Diabetic Ulcer [Assumption: It can happen just once as event. After that, it is considered as history.]
# 9. Death [It can happen just once as event.]

# Multiple events are possible in the same year but only the number of times specified above. 


# The risk factors and the regression coefficients for each complication were taken from the UKPDS paper and are defined below.

# Macrovascular complications: UKPDS paper ESM Table 4. Macrovascular complications include CHF, IHD, MI and stroke.
# The risk factors used to predict macrovascular complications are the following (notation used in the model, definition from 
# UKPDS paper ESM Table 2):

# 1. AFRO: 1 == afro-caribbean ethnicity, 0 == Otherwise. HR referent == Caucasian
# 2. AGE.DIAG: age in years at diagnosis of diabetes. HR per year increase in age at diagnosis
# 3. FEMALE: 1 == female, 0 == male. HR referent == male
# 4. INDIAN: 1 == asian-indian ethnicity, 0 == Otherwise. HR referent == Caucasian
# 5. ATFIB: 1 == attrial fibrilation, 0 == otherwise. Defined from Minnesota codes 831 and 833. HR referent == no ATFIB

# QUESTION: look for definiiton of Minessota codes 831 and 833

# 6. BMI: Body mass index (m/kg^2) measured continuosly. HR per unit increase in BMI
# 7. eGFR: Estimated glomerular filtration rate (ml/min/1.73m^2) from modification of diet in renal disease (MDRD) formula.
#          Continuosly measured and further DIVIDED by 10. HR per 10 ml/min/1.73m^2 increase 
# 8. eGFR60less: Same as eGFR. Continuous spline (knot at 60) and further DIVIDED by 10. HR per 10 ml/min/1.73m^2 increase if < 60

# QUESTION: look for exact definition of the spline funciton used. There is also a variable eGFR60more. 

# 9. HbA1c: percentage HbA1c measured continuously. HR per 1% increase in HbA1c
# 10. HDL: high density lipoprotein cholesterol (mmol/l) measured continuosly and MULTIPLIED by 10. HR per 0.1 mg/dL increase in HDL
# 11. LDL: low density lipoprotein cholesterol (mmol/l) measured continuosly and MULTIPLIED by 10. HR per 0.1 mg/dL increase in LDL
# 12. LDL35more: same as LDL. Continuos spline (knot at 3.5) and MULTIPLIED by 10. HR per 0.1 mg/dL increase in LDL if > 3.5 mg/dL
# 13. MMALB: presence of micro- or macro-albuminuria. 1 == unrine albumin >= 50 mg/l, 0 == otherwise. HR referent == no albuminuria
# 14. PVD: 1 == peripheral vascular disease, 0 == otherwise. Defined from presence of intermittent claudication or ankle brachial
#          pressure index < 0.9. HR referent == no presenceevidence of PVD
# 15. SBP: systolic blood pressure (mm Hg) measured continuously and further DIVIDED by 10. HR per 10 mm Hg increase in SBP
# 16. SMOKER: 1 == current smoker, 0 == otherwise. HR referent == not current smoker
# 17. WBC: white blood cell count measured continuously. HR per 1x10^6/ml increase
# 18. AMP.HIST: 1 == history of amputation, 0 == otherwise. HR referent == no prior amputation
# 19. CHF.HIST: 1 == history of CHF, 0 == otherwise. HR referent == no prior CHF
# 20. IHD.HIST: 1 == history of IHD, 0 == otherwise. HR referent == no prior IHD
# 21. STROKE.HIST: 1 == history of stroke, 0 == otherwise. HR referent == no prior stroke
# 22. ULCER.HIST: 1 == history of diabetic ulcer, 0 == otherwise. HR referent == no prior diabetic ulcer

# Note: pay extra attention to the data transformation. Otherwise, the equations will not make sense!

# The vector below contains the names of all risk factors used to predict macrovacular complications in the model.
risk_factors_macrovascular <- c("AFRO", "AGE.DIAG", "FEMALE", "INDIAN", "ATFIB", "BMI", "eGFR", "eGFR60less", "HbA1c", "HDL",
                                "LDL", "LDL35more", "MMALB", "PVD", "SBP", "SMOKER", "WBC", "AMP.HIST", "CHF.HIST", "IHD.HIST", 
                                "STROKE.HIST", "ULCER.HIST")

# Besides the risk factors, there are two other parameters used in the UKPDS equations called "lambda" and "ro". We added them below.
# Now we have the names of all the parameters (lambda, ro and risk factors) used to predict macrovascular risks in the model.
parameters_macrovascular <- c("lambda", "ro", risk_factors_macrovascular)

# And now, below we have the regression coefficients as reported in the UKPDS paper ESM Table 4.
# Note that not all risk factors are used to predict all complications. For that reason, you see some 0's in the regression coefficients below.

# The format for each vector below is c(lambda, ro, risk factors) where the values are taken from the UKPDS paper ESM Table 4.
CHF <- c(-12.332, 1.514, 0, 0.068,  0, 0, 1.562, 0.072, 0, -0.220, 0, 0, 0.012, 0, 0.771, 0.479, 0, 0, 0, 0.658, 0, 0, 0, 0.654)
IHD <- c(-6.709, 1.276, 0, 0.016, -0.532, 0, 0, 0, -0.053, 0, 0, -0.065, 0.023, 0, 0, 0.486, 0.058, 0, 0, 0.526, 0.824, 0, 0, 0)
# first MI for males: exponential distirbution ro = 1 (2nd element in vector) --> explanatoin below in "annual_p_weibull" function
FMIMALE <- c(-8.791, 1, -0.83, 0.045, 0 , 0.279, 0, 0, 0, 0, 0.108, -0.049, 0.023, 0, 0.203, 0.340, 0.046, 0.277, 0.026, 0.743, 0.814, 0.846, 0.448, 0) 
# first MI for females
FMIFEMALE <- c(-8.708, 1.376, -1.684, 0.041, 0, 0, 0, 0, 0, -0.28, 0.078, 0, 0, 0.035, 0.277, 0.469, 0.056, 0.344, 0.07, 0, 0.853, 0.876, 0, 0)
# second MI: exponential distirbution ro = 1 (2nd element in vector)
SMI <- c(-4.179, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.021, 0, 0.344, 0, 0, 0, 0, 0, 0, 0, 0, 0) 
# first stroke
FSTROKE <- c(-13.053, 1.466, 0, 0.066, -0.42, 0, 1.476, 0, 0, -0.190, 0.092, 0, 0.016, 0, 0.420, 0, 0.170, 0.331, 0.040, 1.090, 0, 0.481, 0, 0)
# second stroke
SSTROKE <- c(-9.431, 1.956, 0, 0.046, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.537, 0, 0, 0.656, 0, 0, 0, 0, 0, 0)

# QUESTION: check that all values have been copied properly!

# Below we simply create a table (R data frame) with all the coefficients of the regresison equations used to predict macrovascular complications.
# When we define the risk equations below in the code, these will read the coefficients from this table in order to predict the annual 
# probability of experiencing a mcrovascular complication.
macrovascular_risk_equations <- data.frame(CHF, IHD, FMIMALE, FMIFEMALE, SMI, FSTROKE, SSTROKE, row.names = parameters_macrovascular)
# macrovascular_risk_equations


# Microvascular complications: UKPDS paper ESM Table 5. Microvascular complications include blindness, diabetic ulcer amputation and renal failure.
# The risk factors (the ones that are not defined above for macrovascular complications) used to predict microvascular complications 
# are the following (notation used in the model, definition from UKPDS paper ESM Table 2):

# 23. eGFR60more: Same as eGFR. Continuous spline (knot at 60) and further DIVIDED by 10. HR per 10 ml/min/1.73m^2 increase if > 60
# 24. HAEM: haemoglobin g/dL measured continuously. HR per 1 g/dL increase
# 25. HEART.R: heart rate (beats per minute) determined from inspiration/expiration RR on ECG. Continuously measued and further DIVIDED by 10.
#              HR per 10 bpm increase.
# 26. BLIND.HIST: 1 == history of blindness, 0 == otherwise. HR referent == no prior blindness

# Note: pay extra attention to the data transformation. Otherwise, the equations will not make sense!

# The vector below contains the names of all risk factors used to predict microvascular complications in the model.
risk_factors_microvascular <- c("AFRO", "AGE.DIAG", "FEMALE", "ATFIB", "BMI", "eGFR60less", "eGFR60more", "HAEM", "HbA1c", "HDL",
                                "HEART.R", "LDL", "MMALB", "PVD", "SBP", "WBC", "AMP.HIST", "BLIND.HIST", "CHF.HIST", 
                                "IHD.HIST", "STROKE.HIST")

# Besides the risk factors, there are two other parameters used in the UKPDS equations called "lambda" and "ro". We added them below.
# Now we have the names of all the parameters (lambda, ro and risk factors) used to predict microvascular risks in the model.
parameters_microvascular <- c("lambda", "ro", risk_factors_microvascular)

# And now, below we have the regression coefficients as reported in the UKPDS paper ESM Table 5.
# Note that not all risk factors are used to predict all complications. For that reason, you see some 0's in the regression coefficients below.
# The format for each vector below is c(lambda, ro, risk factors) where the values are taken from the UKPDS paper ESM Table 5.

# Blindness: ro == 1  --> exponential distribution
BLIND <- c(-11.607, 1, 0, 0.047, 0, 0, 0, 0, 0, 0, 0.171, 0, 0.08, 0, 0, 0, 0.068, 0.052, 0, 0, 0.841, 0.610, 0) 
# Diabetic ulcer: 
ULCER <- c(-11.295, 0, 0, 0.043, 0, 0, 0.053, 0, 0, 0, 0.160, 0, 0, 0, 0, 0.968, 0, 0, 0, 0, 0, 0, 0) 
# QUESTiON: ro is 2nd element. ULCER is logistic so dont know yet if ro = 1  or 0

# First amputation, no prior ulcer history
FAMPNOULCER <- c(-14.844, 2.067, 0, 0.023, -0.445, 1.088, 0, 0, 0, 0, 0.248, -0.059, 0.098, 0, 0.602, 1.010, 0.086, 0.040, 0, 0, 0, 0, 1.299)
# First amputation, prior ulcer history
FAMPULCER <- c(-0.811, 1, 0, -0.065, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.769, 0, 0, 0, 0, 0, 0, 0)
# Second amputation
SAMP <- c(-3.455, 1, rep(0,8), 0.127, rep(0,12))
# Renal failure
RENALF <- c(3.549, 1, 0.686, -0.029, -0.869, 0, -0.054, -1.031, -0.487, -0.268, 0, 0, 0, 0.027, 1.373, 0, 0.085, 0.029, 1.108, 0.732, 0, 0, 0)

# QUESTION: check that all values have been copied properly!

# Below we simply create a table (R data frame) with all the coefficients of the regresison equations used to predict microvascular complications.
# When we define the risk equations below in the code, these will read the coefficients from this table in order to predict the annual 
# probability of experiencing a microvascular complication.
microvascular_risk_equations <- data.frame(BLIND, ULCER, FAMPNOULCER, FAMPULCER, SAMP, RENALF, row.names = parameters_microvascular)
#microvascular_risk_equations


# Risk of death: UKPDS paper ESM Table 6. Four equations depending on events and history. The risk factors (the ones that are not 
# defined above for macro/microvascular complications) are the following (notation used in the model, definition from UKPDS paper ESM Table 2):

# 27. YEAR: duration of diabetes in years measured continuously. HR per year increase in duration of diabetes. 
#           Note it increases + 1 after each simulated year.
# 28. BMI1: 1 == BMI < 18.5 m/kg^2, 0 == otherwise. HR referent 18.5 m/kg^2 <= BMI < 25 m/kg^2. It depends thus on the variable BMI.
# 29. BMI3: 1 == BMI >= 25 m/kg^2, 0 == otherwise. HR referent 18.5 m/kg^2 <= BMI < 25 m/kg^2. It depends thus on the variable BMI.
# 30. CURR.AGE: current age in years measured continuously. HR per year increase in current age. Note CURR.AGE = AGE.DIAG + YEAR at the 
#               beginning of the simulation. Also after each year in the simulation if YEAR is propoerly updated.
# 31. AMP.EVENT: 1 == first amputation, 0 == otherwise. HR referent == no first amputation event

# QUESTION: consider renaming this variable to AMP1.EVENT to make it clearer!

# 32. AMP2.EVENT: 1 == second amputation, 0 == otherwise. HR referent == no second amputation event
# 33. IHD.EVENT: 1 == IHD, 0 == otherwise. HR referent == no IHD event
# 34. MI.EVENT: 1 == MI, 0 == otherwise. HR referent == no MI event

# QUESTION: given that it seems that only 2 MI events are possible, consider creating two variables to make it clearer! These variables would
# be used to control that no more than 2 MIs are possible but they should be treated as one variable to predict death.

# 35. MI.HIST: 1 == history of MI, 0 == otherwise. HR referent == no prior MI
# 36. RENAL.EVENT: 1 == renal failure event, 0 == otherwise. HR referent == no renal failure event
# 37. RENAL.HIST: 1 == history of renal failure, 0 == otherwise. HR referent == no prior renal failure
# 38. STROKE.EVENT: 1 == stroke, 0 == otherwise. HR referent == no stroke event

# QUESTION: given that it seems that only 2 stroke events are possible, consider creating two variables to make it clearer! These variables would
# be used to control that no more than 2 strokes are possible but they should be treated as one variable to predict death.



# Note: pay extra attention to the data transformation. Otherwise, the equations will not make sense!

# The vector below contains the names of all risk factors used to predict death in the model. 
risk_factors_mortality <- c("FEMALE", "INDIAN", "YEAR", "ATFIB", "BMI1", "BMI3", "CURR.AGE", "HDL", "HEART.R", 
                            "MMALB", "PVD", "SMOKER", "WBC", "AMP.EVENT", "AMP.HIST", "AMP2.EVENT", "CHF.HIST", 
                            "IHD.EVENT", "IHD.HIST", "MI.EVENT", "MI.HIST", "RENAL.EVENT", "RENAL.HIST", "STROKE.EVENT", "STROKE.HIST")

# Besides the risk factors, there are two other parameters used in the UKPDS equations called "lambda" and "phi". We added them below.
# Now we have the names of all the parameters (lambda, phi and risk factors) used to predict death in the model.
parameters_mortality <- c("lambda", "phi", risk_factors_mortality)


# And now, below we have the regression coefficients as reported in the UKPDS paper ESM Table 6.
# Note that not all risk factors are used to predict all types of death. For that reason, you see some 0's in the regression coefficients below.
# The format for each vector below is c(lambda, phi, risk factors) where the values are taken from the UKPDS paper ESM Table 6.

# Death in years with no history or events: Gompertz distribution
DEATHNOHIST <- c(-10.908, 0.098, -0.229, rep(0,10), 0.379, rep(0,13))

# Death in first year of events (excluding blindness or ulcer): Logistic distribution
DEATH1YEVENT <- c(-6.916, 0, 0, -0.540, 0.042, rep(0,3), 0.058, 0, 0.124, 0, 0.367, 0.444, 0, 0.321, rep(0,3), 0.423,0, 1.309, 0, 0.584,0, 0.547, 0) 

# QUESTION: 2nd element in logistic not sure if 0 or 1

# Death in years with history but no events (when there is history of ANY event): Gompertz distribution
DEATHHISTNOEVENT <- c(-9.207, 0.073, rep(0,4), 1.083, -0.293, rep(0,3), 0.348, 0, 0.374, 0.048, 0, 0.539, 0, 0.632, rep(0, 5), 1.150, 0, 0.473)

# Death in subsequent years of events(excluding blindness or ulcer): Logistic distribution
DEATHYSEVENT <- c(-4.868, 0, rep(0,3), 1.081, 0, 0, 0.050, 0.068, 0, 0, 0.352, 0, 0.089, -1.267, 0.753, -1.727, 0, 0.583, -0.507, 0.982, 0.440, 0, 0.961, -0.619, 0) 

# QUESTION: 2nd element in logistic not sure if 0 or 1 --> Think it is 0 if we look at paper! Not sure why I wrote that then...

# QUESTION: check that all values have been copied properly!



# Below we simply create a table (R data frame) with all the coefficients of the regresison equations used to predict death.
# When we define the risk equations below in the code, these will read the coefficients from this table in order to predict the annual 
# probability of death.
mortality_risk_equations <- data.frame(DEATHNOHIST, DEATH1YEVENT, DEATHHISTNOEVENT, DEATHYSEVENT, row.names = parameters_mortality)
#mortality_risk_equations



# Below we consider all the patient characteristics (38) that "define" a patient in this model. 
risk_factors_simulation <- unique(sort(c(risk_factors_macrovascular, risk_factors_microvascular, risk_factors_mortality)))

# risk_factors_simulation
# [1] "AFRO"         "AGE.DIAG"     "AMP.EVENT"    "AMP.HIST"     "AMP2.EVENT"   "ATFIB"        "BLIND.HIST"   "BMI"         
# [9] "BMI1"         "BMI3"         "CHF.HIST"     "CURR.AGE"     "eGFR"         "eGFR60less"   "eGFR60more"   "FEMALE"      
# [17] "HAEM"         "HbA1c"        "HDL"          "HEART.R"      "IHD.EVENT"    "IHD.HIST"     "INDIAN"       "LDL"         
# [25] "LDL35more"    "MI.EVENT"     "MI.HIST"      "MMALB"        "PVD"          "RENAL.EVENT"  "RENAL.HIST"   "SBP"         
# [33] "SMOKER"       "STROKE.EVENT" "STROKE.HIST"  "ULCER.HIST"   "WBC"          "YEAR"        


# QUESTION: Not sure if we want to add a patient id (PTID) to allow sampling the same patient multiple times.

# The simulation model relies on the calculations of annual probabilities for the above mentioned complications.
# For example, the unconditional probability of CHF in the interval t to t+1 is calculated as a function of the cumulative hazard as  
# P = 1 - exp{H(t|x_j) - H(t+1|x_j)}, where H is the cumulative hazard function, t is the duration of diabetes (years) and x_j are 
# the covariates in the equation (also called risk factors). 

# For CHF, IHD, 1st MI for females, 1st and 2nd stroke and 1st amputation with no prior ulcer, 
# H(t|x_j) is assumed to follow a (proportional hazards) Weibull distribution where t is the duration of diabetes.
# The functional form of H(t|x_j) is thus the same, only the regression coefficients and the risk factors will be 
# different for each complication. The function used to calculate the annual probabilities of experiencing those events
# is called "annual_p_weibull" and it is defined below. Note that, with the notation used in the UKPDS paper, if ro = 1 then 
# weibull distribution == exponential distribution. Thus, the same function "annual_p_weibull" can be used when an exponential
# distribution is assumed.

# Notation: as a general rule, names for function input parameters end with "_input".

annual_p_weibull <- function(regression_coefficents_input, risk_factors_input, duration_diabetes_input){
  
  risk_factors_input <- as.numeric(risk_factors_input) # delete if not needed
  
  # Note: mind the format of the regression coefficients. It should be a list where the first two coefficients are
  # lambda and ro (with the notation from UKPDS paper - hence the tail(,-2) command) and the rest are the coefficients 
  # for the risk factors associated to each complication.
  linear_predictor <- sum(tail(regression_coefficents_input,n = -2)*c(risk_factors_input))
  
  # Then H returns the value of the cumulative hazard function 
  H_t1 <- exp(regression_coefficents_input[1] + linear_predictor)*duration_diabetes_input^regression_coefficents_input[2]
  H_t2 <- exp(regression_coefficents_input[1] + linear_predictor)*(1+duration_diabetes_input)^regression_coefficents_input[2]
  p = 1 - exp(H_t1 -H_t2)
  
  return(list(H_t1 = H_t1, # no need to return the H's, keep them for now just for validation purposes
              H_t2 = H_t2,
              p = p))
}


# For death in years with no history or events and in years with history but no events, H(t|x_j) is assumed to follow a 
# Gompertz distribution where t is the current age. The functional form of H(t|x_j) is the same for both types of death. 
# Only the regression coefficients and the risk factors are different for each type of death. 
# The function used to calculate the annual death probabilities is called "annual_p_gompertz" and it is defined below. 
annual_p_gompertz <- function(regression_coefficents_input, risk_factors_input, current_age_input){
  
  risk_factors_input <- as.numeric(risk_factors_input) # delete if not needed
  
  # Note: mind the format of the regression coefficients. It should be a list where the first two coefficients are
  # lambda and ro (with the notation from UKPDS paper - hence the tail(,-2) command) and the rest are the coefficients 
  # for the risk factors associated to each complication.
  linear_predictor <- sum(tail(regression_coefficents_input,n = -2)*c(risk_factors_input))
  
  # Then H returns the value of the cumulative hazard function 
  H_t1 <- ((regression_coefficents_input[2])^-1)*exp(regression_coefficents_input[1] + linear_predictor)*(exp(current_age_input*regression_coefficents_input[2])-1)
  H_t2 <- ((regression_coefficents_input[2])^-1)*exp(regression_coefficents_input[1] + linear_predictor)*(exp((1+current_age_input)*regression_coefficents_input[2])-1)
  p = 1 - exp(H_t1 - H_t2)
  
  return(list(H_t1 = H_t1, # no need to return the H's, keep them for now just for validation purposes
              H_t2 = H_t2,
              p = p))
}

# The annual probability of experiencing diabetic ulcer, death in first year of event and death in subsequent years of events, 
# is assumed to follow a logistic distribution where t is the current age. The functional form of the probability function is the same for 
# diabetic ulcer and both types of death. Only the regression coefficients and the risk factors are different. 
# The function used to calculate the annual death or ulcer probabilities is called "annual_p_logistic" and it is defined below. 
annual_p_logistic <- function(regression_coefficents_input, risk_factors_input){
  
  risk_factors_input <- as.numeric(risk_factors_input) # delete if not needed
  
  # Note: mind the format of the regression coefficients. It should be a list where the first two coefficients are
  # lambda and ro (with the notation from UKPDS paper - hence the tail(,-2) command) and the rest are the coefficients 
  # for the risk factors associated to each complication. Note for logistic, ro (also called phi in ESM Table 6) == 0.
  
  linear_predictor <- regression_coefficents_input[1] + sum(tail(regression_coefficents_input,n = -2)*c(risk_factors_input))
  
  # Then p returns the annual probability 
  p = 1 - (exp(-linear_predictor)/(1+exp(-linear_predictor)))
  
  return(list(p = p))
}


#####

# Small validation exercise. Delete when ok.

# Reproduce examples in statistical appendix of UKPDS paper:
# annual_p_weibull is a function of (regression_coefficents_input, risk_factors_input, duration_diabetes_input)
# risk_factors_input for CHF are risk_factors_macrovascular so c("AFRO", "AGE.DIAG", "FEMALE", "INDIAN", "ATFIB", "BMI", "eGFR", "eGFR60less", "HbA1c", "HDL",
# "LDL", "LDL35more", "MMALB", "PVD", "SBP", "SMOKER", "WBC", "AMP.HIST", "CHF.HIST", "IHD.HIST", 
# "STROKE.HIST", "ULCER.HIST")

# Probability of heart failure in the current year for a male, 70 years old, with 8 years diabetes, LDL 3.0 mmol/l, BMI = 32, eGFR 50, with microalbuminuria and history of amputation:

# annual_p_weibull(macrovascular_risk_equations$CHF,c(0,62,0,0,0,32,0,50/10,0,0,3*10,0,1,0,0,0,0,1,0,0,0,0),8)
# $`H_t1`
# [1] 0.1387994
# 
# $H_t2
# [1] 0.1658947
# 
# $p
# [1] 0.02673152

# Same as in the paper.



# Reproduce examples in statistical appendix of UKPDS paper:
# annual_p_logistic is a function of (regression_coefficents_input, risk_factors_input) 
# risk factors inputs for mortality are ("FEMALE", "INDIAN", "YEAR", "ATFIB", "BMI1", "BMI3", "CURR.AGE", "HDL", "HEART.R", 
# "MMALB", "PVD", "SMOKER", "WBC", "AMP.EVENT", "AMP.HIST", "AMP2.EVENT", "CHF.HIST", 
# "IHD.EVENT", "IHD.HIST", "MI.EVENT", "MI.HIST", "RENAL.EVENT", "RENAL.HIST", "STROKE.EVENT", "STROKE.HIST") 
# Probability of death in the current year for a male, 70 years old, smoker, with 12 years of diabetes, heart rate 80 bpm, has an MI but no history of other events: 
# annual_p_logistic(mortality_risk_equations$DEATH1YEVENT,c(0,0,12,0,0,0,70,0,80/10,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0))
# $`p`
# [1] 0.5970047
# Same as in the paper.

# Same person with heart failure (note this probability does not depend on heart failure as predicting factor)
# annual_p_logistic(mortality_risk_equations$DEATH1YEVENT,c(0,0,12,0,0,0,70,0,80/10,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0))
# $`p`
# [1] 0.2857736
# Same as in the paper.


# Reproduce examples in statistical appendix of UKPDS paper:
# annual_p_gompertz is a function of (regression_coefficents_input, risk_factors_input, current_age_input)
# risk factors inputs for mortality are ("FEMALE", "INDIAN", "YEAR", "ATFIB", "BMI1", "BMI3", "CURR.AGE", "HDL", "HEART.R", 
# "MMALB", "PVD", "SMOKER", "WBC", "AMP.EVENT", "AMP.HIST", "AMP2.EVENT", "CHF.HIST", 
# "IHD.EVENT", "IHD.HIST", "MI.EVENT", "MI.HIST", "RENAL.EVENT", "RENAL.HIST", "STROKE.EVENT", "STROKE.HIST") 
# Probability of death in the current year for a male, 70 years old, smoker, with 12 years of diabetes, overweight, white blood cell count 6*10^6/ml, history of heart failure but no events in current year:
# annual_p_gompertz(mortality_risk_equations$DEATHHISTNOEVENT,c(0,0,12,0,0,1,70,0,0,0,0,1,6,0,0,0,1,0,0,0,0,0,0,0,0),70)
# $`H_t1`
# [1] 0.6158459
# 
# $H_t2
# [1] 0.6627674
# 
# $p
# [1] 0.04583776
# I do not get the same numbers as in the paper. But also with the calculations in the paper, I do not get the same numbers.

# QUESTION: double-check the "annual_p_gompertz" function

#####


##############################################
########## PART II - SIMULATION ##############
##############################################


# --> Work in progress...

# The simulation ocnsists of one main function called "SMDMII_model_simulation". 
# This function is used to 1) simulate patients' clinical history, 2) calculate costs and 3) calculate QALYs. 

# QUESTION: Focus for now in (deterministic) clinical history only. Costs, QALYs and PSA to be calculated later.

# The input parameters of the "SMDMII_model_simulation" function are the following:

# 1. patient_size_input: number of patients included in the simulation.
# 2. run_PSA_input: 1 == runs the model in probabilistic mode, 0 == deterministic.
# 3. treatment_effect_input: variable to indicate treatment effect. Default should be 1.

# QUESTION: Treatment effect parameter(s) to be defined. For the moment just a dummy.
#           Note there might be probably one input for each effect modifier: e.g. on events CHF, MI, stroke, etc. or risk factors.

# 4. seed_input: A random seed that it is used to ensure consistency in the model results. 

# QUESTION: Not sure if needed but keep it for the moment. 

SMDMII_model_simulation <- function(patient_size_input,
                                    run_PSA_input,
                                    treatment_effect_input, 
                                    seed_input){
  
  # Step 1: Read the data with the patient characteristics.
  
  # QUESTION: not sure how we will do this yet. In the COPD model we read them from external excel files.
  
  
  all_baseline_patients <- matrix(c(rep(0,length(risk_factors_simulation)),1,0),ncol = 2 + length(risk_factors_simulation)) 
  
  # 1 = value for SDURATION AT THE BEGINNING
  
  # Step 2: indicate the patient characteristics that we will save during the simulation. 
  #         Besides the risk factors, we have a simulation ID, may have a patient ID,
  #         duration of diabetes (not necesarily 0 at baseline?) and indicator variable for dead
  
  #history_characteristics <- c("SIMID","PTID",risk_factors_simulation, "SDURATION","dead")
  history_characteristics <- c("SIMID", risk_factors_simulation, "SDURATION","dead")
  
  # At this moment "simulation_baseline_patients" is simply an empty data frame
  simulation_baseline_patients <- data.frame(matrix(vector(), 0, length(history_characteristics), dimnames=list(c(), history_characteristics)),stringsAsFactors=F)
  
  
  ##################################################
  ########## MAIN PART I: simulate events ##########
  ##################################################
  
  
  ### Regression coefficients for PSA -- once per PSA
  if(run_PSA_input == 1){} # end if regression coef PSA
  
  
  ### IMPORTANT: Set a random seed to be able to replicate the results.
  ### This first seed is used to draw the same pool of patients when the function is called multiple times.
  ### Decide later if this has to go inside or outside the function. It depends on whether this seed is going to be 
  ### always fixed (inside) or if it can change (outside - and make it random or user defined)
  
  set.seed(seed_input) 
  
  ### Select the sample used in the simulation: sample with replacement from the dataset.
  ### At this moment "all_baseline_patients" consists of only one patient with "dummy" values as risk factors
  simulation_baseline_patients[1:patient_size_input,history_characteristics[-1]] <- all_baseline_patients[sample(nrow(all_baseline_patients), patient_size_input, replace = TRUE), ]
  
  ### Create the simulation patient history table (for now is just empty)
  simulation_patients_history <- simulation_baseline_patients[FALSE,c(history_characteristics)]
  
  ### Choose the 1st patient
  patient_index <- 1
  
  ### Begin the loop on the simulation size (i.e. the number of patients you want to simulate)
  for(patient_index in 1:patient_size_input){
    
    # Print the patient index to know how advanced is the simulation.
    # Try to show this in the interface.
    if(run_PSA_input == 0){print(patient_index)}
    
    # Pick the current patient from those selected at baseline and set simulation ID.
    current_patient <- simulation_baseline_patients[patient_index,]
    current_patient$SIMID    <- patient_index
    
    # QUESTION: At this moment we only have one patient with "dummy" values for the risk factors. We will update here some of the
    #           risk factors but we should delete these later. Since we only have one patient, the same patient will be sampled multiple times
    #           in the current simulation.
    current_patient$AGE.DIAG <- 52 
    current_patient$YEAR     <- 8
    current_patient$CURR.AGE <- current_patient$AGE.DIAG + current_patient$YEAR
    
    current_patient$LDL <- 3*10
    current_patient$BMI <- 32
    current_patient$BMI1 <- 0
    current_patient$BMI3 <- 1
    current_patient$FEMALE <- 1
    current_patient$MMALB <- 1
    
    current_patient$eGFR <- 50/10
    current_patient$eGFR60less <- 50/10
    current_patient$SMOKER <- 1
    current_patient$HEART.R <- 80/10
    current_patient$WBC <- 6
    
    current_patient$AMP.HIST <- 0
    current_patient$MI.HIST <- 0
    current_patient$CHF.HIST <- 0
    
    # Initialise tracking variables for second events
    current_SMI_event <- 0
    
    # Save the characteristics to be used in the simulation history 
    # QUESTION: for the moment everything is saved (those that are stable too, but we could keep only those changing) 
    simulation_patients_history <- rbind(simulation_patients_history,current_patient[history_characteristics])
    
    
    ### Start the "timed" simulation (while loop = clock)   #
    current_sim_time   <- 1
    factor_for_seed <- 100 #test = 1 to reporduce results. I normally used 100 but can be anything large enough
    
    while(current_patient$dead==0){
      
      ### These seeds will be used to draw the annual event probabilities while the patient is alive 
      if(run_PSA_input == 0){
        set.seed(factor_for_seed*current_patient$SIMID + current_sim_time)
      }else{
        set.seed(factor_for_seed*seed_input*current_patient$SIMID + current_sim_time)
      }
      
      #####################################
      # Sample annual event probabilities #
      #####################################
      
      # Assumptions: 
      #              1. Multiple events can occur in one year.
      #              2. Occurrence of events only affects death probability in the year they occur.
      #              3. Occurrence of events affects other events probability in the year after they occur (through history variables).
      
      # HOWEVER, IN THE DESCRIPTION OF THE MODEL (UKPDS PAPER) IT SAYS "RANDOMLY ORDER AND RUN EVENT EQUATIONS" -->>
      # ORDER WOULD IMPLY SEQUENTIAL RUN OF EQUATIONS AND UPDATES COULD BE BASED ON CURRENT YEAR. NO FURTHER EXPLANATION TO THIS IS FOUND 
      
      
      ### MACROVASCULAR COMPLICATIONS ###
      
      # Heart Failure is Weibull. This can happen only once; that's why the if condition below is used. 
      if(current_patient$CHF.HIST == 0){
        current_CHF_prob  <- treatment_effect_input*annual_p_weibull(macrovascular_risk_equations$CHF,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
        current_CHF_event <- rbinom(1,1,current_CHF_prob) # Update current_patient$CHF.HIST after the year simulation is finished
      }
      
      # IHD is Weibull. This can happen only once; that's why the if condition below is used.
      if(current_patient$IHD.HIST == 0){
        current_IHD_prob  <- treatment_effect_input*annual_p_weibull(macrovascular_risk_equations$IHD,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
        current_patient$IHD.EVENT <- rbinom(1,1,current_IHD_prob) # Update current_patient$IHD.HIST after the year simulation is finished
      }
      
      
      # MI could be first or second. If it is first, then it is different for male and female. 
      # If it is second, then it is the same for both genders.
      # The variables MI.EVENT and MI.HIST do not distinguish between first and second. However, the model assumption
      # is that no more than 2 MI events are possible. this has to be tracked.
      if(current_patient$MI.HIST == 0){
        
        # If no history of MI, then it is first. Different for males and females
        if(current_patient$FEMALE == 1){
          # First MI for female Weibull
          current_FMIFEMALE_prob  <- treatment_effect_input*annual_p_weibull(macrovascular_risk_equations$FMIFEMALE,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
          current_patient$MI.EVENT <- rbinom(1,1,current_FMIFEMALE_prob) # Update current_patient$MI.HIST after the year simulation is finished
          
        }else{
          # First MI for male Exponential
          current_FMIMALE_prob  <- treatment_effect_input*annual_p_weibull(macrovascular_risk_equations$FMIMALE,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
          current_patient$MI.EVENT <- rbinom(1,1,current_FMIMALE_prob) # Update current_patient$MI.HIST after the year simulation is finished
        } # end if/else for gender
        
      }else{
        # Second MI is Exponential, regardless the gender. MI.HIST will remain == 1.
        # If a patient has a second MI, it is not  possible to have a 3rd. "current_SMI_event" keeps track of this. Not sure whether it has to be initialised. Seems to be 0, which is ok.
        if(current_SMI_event == 0){
          current_SMI_prob <- treatment_effect_input*annual_p_weibull(macrovascular_risk_equations$SMI,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
          current_patient$MI.EVENT <- rbinom(1,1,current_SMI_prob)
          current_SMI_event <- current_patient$MI.EVENT
        }
      }
      
      
      # STROKE could be first or second. 
      # The variables STROKE.EVENT and STROKE.HIST do not distinguish between first and second. However, the model assumption
      # is that no more than 2 STROKE events are possible. this has to be tracked.
      if(current_patient$STROKE.HIST == 0){
        # If no history of STROKE, then it is first and Weibull. 
        current_FSTROKE_prob <- treatment_effect_input*annual_p_weibull(macrovascular_risk_equations$FSTROKE,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
        current_patient$STROKE.EVENT <- rbinom(1,1,current_FSTROKE_prob) #Update current_patient$STROKE.HIST after the year
        
      }else{ 
        #Second STROKE is Weibull. Here current_patient$STROKE.HIST == 1 and remains like that.
        #If a patient has a second STROKE, it is not  possible to have a 3rd. "current_SSTROKE_event" keeps track of this. Not sure whether it has to be initialised. Seems to be 0, which is ok.
        if(current_SSTROKE_event == 0){
          current_SSTROKE_prob <- treatment_effect_input*annual_p_weibull(macrovascular_risk_equations$SSTROKE,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
          current_patient$STROKE.EVENT <- rbinom(1,1,current_SSTROKE_prob) 
          current_SSTROKE_event <- current_patient$STROKE.EVENT
          }
        }
      
      
      ### MICROVASCULAR COMPLICATIONS ###
      
      # BLINDNESS is Exponential. This is assumed to happen only once; that's why the if condition below is used. 
      if(current_patient$BLIND.HIST == 0){
        current_BLIND_prob <- treatment_effect_input*annual_p_weibull(microvascular_risk_equations$BLIND,current_patient %>% select(risk_factors_microvascular),current_patient$YEAR)$p       
        current_BLIND_event <- rbinom(1,1,current_BLIND_prob) #Update current_patient$BLIND.HIST at the end of the year. 
      }

      # ULCER is Logistic. This is assumed to happen only once; that's why the if condition below is used. 
      if(current_patient$ULCER.HIST == 0){
        current_ULCER_prob  <- treatment_effect_input*annual_p_logistic(microvascular_risk_equations$BLIND,current_patient %>% select(risk_factors_microvascular))$p
        current_ULCER_event <- rbinom(1,1,current_ULCER_prob) # Update current_patient$ULCER.HIST at the end of the year.
        }
      

      # AMPUTATION could be first or second. First amputation depends on ULCER history. 
      
      if(current_patient$AMP.HIST == 0){
        # If no history of AMPUTATION, then it is first and depends on ULCER.
        if(current_patient$ULCER.HIST == 0){
          # If no prior ULCER then it is Weibull
          current_FAMPNOULCER_prob  <- treatment_effect_input*annual_p_weibull(microvascular_risk_equations$FAMPNOULCER,current_patient %>% select(risk_factors_microvascular),current_patient$YEAR)$p       
          current_patient$AMP.EVENT <- rbinom(1,1,current_FAMPNOULCER_prob) #AMP.HIST to be updated at the end of the year and does not distinguishes between 1st and 2nd
        }else{
          # If prior ULCER then it is Exponential
          current_FAMPULCER_prob <- treatment_effect_input*annual_p_weibull(microvascular_risk_equations$FAMPULCER,current_patient %>% select(risk_factors_microvascular),current_patient$YEAR)$p       
          current_patient$AMP.EVENT <- rbinom(1,1,current_FAMPULCER_prob) #AMP.HIST to be updated at the end of the year
          }
        }else{ # if there is amputation history it can only be second and a 3rd one is not possible. 
          # Second amputation is exponential: AMP.HIST no need to be updated. But also a patient cannot have more than 2 amputations.
          if(current_AMP2_event == 0){
            current_SAMP_prob <- treatment_effect_input*annual_p_weibull(microvascular_risk_equations$SAMP,current_patient %>% select(risk_factors_microvascular),current_patient$YEAR)$p       
            current_patient$AMP2.EVENT <- rbinom(1,1,current_SAMP_prob) 
            current_AMP2_event <- current_patient$AMP2.EVENT
            }
          }
      
      
      # Renal failure is Exponential. This is assumed to happen only once; that's why the if condition below is used. 
      if(current_patient$RENAL.HIST == 0){
        current_RENALF_prob <- treatment_effect_input*annual_p_weibull(microvascular_risk_equations$RENALF,current_patient %>% select(risk_factors_microvascular),current_patient$YEAR)$p       
        current_patient$RENAL.EVENT <- rbinom(1,1,current_RENALF_prob) #current_patient$RENAL.HIST updated after the year. Check here: I got really large values for the probability so have a second look
        }
      
      
      ### DEATH ###
      
      # There are four equations for death, depending on the events occurred in the current year and the history of previous events. 
      # Assumption (already mentioned above): occurrence of events only affects death probability in the year they occur.
      # The four equations are for:
      #                             1. Years with no history of previous events and no events in the current year
      #                             2. First year of events (so no previous history) excluding blindness or ulcer
      #                             3. Years with history of previous events but no events in the current year
      #                             4. Subsequent years (so there is previous history) of events excluding blindness or ulcer
      
      # Therefore, we need to define variables to determine what equation should be used.  
      
      # If any event happened in the current year, it should be captured with the following variable: 
      current_year_event <- sum(current_CHF_event, 
                                current_patient$IHD.EVENT, 
                                current_patient$MI.EVENT, # could be 1st or 2nd, no distinction
                                current_patient$STROKE.EVENT, # could be 1st or 2nd, no distinction
                                current_BLIND_event, 
                                current_ULCER_event, 
                                current_patient$AMP.EVENT, 
                                current_patient$AMP2.EVENT,
                                current_patient$RENAL.EVENT)
      
      # If any event except blindness and ulcer happened in the current year, it should be captured with the following variable:
      current_year_event_no_blind_no_ulcer <- sum(current_CHF_event, 
                                                  current_patient$IHD.EVENT, 
                                                  current_patient$MI.EVENT, 
                                                  current_patient$STROKE.EVENT, 
                                                  current_patient$AMP.EVENT, 
                                                  current_patient$AMP2.EVENT,
                                                  current_patient$RENAL.EVENT)
      
      # If current_year_event - current_year_event_no_blind_no_ulcer > 0 it means that either blindness or ulcer occurred in the current year.
      
      # .HIST variables are not updated for the current year. But the ones from the previous year are captured in this variable:
      current_hist  <- sum(current_patient$CHF.HIST, 
                           current_patient$IHD.HIST, 
                           current_patient$MI.HIST, 
                           current_patient$STROKE.HIST, 
                           current_patient$BLIND.HIST, 
                           current_patient$ULCER.HIST,
                           current_patient$AMP.HIST, 
                           current_patient$RENAL.HIST)
      
      
      # The four equations for death are then the following: 
      
      # 1. If no history of previous events and no events in the current year, then gompertz distirbution
      if(current_year_event == 0 & current_hist == 0){
        current_DEATH_prob <- annual_p_gompertz(mortality_risk_equations$DEATHNOHIST, current_patient %>% select(risk_factors_mortality),current_patient$AGE.DIAG + current_patient$YEAR)$p       
        }
      
      #2. First year of events (so no previous history) excluding blindness or ulcer, then logistic distribution
      if(current_year_event_no_blind_no_ulcer == 1 & current_hist == 0){
        current_DEATH_prob <- annual_p_logistic(mortality_risk_equations$DEATH1YEVENT, current_patient %>% select(risk_factors_mortality))$p       
        }
      
      #3. Years with history of previous events but no events in the current year, then gompertz distribution
      if(current_year_event == 0 & current_hist == 1){
        current_DEATH_prob <- annual_p_gompertz(mortality_risk_equations$DEATHHISTNOEVENT, current_patient %>% select(risk_factors_mortality),current_patient$AGE.DIAG + current_patient$YEAR)$p       
        }
      
      #4. Subsequent years (so there is previous history) of events excluding blindness or ulcer, then logistic distribution
      if(current_year_event_no_blind_no_ulcer == 1 & current_hist == 1){
        current_DEATH_prob  <- annual_p_logistic(mortality_risk_equations$DEATHYSEVENT, current_patient %>% select(risk_factors_mortality))$p       
        }
      
      current_DEATH_event <- rbinom(1,1,current_DEATH_prob) 
      current_patient$dead <- current_DEATH_event
      
      
      ####################################
      # Update patient characteristics   #
      ####################################
      
      # We first copy all the previous characteristics
      current_patient_update <- current_patient # Do we need this? Or can we update directly current_patient?
      
      # Update first the history characteristics 
      current_patient_update$CHF.HIST    <- ifelse(current_CHF_event + current_patient$CHF.HIST == 0, 0, 1)
      current_patient_update$IHD.HIST    <- ifelse(current_patient$IHD.EVENT + current_patient$IHD.HIST == 0, 0, 1)
      current_patient_update$MI.HIST     <- ifelse(current_patient$MI.EVENT + current_patient$MI.HIST == 0,0,1) 
      current_patient_update$STROKE.HIST <- ifelse(current_patient$STROKE.EVENT + current_patient$STROKE.HIST == 0, 0, 1)
      current_patient_update$BLIND.HIST  <- ifelse(current_BLIND_event + current_patient$BLIND.HIST == 0, 0, 1)
      current_patient_update$ULCER.HIST  <- ifelse(current_ULCER_event + current_patient$ULCER.HIST == 0, 0, 1)
      current_patient_update$AMP.HIST    <- ifelse(current_patient$AMP.EVENT + current_patient$AMP2.EVENT + current_patient$AMP.HIST == 0, 0, 1)
      current_patient_update$RENAL.HIST  <- ifelse(current_patient$RENAL.EVENT + current_patient$RENAL.HIST == 0, 0, 1)
      
      
      # Update risk factors  
      
      ### QUESTION: What remains unclear is whether all  risk factors change with time... and how...
      ### Decide what factors are assumed to be stable and which ones to change with time...
      ### For those that change, we need equations for predicting annual change based on other patient chartacteristics.
      ### WORK IN PROGRESS
      
      #current_patient_update$ATFIB <- current_patient$ATFIB + ?
      current_patient_update$BMI <- current_patient$BMI
      #based on the above BMI, we should update BMI1 and BMI3
      current_patient_update$CURR.AGE <- current_patient$CURR.AGE + 1 
      current_patient_update$eGFR <- current_patient$eGFR
      # based on the above, update eGFR60less and eGFR60more?
      current_patient_update$HAEM  <- current_patient$HAEM
      current_patient_update$HbA1c <- current_patient$HbA1c  
      current_patient_update$HDL   <- current_patient$HDL
      current_patient_update$HEART.R <- current_patient$HEART.R
      current_patient_update$LDL <- current_patient$LDL
      # based on the above update LDL35more
      current_patient_update$MMALB <- current_patient$MMALB
      current_patient_update$PVD <- current_patient$PVD
      current_patient_update$SBP <- current_patient$SBP
      current_patient_update$SMOKER <- current_patient$SMOKER
      current_patient_update$WBC <- current_patient$WBC
      current_patient_update$YEAR <- current_patient$YEAR + 1
      current_patient_update$SDURATION <- current_patient$SDURATION + 1
      
      # Force death after 100 years?
      if(current_patient_update$CURR.AGE > 100){current_patient_update$dead <- 1}
        
      # Force death if updated continuous risk factors take unfeasible value (e.g. too high hba1c)?
      # Consider logical bouns for update continuous attributes, e.g. > 0
      
      ### When all characteristics have been updated, we add these to the patient history
      simulation_patients_history <- rbind(simulation_patients_history,current_patient_update[history_characteristics])
      
      ### And update current patient and go up to while loop
      current_patient <- current_patient_update
      
      ### All the _event and .EVENT variables have to be reset to 0 because for the next year it only counts .HIST
      
      current_CHF_event   <- 0
      current_BLIND_event <- 0 
      current_ULCER_event <- 0
      
      current_patient$AMP.EVENT    <- 0
      current_patient$AMP2.EVENT   <- 0
      current_patient$MI.EVENT     <- 0
      current_patient$IHD.EVENT    <- 0
      current_patient$RENAL.EVENT  <- 0
      current_patient$STROKE.EVENT <- 0
      
      current_sim_time <- current_sim_time + 1 # not sure if we need this. is it the variable .YEAR?
      
    } #end while loop and move to another patient
    
    patient_index <- patient_index + 1
    
  } #end for loop in number of patients
  
  

  #################################################################
  ########## MAIN PART III: Calculate aggregated results ##########
  #################################################################
  
  # To be done: at this moment, only clinical outcomes (e.g. event rates and life expectancy can be calculated).
  
  ### Life expectancy
  mean_life_expectancy <- round(mean(simulation_patients_history[which(simulation_patients_history$dead==1),"SDURATION"]),4)
  
  ### Events calculated differently depending on how were defined: .EVENT or .HIST
  
  # Total number of renal events per patient during lifetime
  cum_RENAL <- aggregate(simulation_patients_history$RENAL.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  # Rate per year
  mean_RENAL_rate <- round(mean(cum_RENAL$x)/mean_life_expectancy, 4)

  
  ### Return model outcomes: at this moment only clinical history is saved and returned
  return(list(simulation_patients_history=simulation_patients_history,
              mean_life_expectancy=mean_life_expectancy))
  
} #end SMDMII_model_simulation function

# To run the model simply call the model function with the appropriate inputs in the correct order. For example,
# the line below will run the model for 500 patients, deterministically, without treatment effects and with a 
# random seed equal to 177.


sim_results <- SMDMII_model_simulation(20, #patient_size_input
                                       0, #run_PSA_input, 0 == no PSA
                                       1, #treatment_effect_input --> note there might be probably one input for each effect modifier: e.g. on CHF, MI, stroke, etc...
                                       177) #seed_input


View(sim_results$simulation_patients_history)
sim_results$mean_life_expectancy

### To discuss:

### MODEL UPDATES: 1. ASSUMPTIONS
###                2. AGGREGATED RESULTS (in progress --> discuss what we want to show: e.g. distinguish between first and second events?)  

### NEXT STEPS:
###                     1. HOW TO MODEL BASELINE POPULATION
###                     2. HOW TO UPDATE RISK FACTORS
