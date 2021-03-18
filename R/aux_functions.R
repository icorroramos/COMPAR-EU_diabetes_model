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
# 6. BMI: Body mass index (m/kg^2) measured continuously. HR per unit increase in BMI
# 7. eGFR: Estimated glomerular filtration rate (ml/min/1.73m^2) from modification of diet in renal disease (MDRD) formula.
#          Continuously measured and further DIVIDED by 10. HR per 10 ml/min/1.73m^2 increase 
# 8. eGFR60less: Same as eGFR. Continuous spline (knot at 60) and further DIVIDED by 10. HR per 10 ml/min/1.73m^2 increase if < 60
# 9. HbA1c: percentage HbA1c measured continuously. HR per 1% increase in HbA1c
# 10. HDL: high density lipoprotein cholesterol (mmol/l) measured continuously and MULTIPLIED by 10. HR per 0.1 mg/dL increase in HDL
# 11. LDL: low density lipoprotein cholesterol (mmol/l) measured continuously and MULTIPLIED by 10. HR per 0.1 mg/dL increase in LDL
# 12. LDL35more: same as LDL. Continuous spline (knot at 3.5) and MULTIPLIED by 10. HR per 0.1 mg/dL increase in LDL if > 3.5 mg/dL
# 13. MMALB: presence of micro- or macro-albuminuria. 1 == urine albumin >= 50 mg/l, 0 == otherwise. HR referent == no albuminuria
# 14. PVD: 1 == peripheral vascular disease, 0 == otherwise. Defined from presence of intermittent claudication or ankle brachial
#          pressure index < 0.9. HR referent == no presence/evidence of PVD
# 15. SBP: systolic blood pressure (mm Hg) measured continuously and further DIVIDED by 10. HR per 10 mm Hg increase in SBP --> 10 or 100?
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
# first MI for males: exponential distribution ro = 1 (2nd element in vector) --> explanation below in "annual_p_weibull" function
FMIMALE <- c(-8.791, 1, -0.83, 0.045, 0 , 0.279, 0, 0, 0, 0, 0.108, -0.049, 0.023, 0, 0.203, 0.340, 0.046, 0.277, 0.026, 0.743, 0.814, 0.846, 0.448, 0) 
# first MI for females
FMIFEMALE <- c(-8.708, 1.376, -1.684, 0.041, 0, 0, 0, 0, 0, -0.28, 0.078, 0, 0, 0.035, 0.277, 0.469, 0.056, 0.344, 0.07, 0, 0.853, 0.876, 0, 0)
# second MI: exponential distribution ro = 1 (2nd element in vector)
SMI <- c(-4.179, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.021, 0, 0.344, 0, 0, 0, 0, 0, 0, 0, 0, 0) 
# first stroke
FSTROKE <- c(-13.053, 1.466, 0, 0.066, -0.42, 0, 1.476, 0, 0, -0.190, 0.092, 0, 0.016, 0, 0.420, 0, 0.170, 0.331, 0.040, 1.090, 0, 0.481, 0, 0)
# second stroke
SSTROKE <- c(-9.431, 1.956, 0, 0.046, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.537, 0, 0, 0.656, 0, 0, 0, 0, 0, 0)

# Below we simply create a table (R data frame) with all the coefficients of the regression equations used to predict macrovascular complications.
# When we define the risk equations below in the code, these will read the coefficients from this table in order to predict the annual 
# probability of experiencing a mcrovascular complication.
macrovascular_risk_equations <- data.frame(CHF, IHD, FMIMALE, FMIFEMALE, SMI, FSTROKE, SSTROKE, row.names = parameters_macrovascular)
# macrovascular_risk_equations


# Microvascular complications: UKPDS paper ESM Table 5. Microvascular complications include blindness, diabetic ulcer amputation and renal failure.
# The risk factors (the ones that are not defined above for macrovascular complications) used to predict microvascular complications 
# are the following (notation used in the model, definition from UKPDS paper ESM Table 2):

# 23. eGFR60more: Same as eGFR. Continuous spline (knot at 60) and further DIVIDED by 10. HR per 10 ml/min/1.73m^2 increase if > 60
# 24. HAEM: haemoglobin g/dL measured continuously. HR per 1 g/dL increase
# 25. HEART.R: heart rate (beats per minute) determined from inspiration/expiration RR on ECG. Continuously measured and further DIVIDED by 10.
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
# QUESTION: ro is 2nd element. ULCER is logistic so don't know yet if ro = 1  or 0

# First amputation, no prior ulcer history
FAMPNOULCER <- c(-14.844, 2.067, 0, 0.023, -0.445, 1.088, 0, 0, 0, 0, 0.248, -0.059, 0.098, 0, 0.602, 1.010, 0.086, 0.040, 0, 0, 0, 0, 1.299)
# First amputation, prior ulcer history
FAMPULCER <- c(-0.811, 1, 0, -0.065, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.769, 0, 0, 0, 0, 0, 0, 0)
# Second amputation
SAMP <- c(-3.455, 1, rep(0,8), 0.127, rep(0,12))
# Renal failure
RENALF <- c(3.549, 1, 0.686, -0.029, -0.869, 0, -0.054, -1.031, -0.487, -0.268, 0, 0, 0, 0.027, 1.373, 0, 0.085, 0.029, 1.108, 0.732, 0, 0, 0)

# Below we simply create a table (R data frame) with all the coefficients of the regression equations used to predict microvascular complications.
# When we define the risk equations below in the code, these will read the coefficients from this table in order to predict the annual 
# probability of experiencing a microvascular complication.
microvascular_risk_equations <- data.frame(BLIND, ULCER, FAMPNOULCER, FAMPULCER, SAMP, RENALF, row.names = parameters_microvascular)

# Risk of death: UKPDS paper ESM Table 6. Four equations depending on events and history. The risk factors (the ones that are not 
# defined above for macro/microvascular complications) are the following (notation used in the model, definition from UKPDS paper ESM Table 2):

# 27. YEAR: duration of diabetes in years measured continuously. HR per year increase in duration of diabetes. 
#           Note it increases + 1 after each simulated year.
# 28. BMI1: 1 == BMI < 18.5 m/kg^2, 0 == otherwise. HR referent 18.5 m/kg^2 <= BMI < 25 m/kg^2. It depends thus on the variable BMI.
# 29. BMI3: 1 == BMI >= 25 m/kg^2, 0 == otherwise. HR referent 18.5 m/kg^2 <= BMI < 25 m/kg^2. It depends thus on the variable BMI.
# 30. CURR.AGE: current age in years measured continuously. HR per year increase in current age. Note CURR.AGE = AGE.DIAG + YEAR at the 
#               beginning of the simulation. Also after each year in the simulation if YEAR is properly updated.
# 31. AMP1.EVENT: 1 == first amputation, 0 == otherwise. HR referent == no first amputation event
# 32. AMP2.EVENT: 1 == second amputation, 0 == otherwise. HR referent == no second amputation event
# 33. IHD.EVENT: 1 == IHD, 0 == otherwise. HR referent == no IHD event
# 34. MI.EVENT: 1 == MI, 0 == otherwise. HR referent == no MI event
# 35. MI.HIST: 1 == history of MI, 0 == otherwise. HR referent == no prior MI
# 36. RENAL.EVENT: 1 == renal failure event, 0 == otherwise. HR referent == no renal failure event
# 37. RENAL.HIST: 1 == history of renal failure, 0 == otherwise. HR referent == no prior renal failure
# 38. STROKE.EVENT: 1 == stroke, 0 == otherwise. HR referent == no stroke event


# QUESTION: pay extra attention to the data transformation. Otherwise, the equations will not make sense!

# The vector below contains the names of all risk factors used to predict death in the model. 
risk_factors_mortality <- c("FEMALE", "INDIAN", "YEAR", "ATFIB", "BMI1", "BMI3", "CURR.AGE", "HDL", "HEART.R", 
                            "MMALB", "PVD", "SMOKER", "WBC", "AMP1.EVENT", "AMP.HIST", "AMP2.EVENT", "CHF.HIST", 
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


# Below we simply create a table (R data frame) with all the coefficients of the regression equations used to predict death.
# When we define the risk equations below in the code, these will read the coefficients from this table in order to predict the annual 
# probability of death.
mortality_risk_equations <- data.frame(DEATHNOHIST, DEATH1YEVENT, DEATHHISTNOEVENT, DEATHYSEVENT, row.names = parameters_mortality)


# Below we consider all the patient characteristics (38) that "define" a patient in this model. 
risk_factors_simulation <- unique(sort(c("CHF.EVENT", "BLIND.EVENT", "ULCER.EVENT",
                                         risk_factors_macrovascular, 
                                         risk_factors_microvascular, 
                                         risk_factors_mortality)))
# Even though not defined as such in the UKPDS equations, we have considered 
# the variables CHF.EVENT, BLIND.EVENT, ULCER.EVENT XXXX to keep the code more consistent.

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
  
  return(list(#H_t1 = H_t1, # no need to return the H's, keep them for now just for validation purposes
              #H_t2 = H_t2,
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
  
  return(list(#H_t1 = H_t1, # no need to return the H's, keep them for now just for validation purposes
              #H_t2 = H_t2,
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

# Small validation exercise. Delete when OK.

# Reproduce examples in statistical appendix of UKPDS paper:
# annual_p_gompertz is a function of (regression_coefficients_input, risk_factors_input, current_age_input)
# risk factors inputs for mortality are ("FEMALE", "INDIAN", "YEAR", "ATFIB", "BMI1", "BMI3", "CURR.AGE", "HDL", "HEART.R", 
# "MMALB", "PVD", "SMOKER", "WBC", "AMP1.EVENT", "AMP.HIST", "AMP2.EVENT", "CHF.HIST", 
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


#################################
# PART II - Define global lists #
#################################

# .EVENT variable names

event_vars <- c("CHF.EVENT", "BLIND.EVENT", "ULCER.EVENT", "AMP1.EVENT", "AMP2.EVENT", 
                "MI.EVENT", "IHD.EVENT", "RENAL.EVENT", "STROKE.EVENT")

# Indicate the patient characteristics that we will save during the simulation. 
# Besides the risk factors, we have a simulation ID, duration of diabetes 
# and an indicator variable for dead. risk_factors_simulation is defined above.

# Thus, the vector below contains the names of all factors used to predict informal care in the model.
risk_factors_informal   <- c("FEMALE", "CURR.AGE.SCALE.INF", "CHF.HIST", "STROKE.HIST", "RENAL.HIST")
risk_factors_prod       <- c("FEMALE", "CURR.AGE.SCALE.PROD", "CHF.HIST", "STROKE.HIST", "RENAL.HIST")
risk_factors_employment <- c("CURR.AGE", "CURR.AGE.2", "FEMALE")

history_characteristics <- c("SIMID", 
                             unique(sort(c(risk_factors_simulation, 
                                           risk_factors_informal, 
                                           risk_factors_prod,
                                           risk_factors_employment))), 
                             "INF.CARE", "EMPLOYED","PROD.LOSS", #Added 29/08/2020
                             "SDURATION",
                             "dead")


#############################################################################
########## PART III - INFORMAL CARE AND PRODUCTIVITY EQUATIONS ##############
#############################################################################

# The equation used to predict informal care use in the Netherlands is the following:

# And now, below we have the regression coefficients as reported above.
# CHECK that female coefficient is actually for females

# UPDATE EQUATIONS WITH MOST RECENT VERSION: 10/12/2020
# Note there is no country specific coefficient now. That's why the last coefficient = 0. 
# Otherwise the code would not work.
informal_care_coef <- c(-1.5355, 0.5036, 0.6852, 0.4231, 0.2245, 0.5498, 0) # Assumed primary analysis central countries for UK
employment_coef    <- c(2.4283, 0.1273, -0.00252, -0.0737, 0) # Assumed primary analysis central countries for UK
prod_costs_coef    <- c(-1.4747, 0.01047, 0.5112, 0.08407, 0.1881, 1.3026, 0) # Assumed primary analysis central countries for UK

# Besides the risk factors, there are two other parameters used in the UKPDS equations called "lambda" and "ro". We added them below.
# Now we have the names of all the parameters (lambda, ro and risk factors) used to predict informal risks in the model.
parameters_informal   <- c("lambda", risk_factors_informal, "UK")
parameters_prod       <- c("lambda", risk_factors_prod, "UK")
parameters_employment <- c("lambda", risk_factors_employment, "UK")

# Below we simply create a table (R data frame) with all the coefficients of the regression equation to predict informal care
informal_care_equations <- data.frame(informal_care_coef,row.names = parameters_informal)
prod_cost_equations     <- data.frame(prod_costs_coef, row.names = parameters_prod)
employment_equations    <- data.frame(employment_coef, row.names = parameters_employment)

# Notation: as a general rule, names for function input parameters end with "_input".

# Note that age should be scaled before passing it to the function to calculate the probability of receiving informal care! 
annual_p_bernoulli <- function(regression_coefficents_input, risk_factors_input){
  risk_factors_input <- as.numeric(risk_factors_input) # delete if not needed
  # Note: mind the order of the regression coefficients. 
  log.oods <- sum(regression_coefficents_input*c(1,risk_factors_input,1))
  p <- exp(log.oods)/(1+exp(log.oods))
  return(list(p = p))
}
