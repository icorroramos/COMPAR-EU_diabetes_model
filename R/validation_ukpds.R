#########################################################
########## VALIDATION UKPDS RISK EQUATIONS ##############
#########################################################

# Patient characteristics
validation_patient <- read.csv("input/UK/baseline_characteristics_UK.csv", sep=",")

retirement_age_input <- 65

validation_patient$AGE.DIAG <- validation_patient$CURR.AGE - validation_patient$YEAR 
validation_patient$CURR.AGE.SCALE.INF  <- (validation_patient$CURR.AGE - 72.5474088)/10.4626624 
validation_patient$CURR.AGE.SCALE.PROD <- (validation_patient$CURR.AGE - 60.2737989)/6.1177269 
validation_patient$CURR.AGE.2 <- (validation_patient$CURR.AGE)^2
validation_patient$INF.CARE   <- 0
ifelse(validation_patient$CURR.AGE >= retirement_age_input, validation_patient$EMPLOYED <- 0,
       {baseline_employed_prob <- apply(validation_patient %>% select(risk_factors_employment), 1, function(x) annual_p_bernoulli(employment_equations$employment_coef,x)$p)
       validation_patient$EMPLOYED <- unlist(lapply(baseline_employed_prob, function(x) rbinom(1,1,x))) #EMPLOYED = yes/no
       })
validation_patient$PROD.LOSS <- 0 
validation_patient$BMI1 <- if_else(validation_patient$BMI < 18.5, 1, 0)
validation_patient$BMI3 <- if_else(validation_patient$BMI >= 25, 1, 0)

#validation_patient[event_vars] <- 0 

validation_patient$eGFR       <- validation_patient$eGFR/10 
validation_patient$eGFR60more <- if_else(validation_patient$eGFR >= 6, validation_patient$eGFR, 0)
validation_patient$eGFR60less <- if_else(validation_patient$eGFR <  6, validation_patient$eGFR, 0)
validation_patient$HDL <- validation_patient$HDL*10
validation_patient$HEART.R <- validation_patient$HEART.R/10 
validation_patient$LDL       <- validation_patient$LDL*10
validation_patient$LDL35more <- if_else(validation_patient$LDL >= 35, validation_patient$LDL, 0)
validation_patient$SBP <- validation_patient$SBP/10

validation_patient


# Complications

# Macrovascular complications: UKPDS paper ESM Table 4. Macrovascular complications include CHF, IHD, MI and stroke.
# The risk factors used to predict macrovascular complications are read from the following file (UKPDS paper ESM Table 2):
macrovascular_risk_equations <- read.csv("input/UKPDS_macrovascular_coef.csv")

# The vector below contains the names of all risk factors used to predict macrovascular complications in the model.
# this is needed below and in the simulation function.
risk_factors_macrovascular <- macrovascular_risk_equations$X[-(1:2)] 

# Microvascular complications: UKPDS paper ESM Table 5. Microvascular complications include blindness, diabetic ulcer amputation and renal failure.
# The risk factors used to predict microvascular complications are read from the following file (UKPDS paper ESM Table 5):
microvascular_risk_equations <- read.csv("input/UKPDS_microvascular_coef.csv")

# The vector below contains the names of all risk factors used to predict microvascular complications in the model.
# this is needed below and in the simulation function.
risk_factors_microvascular <- microvascular_risk_equations$X[-(1:2)]
  
# Risk of death: UKPDS paper ESM Table 6. Four equations depending on events and history:
mortality_risk_equations <- read.csv("input/UKPDS_mortality_coef.csv")

# The vector below contains the names of all risk factors used to predict death in the model. 
risk_factors_mortality <- mortality_risk_equations$X[-(1:2)]

# Below we consider all the patient characteristics (38) that "define" a patient in this model. 
risk_factors_simulation <- unique(sort(c("CHF.EVENT", "BLIND.EVENT", "ULCER.EVENT",
                                         risk_factors_macrovascular, 
                                         risk_factors_microvascular, 
                                         risk_factors_mortality)))

### UKPDS RISK FUNCTIONS ### 

# Note that, with the notation used in the UKPDS paper, if ro = 1 then 
# weibull distribution == exponential distribution. Thus, the same function "annual_p_weibull" can be used when an exponential
# distribution is assumed.

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


# Validation of the Weibull events 

# CHF females
annual_p_weibull(macrovascular_risk_equations$CHF,validation_patient[1,] %>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
# [1] 0.00657647 
# CHF males
annual_p_weibull(macrovascular_risk_equations$CHF,validation_patient[2,] %>% select(risk_factors_macrovascular),validation_patient[2,"YEAR"])$p
# [1] 0.00657647

# Is the predicted probability plausible? --> Gimon
# Is a difference between gender expected? NO


# IHD females
annual_p_weibull(macrovascular_risk_equations$IHD,validation_patient[1,] %>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
# [1] 0.0009983311
# IHD males
annual_p_weibull(macrovascular_risk_equations$IHD,validation_patient[2,] %>% select(risk_factors_macrovascular),validation_patient[2,"YEAR"])$p
# [1] 0.001698897

# Is the predicted probability plausible? --> Gimon
# Is a difference between gender expected? Yes, negative coefficient for females.

# First MI males (for females there is a different equation)
annual_p_weibull(macrovascular_risk_equations$FMIMALE,validation_patient[2,] %>% select(risk_factors_macrovascular),validation_patient[2,"YEAR"])$p
# [1] 0.002285789
# Is the predicted probability plausible? --> Gimon

# First MI females (for males there is a different equation)
annual_p_weibull(macrovascular_risk_equations$FMIFEMALE,validation_patient[1,] %>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
# [1] 0.02571863
# Is the predicted probability plausible? --> Gimon. It looks quite large compared to males. Does it make sense?

# Second MI females
annual_p_weibull(macrovascular_risk_equations$SMI,validation_patient[1,] %>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
# [1] 0.01946934
# Second MI males
annual_p_weibull(macrovascular_risk_equations$SMI,validation_patient[2,] %>% select(risk_factors_macrovascular),validation_patient[2,"YEAR"])$p
# [1] 0.01946934

# Is the predicted probability plausible? --> Gimon
# Is a difference between gender expected? No, gender is not a predictor.

# First stroke females
annual_p_weibull(macrovascular_risk_equations$FSTROKE,validation_patient[1,] %>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
# [1] 0.005297622
# First stroke males
annual_p_weibull(macrovascular_risk_equations$FSTROKE,validation_patient[2,] %>% select(risk_factors_macrovascular),validation_patient[2,"YEAR"])$p
# [1] 0.00805162

# Is the predicted probability plausible? --> Gimon
# Is a difference between gender expected? Yes, negative coefficient for females.

# Second stroke females
annual_p_weibull(macrovascular_risk_equations$SSTROKE,validation_patient[1,] %>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
# [1] 0.01181782
# Second stroke males
annual_p_weibull(macrovascular_risk_equations$SSTROKE,validation_patient[2,] %>% select(risk_factors_macrovascular),validation_patient[2,"YEAR"])$p
# [1] 0.01181782

# Is the predicted probability plausible? --> Gimon
# Is a difference between gender expected? No, gender is not a predictor.


# Blindness females
annual_p_weibull(microvascular_risk_equations$BLIND,validation_patient[1,] %>% select(risk_factors_microvascular),validation_patient[1,"YEAR"])$p
# [1] 0.002737117
# Blindness males
annual_p_weibull(microvascular_risk_equations$BLIND,validation_patient[2,] %>% select(risk_factors_microvascular),validation_patient[2,"YEAR"])$p
# [1] 0.002737117

# Is the predicted probability plausible? --> Gimon
# Is a difference between gender expected? No, gender is not a predictor.


# First amputation no ulcer females
annual_p_weibull(microvascular_risk_equations$FAMPNOULCER,validation_patient[1,] %>% select(risk_factors_microvascular),validation_patient[1,"YEAR"])$p
# [1] 0.0002014471
# First amputation no ulcer males
annual_p_weibull(microvascular_risk_equations$FAMPNOULCER,validation_patient[2,] %>% select(risk_factors_microvascular),validation_patient[2,"YEAR"])$p
# [1] 0.0003143385

# Is the predicted probability plausible? --> Gimon
# Is a difference between gender expected? Yes, negative coefficient for females.


# First amputation ulcer females
annual_p_weibull(microvascular_risk_equations$FAMPULCER,validation_patient[1,] %>% select(risk_factors_microvascular),validation_patient[1,"YEAR"])$p
# [1] 0.01874466
# First amputation ulcer males
annual_p_weibull(microvascular_risk_equations$FAMPULCER,validation_patient[2,] %>% select(risk_factors_microvascular),validation_patient[2,"YEAR"])$p
# [1] 0.01874466

# Is the predicted probability plausible? --> Gimon
# Is a difference between gender expected? No, gender is not a predictor.










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




# ---> Check finally whether in the simulation code, these functions are indeed properly called with the right parameters.
