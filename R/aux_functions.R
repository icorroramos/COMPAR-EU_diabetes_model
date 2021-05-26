#######################################################
########## PART I - UKPDS RISK EQUATIONS ##############
#######################################################

# This model makes use of previously estimated regression equations to predict the risk of experiencing diabetes-related 
# complications. These regression equations were sourced from UKPDS 82 and can be found in the paper by Hayes et al. 2013:

# UKPDS outcomes model 2: a new version of a model to simulate lifetime health outcomes of patients with type 2 diabetes mellitus
# using data from the 30 year United Kingdom Prospective Diabetes Study: UKPDS 82.
# Hayes AJ, Leal J, Gray AM, Holman RR, Clarke PM. Diabetologia. 2013 Sep;56(9):1925-33. Epub 2013 Jun 22.
# DOI: 10.1007/s00125-013-2940-y


### COMPLICATIONS ### 

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

# Multiple events are possible in the same year but only the number of times specified above. The risk factors and the regression 
# coefficients for each complication were taken from the UKPDS paper as explained below.

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
# Even though not defined as such in the UKPDS equations, we have considered the .EVENT variables to keep the code more consistent.

### UKPDS RISK FUNCTIONS ### 

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


### BACKGROUND MORTALITY

background_DEATH_prob <- read.csv("input/UK/nationallifetables3yearuk.csv")


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
