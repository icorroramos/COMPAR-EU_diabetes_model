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
background_DEATH_prob <- read.csv(paste0("input/", country.id, "/background_mortality.csv"))


######################################################
# PART II - Define global input parameters and lists #
######################################################

# Discount factors

discount_factors_all_countries <- data.frame(
  UK = c(0.035, 0.035), #(discount factor costs , discount factor effects)
  NL = c(0.04, 0.015),
  DE = c(0.03, 0.03),
  ES = c(0.03, 0.03),
  GR = c(0.03, 0.03)
)

discount_factors <- discount_factors_all_countries[, country.id]

# Retirement ages

retirement_ages_all_countries <- data.frame(
  UK = 65,
  NL = 66,
  DE = 66,
  ES = 65,
  GR = 62
)

retirement_age_input <- as.numeric(retirement_ages_all_countries[country.id])

# Patient characteristics. 
baseline_characteristics <- read.csv(paste0("input/", country.id, "/baseline_characteristics.csv"), sep=",")

# Direct costs of diabetes-related complications for UK are age-gender dependent.
male_cost_inputs   <- read.csv(paste0("input/", country.id, "/Event_cost_male_2020.csv"), sep=",")
female_cost_inputs <- read.csv(paste0("input/", country.id, "/Event_cost_female_2020.csv"), sep=",")

# When a societal perspective is adopted, we also have future costs. This are obtained from the PAID online tool.
future_medical_cost_inputs    <- read.csv(paste0("input/", country.id, "/Future_medical_costs_2020.csv"), sep=",") # Ingelin: 15/12/2020
future_nonmedical_cost_inputs <- read.csv(paste0("input/", country.id, "/Nonmedical_futurecosts_data.csv"), sep=",") # UK costs updated 11/12/2020 - corrected version from Hamraz

# Informal care costs
inf_care_cost_all_countries <- data.frame(
  UK = 22.79,
  NL = 13.47,
  DE = 13,
  ES = 11.16,
  GR = 7.36
)

inf_care_hour_cost <- as.numeric(inf_care_cost_all_countries[country.id])

inf_care_hours_all_countires <- data.frame(
  UK = 1.25318,
  NL = 1.25318,
  DE = 1.25318,
  ES = 2.26803,
  GR = 2.26803
)

inf_care_hours_input <- as.numeric(inf_care_hours_all_countires[country.id])

inf_care_age_scale_all_countries <- data.frame(
  UK = c(72.5474088, 10.4626624), #(mean, sd)
  NL = c(72.5474088, 10.4626624),
  DE = c(72.5474088, 10.4626624),
  ES = c(76.0286769, 9.6290315),
  GR = c(76.0286769, 9.6290315)
)

inf_care_age_scale_input <- inf_care_age_scale_all_countries[, country.id]

# Productivity costs
worked_hours_all_countries <- data.frame(
  UK = c(28, 38), #(female value, male value)
  NL = c(28, 38),
  DE = c(28, 38),
  ES = c(30, 37.5),
  GR = c(30, 37.5)
)

worked_hours_input <- worked_hours_all_countries[, country.id]

working_days_lost_all_countries <- data.frame(
  UK = c(14, 27.5), #((diabetes value, diabetes & CHF/stroke value)
  NL = c(14, 27.5),
  DE = c(14, 27.5),
  ES = c(10, 20),
  GR = c(10, 20)
)

working_days_lost_input <- working_days_lost_all_countries[, country.id]

cost_hour_sick_all_countries <- data.frame(
  UK = 29, #EUR --> 25.781, #GBP
  NL = 36.8, #EUR
  DE = 36.6, #EUR
  ES = 22.8, #EUR
  GR = 16.9 #EUR
)

cost_hour_sick_input <- as.numeric(cost_hour_sick_all_countries[country.id])

friction_period_all_countries <- data.frame(
  UK = 82.18125,
  NL = 85,
  DE = 69,
  ES = 75,
  GR = 98.6175
)

friction_period_input <- as.numeric(friction_period_all_countries[country.id])

prod.loss_age_scale_all_countries <- data.frame(
  UK = c(60.2737989, 60.2737989), #(mean, sd)
  NL = c(60.2737989, 60.2737989),
  DE = c(60.2737989, 60.2737989),
  ES = c(62.5992071, 6.6265962),
  GR = c(62.5992071, 6.6265962)
)

prod_loss_age_scale_input <- prod.loss_age_scale_all_countries[, country.id]


# Utilities are UK-based and age/gender dependent
qol_inputs <- read.csv(paste0("input/", country.id, "/qol_inputs.csv"), sep=",")

# Utility decrements associated to diabetes-related events
qol_events_inputs <- read.csv(paste0("input/", country.id, "/qol_events_inputs.csv"), sep=",")


# .EVENT variable names

event_vars <- c("CHF.EVENT", "BLIND.EVENT", "ULCER.EVENT", "AMP1.EVENT", "AMP2.EVENT", 
                "MI.EVENT", "IHD.EVENT", "RENAL.EVENT", "STROKE.EVENT")


event_list <- c("CHF", "MI", "IHD", "STROKE", "BLIND", "ULCER", "AMP", "RENAL")

# Indicate the patient characteristics that we will save during the simulation. 
# Besides the risk factors, we have a simulation ID, duration of diabetes 
# and an indicator variable for dead. risk_factors_simulation is defined above.

# Thus, the vector below contains the names of all factors used to predict informal care in the model.
risk_factors_informal   <- c("FEMALE",   "CURR.AGE.SCALE.INF", "CHF.HIST", "STROKE.HIST", "RENAL.HIST")
risk_factors_prod       <- c("FEMALE",   "CURR.AGE.SCALE.PROD", "CHF.HIST", "STROKE.HIST", "RENAL.HIST")
risk_factors_employment <- c("CURR.AGE", "CURR.AGE.2", "FEMALE")

history_characteristics <- c("SIMID", 
                             unique(sort(c(risk_factors_simulation, 
                                           risk_factors_informal, 
                                           risk_factors_prod,
                                           risk_factors_employment))), 
                             "INF.CARE", "EMPLOYED","PROD.LOSS", #Added 29/08/2020
                             "SDURATION",
                             "dead")


# Export tables into Excel using this function:
export_csv <- function(object_input){
  write.csv(object_input, paste(results_dir, substitute(object_input),"_", tx_label,"_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".csv", sep = ""))
}


#############################################################################
########## PART III - INFORMAL CARE AND PRODUCTIVITY EQUATIONS ##############
#############################################################################

# The equation used to predict informal care use in the Netherlands is the following:

# And now, below we have the regression coefficients as reported above.
# CHECK that female coefficient is actually for females

# UPDATE EQUATIONS WITH MOST RECENT VERSION: 10/12/2020
# Note there is no country specific coefficient now. That's why the last coefficient = 0. 
# Otherwise the code would not work.
# Updated 15/11/2021: country-specific equations
informal_care_coef_all_countries <- data.frame(
  UK = c(-1.5355, 0.5036, 0.6852, 0.4231, 0.2245, 0.5498, 0),
  NL = c(-1.5355, 0.5036, 0.6852, 0.4231, 0.2245, 0.5498, 0),
  DE = c(-1.5355, 0.5036, 0.6852, 0.4231, 0.2245, 0.5498, 0),
  ES = c(-1.5259, 0.5992, 0.9210, 0.3259, 1.1048, 0.2057, 0),
  GR = c(-1.5259, 0.5992, 0.9210, 0.3259, 1.1048, 0.2057, 0)
)

informal_care_coef_input <- informal_care_coef_all_countries[, country.id]

employment_coef_all_countries <- data.frame(
  UK = c(2.4283,   0.1273, -0.00252, -0.0737, 0),
  NL = c(2.4283,   0.1273, -0.00252, -0.0737, 0),
  DE = c(2.4283,   0.1273, -0.00252, -0.0737, 0),
  ES = c(-11.0097, 0.1664, 0, 0.2843, 0),
  GR = c(-11.0097, 0.1664, 0, 0.2843, 0)
)

employment_coef_input <- employment_coef_all_countries[, country.id]

prod_costs_coef_all_countries <- data.frame(
  UK = c(-1.4747, 0.01047, 0.5112, 0.08407, 0.1881, 1.3026, 0),
  NL = c(-1.4747, 0.01047, 0.5112, 0.08407, 0.1881, 1.3026, 0),
  DE = c(-1.4747, 0.01047, 0.5112, 0.08407, 0.1881, 1.3026, 0),
  ES = c(-1.5503, 0.1943,  0.4042, 0.1543,  0.3804, 0.4103, 0),
  GR = c(-1.5503, 0.1943,  0.4042, 0.1543,  0.3804, 0.4103, 0)
)

prod_costs_coef_input <- prod_costs_coef_all_countries[, country.id]


# FIXME: objects below are only used to construct the data frames 16 lines below. Are they needed for model functionality or can they be removed?
informal_care_coef_UK_NL_DE <- c(-1.5355, 0.5036, 0.6852, 0.4231, 0.2245, 0.5498, 0) 
informal_care_coef_ES_GR    <- c(-1.5259, 0.5992, 0.9210, 0.3259, 1.1048, 0.2057, 0) 

employment_coef_UK_NL_DE    <- c(2.4283,   0.1273, -0.00252, -0.0737, 0) 
employment_coef_ES_GR       <- c(-11.0097, 0.1664, 0, 0.2843, 0) 

prod_costs_coef_UK_NL_DE    <- c(-1.4747, 0.01047, 0.5112, 0.08407, 0.1881, 1.3026, 0)
prod_costs_coef_ES_GR       <- c(-1.5503, 0.1943,  0.4042, 0.1543,  0.3804, 0.4103, 0)  

# Besides the risk factors, there are two other parameters used in the UKPDS equations called "lambda" and "ro". We added them below.
# Now we have the names of all the parameters (lambda, ro and risk factors) used to predict informal risks in the model.
parameters_informal   <- c("lambda", risk_factors_informal, country.id)
parameters_prod       <- c("lambda", risk_factors_prod, country.id)
parameters_employment <- c("lambda", risk_factors_employment, country.id)

#FIXME: are the objects below used at any point in the model? 
# Below we simply create a table (R data frame) with all the coefficients of the regression equation to predict informal care
informal_care_equations <- data.frame(informal_care_coef_UK_NL_DE, informal_care_coef_ES_GR,row.names = parameters_informal)
prod_cost_equations     <- data.frame(prod_costs_coef_UK_NL_DE, prod_costs_coef_ES_GR, row.names = parameters_prod)
employment_equations    <- data.frame(employment_coef_UK_NL_DE, employment_coef_ES_GR, row.names = parameters_employment)

# Notation: as a general rule, names for function input parameters end with "_input".

# Note that age should be scaled before passing it to the function to calculate the probability of receiving informal care! 
annual_p_bernoulli <- function(regression_coefficents_input, risk_factors_input){
  risk_factors_input <- as.numeric(risk_factors_input) # delete if not needed
  # Note: mind the order of the regression coefficients. 
  log.oods <- sum(regression_coefficents_input*c(1,risk_factors_input,1))
  p <- exp(log.oods)/(1+exp(log.oods))
  return(list(p = p))
}
