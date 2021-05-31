#########################################################
########## VALIDATION UKPDS RISK EQUATIONS ##############
#########################################################

source("setup.R")

# Patient characteristics
validation_patient <- read.csv("input/UK/baseline_characteristics_UK.csv", sep=",")

retirement_age_input <- 65

validation_patient$AGE.DIAG <- validation_patient$CURR.AGE - validation_patient$YEAR 
validation_patient$CURR.AGE.SCALE.INF  <- (validation_patient$CURR.AGE - 72.5474088)/10.4626624 
validation_patient$CURR.AGE.SCALE.PROD <- (validation_patient$CURR.AGE - 60.2737989)/6.1177269 
validation_patient$CURR.AGE.2 <- (validation_patient$CURR.AGE)^2
validation_patient$INF.CARE   <- 0
# ifelse(validation_patient$CURR.AGE >= retirement_age_input, validation_patient$EMPLOYED <- 0,
#        {baseline_employed_prob <- apply(validation_patient %>% select(risk_factors_employment), 1, function(x) annual_p_bernoulli(employment_equations$employment_coef,x)$p)
#        validation_patient$EMPLOYED <- unlist(lapply(baseline_employed_prob, function(x) rbinom(1,1,x))) #EMPLOYED = yes/no
#        })
# validation_patient$PROD.LOSS <- 0 
validation_patient$BMI1 <- if_else(validation_patient$BMI < 18.5, 1, 0)
validation_patient$BMI3 <- if_else(validation_patient$BMI >= 25, 1, 0)


# .EVENT variable names

event_vars <- c("CHF.EVENT", "BLIND.EVENT", "ULCER.EVENT", "AMP1.EVENT", "AMP2.EVENT", 
                "MI.EVENT", "IHD.EVENT", "RENAL.EVENT", "STROKE.EVENT")
validation_patient[event_vars] <- 0 

validation_patient$eGFR       <- validation_patient$eGFR/10 
validation_patient$eGFR60more <- if_else(validation_patient$eGFR > 6,  validation_patient$eGFR,0)
validation_patient$eGFR60less <- if_else(validation_patient$eGFR <= 6, validation_patient$eGFR,0) #defining it like this makes the equation for MI females unrealistically high

validation_patient$HDL        <- validation_patient$HDL*10
validation_patient$HEART.R    <- validation_patient$HEART.R/10 

validation_patient$LDL        <- validation_patient$LDL*10
validation_patient$LDL35more  <- if_else(validation_patient$LDL >= 35, validation_patient$LDL, 0) # not sure either, but I think it does not make sense to cut at 35

validation_patient$SBP        <- validation_patient$SBP/10

validation_patient$CHF.HIST <- 0

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
  
  return(list(H_t1 = H_t1, # no need to return the H's, keep them for now just for validation purposes
              H_t2 = H_t2,
              linear_predictor = linear_predictor,
              p = p))
}



# The annual probability of experiencing diabetic ulcer, death in first year of event and death in subsequent years of events,
# is assumed to follow a logistic distribution where t is the current age.
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


# Validation of the events 

# CHF females and males is the same equation
#male_macro <- unlist(validation_patient[1,] %>% select(risk_factors_macrovascular))
#annual_p_weibull(macrovascular_risk_equations$CHF,male_macro,validation_patient[1,"YEAR"])$p

p_CHF <- annual_p_weibull(macrovascular_risk_equations$CHF,validation_patient[2,] 
                          %>% select(risk_factors_macrovascular),validation_patient[2,"YEAR"])$p
p_CHF
mean(rbinom(1000,1,p_CHF))


# IHD females and males. Is a difference between gender expected? Yes, negative coefficient for females.
p_IHD_female <- annual_p_weibull(macrovascular_risk_equations$IHD,validation_patient[1,] 
                                 %>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
p_IHD_female
mean(rbinom(1000,1,p_IHD_female))

p_IHD_male <- annual_p_weibull(macrovascular_risk_equations$IHD,validation_patient[2,]%>% select(risk_factors_macrovascular),validation_patient[2,"YEAR"])$p
p_IHD_male
mean(rbinom(1000,1,p_IHD_male))


# First MI males
p_MI_male <- annual_p_weibull(macrovascular_risk_equations$FMIMALE,validation_patient[2,] %>% select(risk_factors_macrovascular),validation_patient[2,"YEAR"])$p
p_MI_male
# [1] 0.002285789

# First MI females
p_MI_female <- annual_p_weibull(macrovascular_risk_equations$FMIFEMALE,validation_patient[1,] %>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
p_MI_female
# [1] 0.004844228

# Second MI females and males
p_SMI <- annual_p_weibull(macrovascular_risk_equations$SMI,validation_patient[1,] %>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
#annual_p_weibull(macrovascular_risk_equations$SMI,validation_patient[2,] %>% select(risk_factors_macrovascular),validation_patient[2,"YEAR"])$p
p_SMI
mean(rbinom(1000,1,p_SMI))
# [1] 0.01946934

# First stroke females
p_stroke_female <- annual_p_weibull(macrovascular_risk_equations$FSTROKE,validation_patient[1,]%>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
p_stroke_female

# First stroke males
p_stroke_male <- annual_p_weibull(macrovascular_risk_equations$FSTROKE,validation_patient[2,] %>% select(risk_factors_macrovascular),validation_patient[2,"YEAR"])$p
p_stroke_male


# Second stroke females and males
p_sstroke <- annual_p_weibull(macrovascular_risk_equations$SSTROKE,validation_patient[1,] %>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
p_sstroke


# Blindness females
p_blind <- annual_p_weibull(microvascular_risk_equations$BLIND,validation_patient[1,] %>% select(risk_factors_microvascular),validation_patient[1,"YEAR"])$p
p_blind

# First amputation no ulcer females
p_amp_no_ulcer_female <- annual_p_weibull(microvascular_risk_equations$FAMPNOULCER,validation_patient[1,] %>% select(risk_factors_microvascular),validation_patient[1,"YEAR"])$p
p_amp_no_ulcer_female
# [1] 0.0002014471

# First amputation no ulcer males
p_amp_no_ulcer_male <- annual_p_weibull(microvascular_risk_equations$FAMPNOULCER,validation_patient[2,] %>% select(risk_factors_microvascular),validation_patient[2,"YEAR"])$p
p_amp_no_ulcer_male

# First amputation ulcer females and males
p_amp_ulcer <- annual_p_weibull(microvascular_risk_equations$FAMPULCER,validation_patient[1,] %>% select(risk_factors_microvascular),validation_patient[1,"YEAR"])$p
p_amp_ulcer

# Second amputation 
p_sec_amp <- annual_p_weibull(microvascular_risk_equations$SAMP,validation_patient[1,] %>% select(risk_factors_microvascular),validation_patient[1,"YEAR"])$p
p_sec_amp

# Renal failure females
p_renal_female <- annual_p_weibull(microvascular_risk_equations$RENALF,validation_patient[1,] %>% select(risk_factors_microvascular),validation_patient[1,"YEAR"])$p
p_renal_female

# Renal failure males
p_renal_male <- annual_p_weibull(microvascular_risk_equations$RENALF,validation_patient[2,] %>% select(risk_factors_microvascular),validation_patient[2,"YEAR"])$p
p_renal_male


# Ulcer females
p_ulcer_female <- annual_p_logistic(microvascular_risk_equations$ULCER,validation_patient[1,] %>% select(risk_factors_microvascular))$p
p_ulcer_female

# Ulcer males
p_ulcer_male <- annual_p_logistic(microvascular_risk_equations$ULCER,validation_patient[2,] %>% select(risk_factors_microvascular))$p
p_ulcer_male


# Death in 1st year of events females
annual_p_logistic(mortality_risk_equations$DEATH1YEVENT,validation_patient[1,] %>% select(risk_factors_mortality))$p
annual_p_logistic(mortality_risk_equations$DEATH1YEVENT,validation_patient[2,] %>% select(risk_factors_mortality))$p
# [1] 0.09313699

# Death in subsequent year of events females
annual_p_logistic(mortality_risk_equations$DEATHYSEVENT,validation_patient[1,] %>% select(risk_factors_mortality))$p
annual_p_logistic(mortality_risk_equations$DEATHYSEVENT,validation_patient[2,] %>% select(risk_factors_mortality))$p
# [1] 0.6601135


# Death in year with no history or events females
annual_p_gompertz(mortality_risk_equations$DEATHNOHIST, validation_patient[1,] %>% select(risk_factors_mortality),validation_patient[1,]$AGE.DIAG + validation_patient[1,]$YEAR)$p
annual_p_gompertz(mortality_risk_equations$DEATHNOHIST, validation_patient[1,] %>% select(risk_factors_mortality),validation_patient[1,]$CURR.AGE)$p # not sure why in the simulation code we used the line above, it might give problems when updating and it's slower? updating seems fine because YEAR is updating
# [1] 0.004071926

# Death in year with no history or events males
annual_p_gompertz(mortality_risk_equations$DEATHNOHIST, validation_patient[2,] %>% select(risk_factors_mortality),validation_patient[1,]$AGE.DIAG + validation_patient[1,]$YEAR)$p
annual_p_gompertz(mortality_risk_equations$DEATHNOHIST, validation_patient[2,] %>% select(risk_factors_mortality),validation_patient[1,]$CURR.AGE)$p # not sure why in the simulation code we used the line above, it might give problems when updating and it's slower? updating seems fine because YEAR is updating
# [1] 0.005117118

# Death in year with history but no events females
annual_p_gompertz(mortality_risk_equations$DEATHHISTNOEVENT, validation_patient[1,] %>% select(risk_factors_mortality),validation_patient[1,]$AGE.DIAG + validation_patient[1,]$YEAR)$p
annual_p_gompertz(mortality_risk_equations$DEATHHISTNOEVENT, validation_patient[1,] %>% select(risk_factors_mortality),validation_patient[1,]$CURR.AGE)$p # not sure why in the simulation code we used the line above, it might give problems when updating and it's slower? updating seems fine because YEAR is updating
# [1] 0.006813194

# Death in year with history but no events males
annual_p_gompertz(mortality_risk_equations$DEATHHISTNOEVENT, validation_patient[2,] %>% select(risk_factors_mortality),validation_patient[1,]$AGE.DIAG + validation_patient[1,]$YEAR)$p
annual_p_gompertz(mortality_risk_equations$DEATHHISTNOEVENT, validation_patient[2,] %>% select(risk_factors_mortality),validation_patient[1,]$CURR.AGE)$p # not sure why in the simulation code we used the line above, it might give problems when updating and it's slower? updating seems fine because YEAR is updating
# [1] 0.006813194




### Include order as in UKPDS model

# Validation of the events 

validation_patient <- validation_patient[1,]

CHF_event_function <- function(){
  
  if(validation_patient$CHF.HIST == 0){
    current_CHF_prob  <- 0.5#annual_p_weibull(macrovascular_risk_equations$CHF,validation_patient%>% select(risk_factors_macrovascular),validation_patient[,"YEAR"])$p
    validation_patient$CHF.EVENT <<- rbinom(1,1,current_CHF_prob)
    if(validation_patient$CHF.EVENT == 1){validation_patient$CHF.HIST <<- 1}
  } 
  
  print(current_CHF_prob)
  
  }


IHD_event_function <- function(){
  
  if(validation_patient$IHD.HIST == 0){
    current_IHD_prob  <- annual_p_weibull(macrovascular_risk_equations$IHD,validation_patient%>% select(risk_factors_macrovascular),validation_patient[,"YEAR"])$p
    validation_patient$IHD.EVENT <<- rbinom(1,1,current_IHD_prob) 
    if(validation_patient$IHD.EVENT == 1){validation_patient$IHD.HIST <<- 1}
  } 
  
  print(current_IHD_prob)
}




# First MI males
p_MI_male <- annual_p_weibull(macrovascular_risk_equations$FMIMALE,validation_patient[2,] %>% select(risk_factors_macrovascular),validation_patient[2,"YEAR"])$p
p_MI_male
# [1] 0.002285789

# First MI females
p_MI_female <- annual_p_weibull(macrovascular_risk_equations$FMIFEMALE,validation_patient[1,] %>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
p_MI_female
# [1] 0.004844228

# Second MI females and males
p_SMI <- annual_p_weibull(macrovascular_risk_equations$SMI,validation_patient[1,] %>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
#annual_p_weibull(macrovascular_risk_equations$SMI,validation_patient[2,] %>% select(risk_factors_macrovascular),validation_patient[2,"YEAR"])$p
p_SMI
mean(rbinom(1000,1,p_SMI))
# [1] 0.01946934

# First stroke females
p_stroke_female <- annual_p_weibull(macrovascular_risk_equations$FSTROKE,validation_patient[1,]%>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
p_stroke_female

# First stroke males
p_stroke_male <- annual_p_weibull(macrovascular_risk_equations$FSTROKE,validation_patient[2,] %>% select(risk_factors_macrovascular),validation_patient[2,"YEAR"])$p
p_stroke_male


# Second stroke females and males
p_sstroke <- annual_p_weibull(macrovascular_risk_equations$SSTROKE,validation_patient[1,] %>% select(risk_factors_macrovascular),validation_patient[1,"YEAR"])$p
p_sstroke


# Blindness females
p_blind <- annual_p_weibull(microvascular_risk_equations$BLIND,validation_patient[1,] %>% select(risk_factors_microvascular),validation_patient[1,"YEAR"])$p
p_blind

# First amputation no ulcer females
p_amp_no_ulcer_female <- annual_p_weibull(microvascular_risk_equations$FAMPNOULCER,validation_patient[1,] %>% select(risk_factors_microvascular),validation_patient[1,"YEAR"])$p
p_amp_no_ulcer_female
# [1] 0.0002014471

# First amputation no ulcer males
p_amp_no_ulcer_male <- annual_p_weibull(microvascular_risk_equations$FAMPNOULCER,validation_patient[2,] %>% select(risk_factors_microvascular),validation_patient[2,"YEAR"])$p
p_amp_no_ulcer_male

# First amputation ulcer females and males
p_amp_ulcer <- annual_p_weibull(microvascular_risk_equations$FAMPULCER,validation_patient[1,] %>% select(risk_factors_microvascular),validation_patient[1,"YEAR"])$p
p_amp_ulcer

# Second amputation 
p_sec_amp <- annual_p_weibull(microvascular_risk_equations$SAMP,validation_patient[1,] %>% select(risk_factors_microvascular),validation_patient[1,"YEAR"])$p
p_sec_amp

# Renal failure females
p_renal_female <- annual_p_weibull(microvascular_risk_equations$RENALF,validation_patient[1,] %>% select(risk_factors_microvascular),validation_patient[1,"YEAR"])$p
p_renal_female

# Renal failure males
p_renal_male <- annual_p_weibull(microvascular_risk_equations$RENALF,validation_patient[2,] %>% select(risk_factors_microvascular),validation_patient[2,"YEAR"])$p
p_renal_male


# Ulcer females
p_ulcer_female <- annual_p_logistic(microvascular_risk_equations$ULCER,validation_patient[1,] %>% select(risk_factors_microvascular))$p
p_ulcer_female

# Ulcer males
p_ulcer_male <- annual_p_logistic(microvascular_risk_equations$ULCER,validation_patient[2,] %>% select(risk_factors_microvascular))$p
p_ulcer_male


# Death in 1st year of events females
annual_p_logistic(mortality_risk_equations$DEATH1YEVENT,validation_patient[1,] %>% select(risk_factors_mortality))$p
annual_p_logistic(mortality_risk_equations$DEATH1YEVENT,validation_patient[2,] %>% select(risk_factors_mortality))$p
# [1] 0.09313699

# Death in subsequent year of events females
annual_p_logistic(mortality_risk_equations$DEATHYSEVENT,validation_patient[1,] %>% select(risk_factors_mortality))$p
annual_p_logistic(mortality_risk_equations$DEATHYSEVENT,validation_patient[2,] %>% select(risk_factors_mortality))$p
# [1] 0.6601135


# Death in year with no history or events females
annual_p_gompertz(mortality_risk_equations$DEATHNOHIST, validation_patient[1,] %>% select(risk_factors_mortality),validation_patient[1,]$AGE.DIAG + validation_patient[1,]$YEAR)$p
annual_p_gompertz(mortality_risk_equations$DEATHNOHIST, validation_patient[1,] %>% select(risk_factors_mortality),validation_patient[1,]$CURR.AGE)$p # not sure why in the simulation code we used the line above, it might give problems when updating and it's slower? updating seems fine because YEAR is updating
# [1] 0.004071926

# Death in year with no history or events males
annual_p_gompertz(mortality_risk_equations$DEATHNOHIST, validation_patient[2,] %>% select(risk_factors_mortality),validation_patient[1,]$AGE.DIAG + validation_patient[1,]$YEAR)$p
annual_p_gompertz(mortality_risk_equations$DEATHNOHIST, validation_patient[2,] %>% select(risk_factors_mortality),validation_patient[1,]$CURR.AGE)$p # not sure why in the simulation code we used the line above, it might give problems when updating and it's slower? updating seems fine because YEAR is updating
# [1] 0.005117118

# Death in year with history but no events females
annual_p_gompertz(mortality_risk_equations$DEATHHISTNOEVENT, validation_patient[1,] %>% select(risk_factors_mortality),validation_patient[1,]$AGE.DIAG + validation_patient[1,]$YEAR)$p
annual_p_gompertz(mortality_risk_equations$DEATHHISTNOEVENT, validation_patient[1,] %>% select(risk_factors_mortality),validation_patient[1,]$CURR.AGE)$p # not sure why in the simulation code we used the line above, it might give problems when updating and it's slower? updating seems fine because YEAR is updating
# [1] 0.006813194

# Death in year with history but no events males
annual_p_gompertz(mortality_risk_equations$DEATHHISTNOEVENT, validation_patient[2,] %>% select(risk_factors_mortality),validation_patient[1,]$AGE.DIAG + validation_patient[1,]$YEAR)$p
annual_p_gompertz(mortality_risk_equations$DEATHHISTNOEVENT, validation_patient[2,] %>% select(risk_factors_mortality),validation_patient[1,]$CURR.AGE)$p # not sure why in the simulation code we used the line above, it might give problems when updating and it's slower? updating seems fine because YEAR is updating
# [1] 0.006813194

