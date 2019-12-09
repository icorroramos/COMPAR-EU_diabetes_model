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


#######################################################
########## PART I - UKPDS RISK EQUATIONS ##############
#######################################################

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

risk_factors_macrovascular <- c("AFRO", "AGE.DIAG", "FEMALE", "INDIAN", "ATFIB", "BMI", "eGFR", "eGFR60less", "HbA1c", "HDL",
                                "LDL", "LDL35more", "MMALB", "PVD", "SBP", "SMOKER", "WBC", "AMP.HIST", "CHF.HIST", "IHD.HIST", 
                                "STROKE.HIST", "ULCER.HIST")

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

risk_factors_microvascular <- c("AFRO", "AGE.DIAG", "FEMALE", "ATFIB", "BMI", "eGFR60less", "eGFR60more", "HAEM", "HbA1c", "HDL",
                                "HEART.R", "LDL", "MMALB", "PVD", "SBP", "WBC", "AMP.HIST", "BLIND.HIST", "CHF.HIST", 
                                "IHD.HIST", "STROKE.HIST")

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

risk_factors_mortality <- c("FEMALE", "INDIAN", "YEAR", "ATFIB", "BMI1", "BMI3", "CURR.AGE", "HDL", "HEART.R", 
                            "MMALB", "PVD", "SMOKER", "WBC", "AMP.EVENT", "AMP.HIST", "AMP2.EVENT", "CHF.HIST", 
                            "IHD.EVENT", "IHD.HIST", "MI.EVENT", "MI.HIST", "RENAL.EVENT", "RENAL.HIST", "STROKE.EVENT", "STROKE.HIST")

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

cum_haz_weibull(mortality_risk_equations$DEATHNOHIST,c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),8)


### Below all the patient characteristics (38) defining a patient. Not sure if we want a patient id (PTID)
risk_factors_simulation <- unique(sort(c(risk_factors_macrovascular, risk_factors_microvascular, risk_factors_mortality)))

# > risk_factors_simulation
# [1] "AFRO"         "AGE.DIAG"     "AMP.EVENT"    "AMP.HIST"     "AMP2.EVENT"   "ATFIB"        "BLIND.HIST"   "BMI"         
# [9] "BMI1"         "BMI3"         "CHF.HIST"     "CURR.AGE"     "eGFR"         "eGFR60less"   "eGFR60more"   "FEMALE"      
# [17] "HAEM"         "HbA1c"        "HDL"          "HEART.R"      "IHD.EVENT"    "IHD.HIST"     "INDIAN"       "LDL"         
# [25] "LDL35more"    "MI.EVENT"     "MI.HIST"      "MMALB"        "PVD"          "RENAL.EVENT"  "RENAL.HIST"   "SBP"         
# [33] "SMOKER"       "STROKE.EVENT" "STROKE.HIST"  "ULCER.HIST"   "WBC"          "YEAR"        

# Describe how are these variables defined: categorial or numerical

# Besides the risk factors, we need to add time, which is duration of diabetes (DDURATION) and dead status

##############################################
########## PART II - SIMULATION ##############
##############################################


# The code below is adapted from the COPD model. --> Work in progress...

# The main simulation starts below. The code is used to 1) simulate patients’ clinical history, 2) calculate costs and 
# 3) calculate QALYs. 

# Focus for now in clinical history only. 
# PSA to be calculated later.

# The input parameters of the SMDMII_model_simulation function are the following:
# 1. patient_size_input = number of patients included in the simulation.
# 2. run_PSA_input = runs the model in probabilistic mode. Otherwise, deterministic.

# Treatment effect parameters to be defined. For the moment just a dummy:
# 3. treatment_effect_input = variable to indicate treatment effect. Default should be 1.

# A random seed that it is used to ensure consistency in the model results. Not sure if needed but keep it for the moment. 

# 4. seed_input = 

SMDMII_model_simulation <- function(patient_size_input,
                                    run_PSA_input,
                                    treatment_effect_input,
                                    seed_input){
  
  ### Step 1: Read the data with the patient characteristics: not sure how we will do this yet
  
  ### Complete cases in the COPD model was a subset of the available data: patients with all predictors observed
  ### Keep it for now but consider changing it later.
  
  all_baseline_patients <- matrix(c(rep(0,length(risk_factors_simulation)),8,0),ncol = 2+ length(risk_factors_simulation)) # dummy value
  
  ### Step 2: indicate the patient characteristics that we will save during the simulation. 
  ###         Besides the risk factors, we have a simulation ID, may have a patient ID,
  ###         duration of diabetes (not necesarily 0 at baseline?) and indicator variable for dead
  
  #history_characteristics <- c("SIMID","PTID",risk_factors_simulation, "DDURATION","dead")
  history_characteristics <- c("SIMID", risk_factors_simulation, "DDURATION","dead")
  
  
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
  
  ### Select the sample used in the simulation: sample with replacement from the dataset
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
    
    # Pick the current patient from those selected fomr baseline
    current_patient <- simulation_baseline_patients[patient_index,]
    current_patient$SIMID <- patient_index
    
    # Save the characteristics to be used in the simulation history (not those that are stable, only those changing) 
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
      
      # NOT SURE WHETHER THERE IS ANY ORDER ON THIS: ALL EQUATIONS AT THE SAME TIME OR IN SEQUENTIAL ORDER? -->>
      # -->> DEATH IS CALCULATED AFTERWARDS
      # IT IS CLEAR THAT MORE THAN ONE EVENT CAN OCCUR IN ONE YEAR. WHAT IT IS NOT SO CLEAR IS WHETHER THESE 
      # EVENTS CAN BE ANY OF THEM OR THOSE DEPENDING OR NOT ON THE HISTORY
      # FOR EXAMPLE, IF ULCER HAPPENS IN ONE YEAR, THE PROBABILITY OF AMPUTATION DEPENDS ON THE ULCER HISTORY
      # THE QUESTION IS: WILL THE HISTORY BE UPDATED FOR THE CURRENT YEAR OR THE NEXT YEAR? IN THE CURRENT YEAR,
      # THE PROBABILITY OF AMPUTATION DEPENDS ON THE YEAR BEFORE OR THE CURRENT YEAR? THE LATTER WOULD IMPLY 
      # SOME ORDER IN THE OCCURRENCE OF EVENTS. HERE THERE IS NO TIME TO EVENT, SO THE ORDER CANNOT BE ESTABLISHED
      # HOWEVER, IN THE DESCRIPTION OF THE MODEL (UKPDS PAPER) IT SAYS "RANDOMLY ORDER AND RUN EVENT EQUATIONS" -->>
      # ORDER WOULD IMPLY SEQUENTIAL RUN OF EQUATIONS AND UPDATES COULD BE BASED ON CURRENT YEAR 
      
      
      ### MACROVASCULAR COMPLICATIONS
      
      # Heart Failure is Weibull (cum_haz_weibull)
      
      #current_CHF_parameters <- predicted_exacerbation_weibull(exacerbation_weibull_regression_coef,current_patient[exacerbation_predictors])
      #current_CHF_prob       <- rweibull(1,current_CHF_parameters$shape_exacerbation_weibull,current_CHF_parameters$scale_exacerbation_weibull)
      #current_CHF_prob       <- treatment_effect_input*current_CHF_prob
      
      # IHD is Weibull
      #current_IHD_parameters <- predicted_pneumonia_weibull(pneumonia_weibull_regression_coef,current_patient[pneumonia_predictors])
      #current_IHD_prob       <- rweibull(1,current_IHD_parameters$shape_pneumonia_weibull,current_IHD_parameters$scale_pneumonia_weibull)
      
      
      # MI could be first or second. If it is first then it is different for male and female. If it is second, then it is the same for both genders.
      
      if(current_patient$MI.HIST == 0){
        
        # If no history of MI, then it is first. Different for males and females
        
        if(current_patient$FEMALE == 1){
          # First MI for female Weibull
          #current_FMIFEMALE_parameters <-
          #current_FMIFEMALE_prob <-
          
        }else{
          # First MI for male Exponential
          #current_FMIMALE_parameters <-
          #current_FMIMALE_prob <-
          
        }
        
      }else{
        # Second MI is Exponential, regardless the gender 
        #current_SMI_parameters <-
        #current_SMI_prob <-
      }
      
       
      # STROKE could be first or second. 
      
      if(current_patient$STROKE.HIST == 0){
        
        # If no history of STROKE, then it is first and Weibull. 
        #current_FSTROKE_parameters <-
        #current_FSTROKE_prob <-
        
      }else{
        # Second STROKE is Weibull
        #current_SSTROKE_parameters <-
        #current_SSTROKE_prob <-
      }
      
      
      ### MICROVASCULAR COMPLICATIONS
      
      # BLINDNESS is Exponential 

      #current_BLIND_parameters
      #current_BLIND_prob      
  
      # ULCER is Logistic 
      
      #current_ULCER_parameters
      #current_ULCER_prob      

      #CONTINUE FROM HERE
      
      FAMPNOULCER <- 
      
      FAMPULCER <- 
      
      SAMP <- 
      
      RENALF <- 
      
      
      
      
      
      ##############################################
      # STEP 3.2: Update patient characteristics   #
      ##############################################
      
      # We first copy all the previous characteristics
      current_patient_update <- current_patient
      
      #######################################################
      # STEP 3.2.1: Update first the lagged characteristics #
      #######################################################
      
      current_patient_update$lag_FEVPPA        <- current_patient$FEVPPA
      
      current_patient_update$lag_SGACT         <- current_patient$SGACT
      
      
      current_patient_update$lag_CWE_TOT       <- current_patient$CWE_TOT
      current_patient_update$lag_BREATHLESS_yn  <- current_patient$BREATHLESS_yn
      current_patient_update$lag_COUGHSPUTUM_yn <- current_patient$COUGHSPUTUM_yn
      
      current_patient_update$lag_SGTOT         <- current_patient$SGTOT
      
      
      current_patient_update$PREV_SEVEXAC_yn <- current_patient$SEVEXAC_yn
      current_patient_update$PREV_TOTEXAC_yn <- ifelse(current_patient$SEVEXAC_yn==1 | current_patient$MODEXAC_yn==1,1,0)
      
      
      #######################################################
      # STEP 3.2.2: Update depending on the event occurred  #
      #######################################################
      
      # If exacerbation happened, then update at exacerbation time
      if(min(current_remaining_life_exp,current_CHF_prob,current_IHD_prob)==current_CHF_prob){
        
        # Patient is still alive
        current_patient_update$dead <- 0
        
        # Update time-relatedcharacteristics
        current_patient_update$ANLYEAR         <- current_patient$ANLYEAR  + current_CHF_prob
        current_patient_update$ANLYEAR_SCALED  <- (current_patient_update$ANLYEAR - 1.529411)/1.376549 # hardcoded based on dataset
        current_patient_update$lag_ANLYEAR     <- current_CHF_prob
        current_patient_update$AGE_TIME        <- current_patient$AGE_TIME + current_CHF_prob
        
        
        # If age + current life expectancy is > 100 years then we force death
        if(current_patient_update$AGE_TIME>100){current_patient_update$dead <- 1}
        
        # Update exacerbation status. Decide first whether the exacerbation is moderate or severe.
        
        current_severity_prob   <- predicted_exacerbation_severity(exacerbation_severity_regression_coef,current_patient_update[exacerbation_severity_predictors])$p.exacerbation_severity
        current_severity_sample <- rbinom(1,1,current_severity_prob)
        if(current_severity_sample==1){
          current_patient_update$MODEXAC_yn <- 0
          current_patient_update$SEVEXAC_yn <- 1
          
          # Additional death risk because of severe exacerbation 
          current_sevexa_death_prob   <- 0.063 ### hardcoded for now!
          current_sevexa_death_sample <- rbinom(1,1,current_sevexa_death_prob)
          
          if(current_sevexa_death_sample==1){
            current_patient_update$dead <- 1
          }else{
            current_patient_update$dead <- 0
          }
          
          
        }else{
          current_patient_update$MODEXAC_yn <- 1
          current_patient_update$SEVEXAC_yn <- 0
        }
        
        # Update pneumonia status 
        current_patient_update$PNEU_yn      <- 0
        current_patient_update$pneu_hosp_yn <- 0
        
      }
      
      
      # If pneumonia happened, then update at pneumonia time
      if(min(current_remaining_life_exp,current_CHF_prob,current_IHD_prob)==current_IHD_prob){
        
        # Update time-relatedcharacteristics
        current_patient_update$ANLYEAR         <- current_patient$ANLYEAR  + current_IHD_prob
        current_patient_update$ANLYEAR_SCALED  <- (current_patient_update$ANLYEAR - 1.529411)/1.376549 # hardcoded based on dataset
        current_patient_update$lag_ANLYEAR     <- current_IHD_prob
        current_patient_update$AGE_TIME        <- current_patient$AGE_TIME + current_IHD_prob
        
        
        
        # Force death at 100 years
        if(current_patient_update$AGE_TIME>100){current_patient_update$dead <- 1}
        
        # Update pneumonia status. 
        current_patient_update$PNEU_yn      <- 1
        current_patient_update$pneu_hosp_yn <- rbinom(1,1,predicted_pneumonia_hosp(pneumonia_hosp_regression_coef,current_patient_update[pneumonia_hosp_predictors])$p.pneumonia.hosp)
        
        # Patient might die because of pneumonia after hospitalisation
        if(current_patient_update$pneu_hosp_yn==1){
          
          current_pneu_death_prob   <- 1607/19786 ### hardcoded based on dataset
          current_pneu_death_sample <- rbinom(1,1,current_pneu_death_prob)
          
          if(current_pneu_death_sample==1){
            current_patient_update$dead <- 1
          }else{
            current_patient_update$dead <- 0
          }
        }
        
      }
      
      
      # If death  happened then update age and finish the simulation for this patient
      if(min(current_remaining_life_exp,current_CHF_prob)==current_remaining_life_exp){
        
        # Patient is dead
        current_patient_update$dead <- 1
        
        # Update time-relatedcharacteristics
        current_patient_update$ANLYEAR          <- current_patient$ANLYEAR  + current_remaining_life_exp
        current_patient_update$ANLYEAR_SCALED   <- (current_patient_update$ANLYEAR - 1.529411)/1.376549 # hardcoded based on dataset
        current_patient_update$lag_ANLYEAR      <- current_remaining_life_exp
        current_patient_update$AGE_TIME         <- current_patient$AGE_TIME + current_remaining_life_exp
        
        
        # Update exacerbation status. 
        current_patient_update$MODEXAC_yn <- 0
        current_patient_update$SEVEXAC_yn <- 0
        
        # Update pneumonia status. 
        current_patient_update$PNEU_yn      <- 0
        current_patient_update$pneu_hosp_yn <- 0
      }
      
      
      ############################################################################
      # STEP 3.2.3: Update continuous variables depending on the event occurred  #
      ############################################################################
      
      # Update FEV1, CWE, SGACT, symptoms and SGTOT. The order is important here:
      
      # Update FEV1
      current_patient_update$FEVA          <- max(0,predicted_fev1(lung_function_regression_coef,current_patient_update[fev1_predictors],fev1_treatment_effect_input)$fev_1) 
      current_FEVPPA_calc                  <- FEVPPA_calc(current_patient_update$FEMALE,current_patient_update$HTSTD,current_patient_update$AGE_TIME,current_patient_update$FEVA)
      current_patient_update$FEVPPA        <- current_FEVPPA_calc$FEVPPA
      current_patient_update$FEV1pred      <- current_FEVPPA_calc$FEV1_pred
      
      
      # Force death if FEV1 < 0.2
      if(current_patient_update$FEVA < 0.2){current_patient_update$dead <- 1}
      
      
      # Update CWE (min observed is 42. here for now we truncate at 0)
      current_patient_update$CWE_TOT        <- max(0,predicted_cwe_tot(cwe_tot_regression_coef,current_patient_update[cwe_tot_predictors])$cwe_tot) #max(0,predicted_cwe_tot(cwe_tot_regression_coef$Value,current_patient_update[cwe_tot_predictors])$cwe_tot)
      
      
      # Update SGACT
      current_patient_update$SGACT        <- min(100,max(0,predicted_SGACT(SGACT_regression_coef,current_patient_update[SGACT_predictors])$SGACT))
      
      
      # Update breathlesness and coughsputum
      current_patient_update$BREATHLESS_yn  <- rbinom(1,1,predicted_breathless(breathless_regression_coef,current_patient_update[breathless_predictors])$p.breathless)
      current_patient_update$COUGHSPUTUM_yn <- rbinom(1,1,predicted_coughsputum(coughsputum_regression_coef,current_patient_update[coughsputum_predictors])$p.coughsputum)
      
      # Update SGTOT
      current_patient_update$SGTOT        <- min(100,max(0,predicted_SGTOT(SGTOT_regression_coef,current_patient_update[SGTOT_predictors])$SGTOT))
      
      
      ### When all characteristics have been updated, we add these to the patient history
      simulation_patients_history <- rbind(simulation_patients_history,current_patient_update[history_characteristics])
      
      ### And update current patient and go up to while loop
      current_patient <- current_patient_update
      
      
      ######################################################################
      # STEP 3.2.4: Adjust Mortality according to updated characteristics  #
      ######################################################################
      
      ### Fixed to Weibull
      current_mortality_parameters <- predicted_mortality_weibull(mortality_weibull_regression_coef,current_patient[mortality_predictors])
      current_mortality            <- current_mortality_parameters$scale_mortality_weibull*gamma(1+1/current_mortality_parameters$shape_mortality_weibull)/365
      
      ### Adjust remaining life expectancy according to improvement or worsened in condition wrt baseline or previous period
      current_remaining_life_exp <- max(0,(current_remaining_life_exp - current_patient$lag_ANLYEAR)*(current_mortality/lag_current_mortality))
      
      ### Update mortality factor and event index
      lag_current_mortality    <- current_mortality 
      current_sim_time            <- current_sim_time + 1
      
    } #end while loop and move to another patient
    
    patient_index <- patient_index + 1
    
  } #end for loop in number of patients
  
  

  #################################################################
  ########## MAIN PART III: Calculate aggregated results ##########
  #################################################################
  
  ### We first need to make additional columns for the "diff" variables
  patient_event_history_update$diff_ANLYEAR <- "NA"
  patient_event_history_update$diff_FEVA    <- "NA"
  patient_event_history_update$diff_SGTOT   <- "NA"
  patient_event_history_update$diff_SGACT   <- "NA"
  patient_event_history_update$diff_CWE_TOT <- "NA"
  
  ### Calculate the "diff" variables to calculate the change in outcomes per year
  diff_ANLYEAR <- ddply(patient_event_history_update, "SIMID", summarize, diff_ANLYEAR = c(0,diff(ANLYEAR)))
  diff_FEVA    <- ddply(patient_event_history_update, "SIMID", summarize, diff_FEVA    = c(0,diff(FEVA)))
  diff_SGTOT   <- ddply(patient_event_history_update, "SIMID", summarize, diff_SGTOT   = c(0,diff(SGTOT)))
  diff_SGACT   <- ddply(patient_event_history_update, "SIMID", summarize, diff_SGACT   = c(0,diff(SGACT)))
  diff_CWE_TOT <- ddply(patient_event_history_update, "SIMID", summarize, diff_CWE_TOT = c(0,diff(CWE_TOT)))
  
  patient_event_history_update$diff_ANLYEAR <- diff_ANLYEAR$diff_ANLYEAR
  patient_event_history_update$diff_FEVA    <- diff_FEVA$diff_FEVA
  patient_event_history_update$diff_SGTOT   <- diff_SGTOT$diff_SGTOT
  patient_event_history_update$diff_SGACT   <- diff_SGACT$diff_SGACT
  patient_event_history_update$diff_CWE_TOT <- diff_CWE_TOT$diff_CWE_TOT
  
  ### Calculate model outcomes ###
  
  ### Adjust the number of exacerbations (no exacerbation when pneumonia happened)
  patient_event_history_update[which(patient_event_history_update$PNEU_yn==1),"MODEXAC_yn"] <- 0
  patient_event_history_update[which(patient_event_history_update$PNEU_yn==1),"SEVEXAC_yn"] <- 0
  
  
  ### Lung function: Mean FEV1 decline per year
  mean_annual_fev1_decline <- round(mean(patient_event_history_update$diff_FEVA)/mean(patient_event_history_update$diff_ANLYEAR),4)
  
  ### Life expectancy
  mean_life_expectancy <- round(mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  
  # Moderate
  modexac_year <-function(low_limit,upper_limit){
    sum(patient_event_history_update[which(patient_event_history_update$ANLYEAR>low_limit & patient_event_history_update$ANLYEAR<=upper_limit),"MODEXAC_yn"])
  }
  mean_mod_exa_rate <- round(modexac_year(0,1000)/(patient_size_input*mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"])),4)
  
  # Severe 
  sevexac_year <-function(low_limit,upper_limit){
    sum(patient_event_history_update[which(patient_event_history_update$ANLYEAR>low_limit & patient_event_history_update$ANLYEAR<=upper_limit),"SEVEXAC_yn"])
  }
  mean_sev_exa_rate <- round(sevexac_year(0,1000)/(patient_size_input*mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"])),4)
  
  
  ### Exercise capacity
  mean_CWE_TOT_change <- round(mean(patient_event_history_update$diff_CWE_TOT)/mean(patient_event_history_update$diff_ANLYEAR),4)
  
  ### SGRQ activity score 
  mean_SGACT_change <- round(mean(patient_event_history_update$diff_SGACT)/mean(patient_event_history_update$diff_ANLYEAR),4)
  
  ### Cough/sputum
  # Total number of cough/sputum per patient during lifetime
  cum_coughsputum <- aggregate(patient_event_history_update$COUGHSPUTUM_yn, list(Patient = patient_event_history_update$SIMID), sum)
  # Rate per year
  mean_cough_rate <- round(mean(cum_coughsputum$x)/mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  
  # July 2018 added
  # create new time variable
  diff_ANLYEAR_sym <- ddply(patient_event_history_update, "SIMID", summarize, diff_ANLYEAR_sym = c(diff(ANLYEAR),0))
  
  # create data frame just for symptoms
  simulation_clinical_history_sym <- cbind(patient_event_history_update[,c("ANLYEAR", "COUGHSPUTUM_yn", "BREATHLESS_yn","dead")], diff_ANLYEAR_sym)
  simulation_dead                 <- simulation_clinical_history_sym[which(simulation_clinical_history_sym$dead == 1),c("SIMID", "ANLYEAR")]
  simulation_cough                <- simulation_clinical_history_sym[,c("SIMID", "diff_ANLYEAR_sym")]
  simulation_cough_time           <- simulation_clinical_history_sym$COUGHSPUTUM_yn * simulation_clinical_history_sym$diff_ANLYEAR_sym
  simulation_cough                <- cbind(simulation_cough, simulation_cough_time)
  mean_time_cough                 <- mean(aggregate(simulation_cough$simulation_cough_time, list(Patient = simulation_cough$SIMID), sum)$x/simulation_dead$ANLYEAR)
  #ci_time_cough         <- quantile(aggregate(simulation_cough$simulation_cough_time, list(Patient = simulation_cough$SIMID), sum)$x/simulation_dead$ANLYEAR, c(0.025,0.975))
  # End - July 2018 added
  
  ### Shortness of breath
  # Total number of cough/sputum per patient during lifetime
  cum_breathless <- aggregate(patient_event_history_update$BREATHLESS_yn, list(Patient = patient_event_history_update$SIMID), sum)
  # Rate per year
  mean_breathless_rate <- round(mean(cum_breathless$x)/mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  
  
  # July 2018 added
  simulation_breath      <- simulation_clinical_history_sym[,c("SIMID", "diff_ANLYEAR_sym")]
  simulation_breath_time <- simulation_clinical_history_sym$BREATHLESS_yn * simulation_clinical_history_sym$diff_ANLYEAR_sym
  simulation_breath      <- cbind(simulation_breath, simulation_breath_time)
  mean_time_breath       <- mean(aggregate(simulation_breath$simulation_breath_time, list(Patient = simulation_breath$SIMID), sum)$x/simulation_dead$ANLYEAR)
  # End - July 2018 added
  
  
  
  ### Adverse events (pneumonia)
  # Total number of pneumonias per patient during lifetime
  cum_pneu <- aggregate(patient_event_history_update$PNEU_yn, list(Patient = patient_event_history_update$SIMID), sum)
  # Rate per year
  mean_pneu_rate <- round(mean(cum_pneu$x)/mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  
  # Total number of pneumonias leading to hospitalisation per patient during lifetime
  cum_pneu_hosp <- aggregate(patient_event_history_update$pneu_hosp_yn, list(Patient = patient_event_history_update$SIMID), sum)
  # Rate per year
  mean_pneu_hosp_rate <- round(mean(cum_pneu_hosp$x)/mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  

  ### Return model outcomes 
  
  return(list(mean_annual_fev1_decline       = mean_annual_fev1_decline,
              mean_life_expectancy           = mean_life_expectancy,
              mean_qalys                     = mean_qalys,
              mean_qalys_disc                = mean_qalys_disc,
              mean_total_costs_hc            = mean_total_costs_hc,
              mean_total_costs_hc_disc       = mean_total_costs_hc_disc))
} #end SMDMII_model_simulation function

# To run the model simply call the model function with the appropriate inputs in the correct order. For example,
# the line below will run the model for 500 patients, deterministically, without treatment effects and with a 
# random seed equal to 177.
