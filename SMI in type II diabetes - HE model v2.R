######################################################################################################
########## Simulation model for self-management interventions (SMI) in type II diabetes ##############
######################################################################################################

# Original work Copyright (C) 2020 Isaac Corro Ramos & the iMTA COMPAR-EU team. 
# Institute for Medical Technology Assessment (iMTA) of the Erasmus University Rotterdam.

# This is the R code of the health economic simulation model for SMI interventions in type II diabetes.
# This model is developed as part of the COMPAR-EU project (https://self-management.eu/)

########################################
########## PART I - SETUP ##############
########################################

# For packages needed, etc. see: 
source("setup.R")

# For functions, global lists, etc. used in the simulation (e.g. UKPDS, informal care and productivity costs) see:
source("aux_functions.R")

# Variable defined to keep track of simulation time (delete afterwards)
init <- Sys.time()

##############################################
########## PART II - SIMULATION ##############
##############################################

# The simulation consists of one main function called "SMDMII_model_simulation". 
# This function is used to 1) simulate patients' clinical history, 2) calculate costs and 3) calculate QALYs. 

# A patient-level model make use of patient characteristics. For the current version of the model, these are taken from
# the 2020 Mount Hood diabetes challenge. 
baseline_characteristics <- read.csv("input/baseline_characteristics_MH2020_LOLA.csv", sep=",")

# Direct costs of diabetes-related complications are also based on the Mount Hood challenge. Please change them to NL.
# For the moment, just converted to euro
cost_inputs <- read.csv("input/cost_inputs_MH2020_EURO.csv", sep=",")

# When a societal perspective is adopted, we also have future costs. This are obtained from the PAID online tool.
future_medical_cost_inputs <- read.csv("input/PAID_COMPAR_DM2_future_unrelated_costs_Unrelated_Costs_2020-09-03.csv", sep=",")
future_nonmedical_cost_inputs <- read.csv("input/PAID_COMPAR_DM2_future_unrelated_costs_Non-Medical_Consumption_2020-09-03.csv", sep=",")

# Also, these utilities were used as in the Mount Hood challenge. 
# For Dutch utilities see the file "utilities.R". The input file has been changes and Dutch utilities are used.
qol_inputs_NL <- read.csv("input/qol_inputs_MH2020_NL.csv", sep=",")

# At this moment the input parameters of the "SMDMII_model_simulation" function are the following:
# 1. patient_size_input: number of patients included in the simulation.
# 2. female_input: Gender binary indicator variable (1 == female) --> MH2020 separate analysis for males and females
# 3. tx_cost_input: treatment costs
# 4. Treatment effect parameters: all set to 0, so no treatment effect will be applied
# 5. discount rates for costs and effects 
# 6. run_PSA_input: 1 == runs the model in probabilistic mode, 0 == deterministic. Not working at the moment
# 7. seed_input: A random seed that it is used to ensure consistency in the model results. 

SMDMII_model_simulation <- function(patient_size_input, female_input, tx_cost_input,treatment_effect_SBP_input,
                                    treatment_effect_LDL_input,treatment_effect_BMI_input, cost_disc_rate_input,
                                    qol_disc_rate_input, run_PSA_input, seed_input){
  
  # At this moment "simulation_baseline_patients" is simply an empty data frame that will be used to store patient characteristics at baseline.
  simulation_baseline_patients <- data.frame(matrix(vector(), 0, length(history_characteristics), dimnames=list(c(), history_characteristics)),stringsAsFactors=F)
  # Note: the object "history_characteristics is defined in"the file "aux_functions.R".
  
  ##################################################
  ########## MAIN PART I: simulate events ##########
  ##################################################
  
  # Regression coefficients for PSA -- once per PSA --> To be done
  # if(run_PSA_input == 1){} # end if regression coef. PSA
  
  # IMPORTANT: Set a random seed to be able to replicate the results.
  set.seed(seed_input) 
  
  # Select the sample used in the simulation: sample with replacement from the data set.
  # At this moment "baseline_characteristics" consists of only two patients with exactly the same characteristics, except for gender.
  # Therefore, the sampling "functionality" is not really used here but it's kept for future versions of the model.
  simulation_baseline_patients[1:patient_size_input,colnames(baseline_characteristics)] <- baseline_characteristics[which(baseline_characteristics$FEMALE==female_input),]
  
  # Because at baseline there were less observed variables than those needed in the simulation, we defined first these new variables:
  # Age of diabetes diagnosis
  simulation_baseline_patients$AGE.DIAG <- simulation_baseline_patients$CURR.AGE - simulation_baseline_patients$YEAR 
  
  # Added 29/08/2020: age scaled parameters and informal care and productivity loss indicators
  simulation_baseline_patients$CURR.AGE.SCALE.INF  <- (simulation_baseline_patients$CURR.AGE - 73.75225)/9.938705 #hard-coded values from data
  simulation_baseline_patients$CURR.AGE.SCALE.PROD <- (simulation_baseline_patients$CURR.AGE - 70.34077)/9.578899
  simulation_baseline_patients$CURR.AGE.2 <- (simulation_baseline_patients$CURR.AGE)^2
  simulation_baseline_patients$INF.CARE <- 0
  
  # Note: changed baseline age in Excel (from original 66 to 55) to have meaningful employment/productivity
  
  # Assume retirement age 65 for now: this should also be an input in the model!
  # If age at baseline  >= 65 (retirement age), then we assume patient not employed and no productivity loss
  # Otherwise, age at baseline < 65, we calculate the probability of being employed (at baseline)
  ifelse(simulation_baseline_patients$CURR.AGE >= 65, simulation_baseline_patients$EMPLOYED <- 0,
         {baseline_employed_prob <- apply(simulation_baseline_patients %>% select(risk_factors_employment), 1, function(x) annual_p_bernoulli(employment_equations$employment_coef,x)$p)
         simulation_baseline_patients$EMPLOYED <- unlist(lapply(baseline_employed_prob, function(x) rbinom(1,1,x))) #EMPLOYED = yes/no
         })
  # We set productivity loss to 0 at BASELINE. This will be updated according to employment status as the simulation advances.
  simulation_baseline_patients$PROD.LOSS <- 0 
  
  # BMI categories 1 (BMI < 18.5) and 3 (BMI >= 25) are hard-coded. If the definition of the categories changes, these have to change too.
  simulation_baseline_patients$BMI1 <- if_else(simulation_baseline_patients$BMI < 18.5, 1, 0)
  simulation_baseline_patients$BMI3 <- if_else(simulation_baseline_patients$BMI >= 25, 1, 0)
  
  # Initialize all .EVENT variables. The object "event_vars" is defined in the file "aux_functions.R".
  simulation_baseline_patients[event_vars] <- 0 
  
  # Note: mind the units for all continuous variables since some of them were re-scaled in the UKPDS equations. Double-check the more/less variables!
  # See description in the file "aux_functions.R".
  simulation_baseline_patients$SBP <- simulation_baseline_patients$SBP/100
  simulation_baseline_patients$HDL <- simulation_baseline_patients$HDL*10
  simulation_baseline_patients$LDL <- simulation_baseline_patients$LDL*10
  simulation_baseline_patients$LDL35more <- if_else(simulation_baseline_patients$LDL >= 35, simulation_baseline_patients$LDL, 0)
  simulation_baseline_patients$eGFR <- simulation_baseline_patients$eGFR/10 
  simulation_baseline_patients$eGFR60more <- if_else(simulation_baseline_patients$eGFR >= 6, simulation_baseline_patients$eGFR, 0)
  simulation_baseline_patients$eGFR60less <- if_else(simulation_baseline_patients$eGFR <  6, simulation_baseline_patients$eGFR, 0)
  simulation_baseline_patients$HEART.R <- simulation_baseline_patients$HEART.R/10 
  
  # Create now the simulation patient history table (for now is just empty) to save all simulation results 
  simulation_patients_history <- simulation_baseline_patients[FALSE,c(history_characteristics)]
  
  # Initialize the index for the patients entering the for loop 
  patient_index <- 1
  
  # Begin the loop on the simulation size (i.e. the number of patients to be simulated)
  for(patient_index in 1:patient_size_input){
    
    # Print the patient index to know how advanced is the simulation. Delete later if not needed
    if(run_PSA_input == 0){print(patient_index)}
    
    # Pick the current patient from those selected at baseline and set simulation ID ("SIMID"). This is needed to produce aggregated results.
    current_patient <- simulation_baseline_patients[patient_index,]
    current_patient$SIMID <- patient_index
    
    # Initialize tracking variables for second events. If a patient has a second MI, a second stroke or a second amputation, it is not  possible to have a 3rd. 
    # These variables keep track of this.
    current_SMI_event <- 0
    current_SSTROKE_event <- 0
    current_AMP2_event <- 0
    
    # Save the characteristics to be used in the simulation history 
    simulation_patients_history <- bind_rows(simulation_patients_history,current_patient[history_characteristics])
    
    # Start the "timed" simulation (while loop = clock)
    # Set random seed large enough. Will be used for drawing event probabilities below.
    factor_for_seed <- 100 
    
    while(current_patient$dead==0){
      
      # These seeds will be used to draw the annual event probabilities while the patient is alive 
      if(run_PSA_input == 0){set.seed(factor_for_seed*current_patient$SIMID + current_patient$SDURATION)}
      else{set.seed(factor_for_seed*seed_input*current_patient$SIMID + current_patient$SDURATION)}
      
      #####################################
      # Sample annual event probabilities #
      #####################################
      
      # Assumptions: 
      #              1. Multiple events may occur in one year.
      #              2. Occurrence of events only affects death probability in the year they occur.
      #              3. Occurrence of events affects other events probability in the year after they occur (through history variables).
      
      # QUESTION: HOWEVER, IN THE DESCRIPTION OF THE MODEL (UKPDS PAPER) IT SAYS "RANDOMLY ORDER AND RUN EVENT EQUATIONS" -->>
      # ORDER WOULD IMPLY SEQUENTIAL RUN OF EQUATIONS AND UPDATES COULD BE BASED ON CURRENT YEAR. NO FURTHER EXPLANATION TO THIS IS FOUND 
      
      ### MACROVASCULAR COMPLICATIONS ###
      
      # Heart Failure is Weibull. This can happen only once; that's why the if condition below is used. 
      if(current_patient$CHF.HIST == 0){
        current_CHF_prob  <- annual_p_weibull(macrovascular_risk_equations$CHF,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
        current_patient$CHF.EVENT <- rbinom(1,1,current_CHF_prob) 
      } # Update current_patient$CHF.HIST after the year simulation is finished
      
      # IHD is Weibull. This can happen only once; that's why the if condition below is used.
      if(current_patient$IHD.HIST == 0){
        current_IHD_prob  <- annual_p_weibull(macrovascular_risk_equations$IHD,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
        current_patient$IHD.EVENT <- rbinom(1,1,current_IHD_prob) 
      } # Update current_patient$IHD.HIST after the year simulation is finished
      
      # MI could be first or second. If it is first, then it is different for male and female. 
      # If it is second, then it is the same for both genders.
      # The variables MI.EVENT and MI.HIST do not distinguish between first and second. 
      # However, the model assumption is that no more than 2 MI events are possible. 
      if(current_patient$MI.HIST == 0){
        
        # If no history of MI, then it is first. Different for males and females
        if(current_patient$FEMALE == 1){
          # First MI for female Weibull
          current_FMIFEMALE_prob  <- annual_p_weibull(macrovascular_risk_equations$FMIFEMALE,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
          current_patient$MI.EVENT <- rbinom(1,1,current_FMIFEMALE_prob) # Update current_patient$MI.HIST after the year simulation is finished
        }
        else{
          # First MI for male Exponential
          current_FMIMALE_prob <- annual_p_weibull(macrovascular_risk_equations$FMIMALE,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
          current_patient$MI.EVENT <- rbinom(1,1,current_FMIMALE_prob) # Update current_patient$MI.HIST after the year simulation is finished
        } # end if/else for gender
      }
      else{
        # Second MI is Exponential, regardless the gender. MI.HIST will remain == 1. If a patient has a second MI, it is not  possible to have a 3rd. 
        # The variable "current_SMI_event" keeps track of this. 
        if(current_SMI_event == 0){
          current_SMI_prob <- annual_p_weibull(macrovascular_risk_equations$SMI,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
          current_patient$MI.EVENT <- rbinom(1,1,current_SMI_prob)
          if(current_patient$MI.EVENT == 1){current_SMI_event <- 1}
        }
      }
      
      # STROKE could be first or second. The variables STROKE.EVENT and STROKE.HIST do not distinguish between first and second. 
      # However, the model assumption is that no more than 2 STROKE events are possible. 
      if(current_patient$STROKE.HIST == 0){
        # If no history of STROKE, then it is first and Weibull. 
        current_FSTROKE_prob <- annual_p_weibull(macrovascular_risk_equations$FSTROKE,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
        current_patient$STROKE.EVENT <- rbinom(1,1,current_FSTROKE_prob) #Update current_patient$STROKE.HIST after the year
      }
      else{ 
        #Second STROKE is Weibull. Here current_patient$STROKE.HIST == 1 and remains like that.
        #If a patient has a second STROKE, it is not  possible to have a 3rd. "current_SSTROKE_event" keeps track of this. Not sure whether it has to be initialised. Seems to be 0, which is ok.
        if(current_SSTROKE_event == 0){
          current_SSTROKE_prob <- annual_p_weibull(macrovascular_risk_equations$SSTROKE,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
          current_patient$STROKE.EVENT <- rbinom(1,1,current_SSTROKE_prob) 
          if(current_patient$STROKE.EVENT == 1){current_SSTROKE_event <- 1} 
        }
      }
      
      ### MICROVASCULAR COMPLICATIONS ###
      
      # BLINDNESS is Exponential. This is assumed to happen only once; that's why the if condition below is used. 
      if(current_patient$BLIND.HIST == 0){
        current_BLIND_prob <- annual_p_weibull(microvascular_risk_equations$BLIND,current_patient %>% select(risk_factors_microvascular),current_patient$YEAR)$p       
        current_patient$BLIND.EVENT <- rbinom(1,1,current_BLIND_prob) #Update current_patient$BLIND.HIST at the end of the year. 
      }
      
      # ULCER is Logistic. This is assumed to happen only once; that's why the if condition below is used. 
      if(current_patient$ULCER.HIST == 0){
        current_ULCER_prob  <- annual_p_logistic(microvascular_risk_equations$BLIND,current_patient %>% select(risk_factors_microvascular))$p
        current_patient$ULCER.EVENT <- rbinom(1,1,current_ULCER_prob) # Update current_patient$ULCER.HIST at the end of the year.
      }
      
      # AMPUTATION could be first or second. First amputation depends on ULCER history. 
      if(current_patient$AMP.HIST == 0){
        # If no history of AMPUTATION, then it is first and depends on ULCER.
        if(current_patient$ULCER.HIST == 0){
          # If no prior ULCER then it is Weibull
          current_FAMPNOULCER_prob  <- annual_p_weibull(microvascular_risk_equations$FAMPNOULCER,current_patient %>% select(risk_factors_microvascular),current_patient$YEAR)$p       
          current_patient$AMP1.EVENT <- rbinom(1,1,current_FAMPNOULCER_prob) #AMP.HIST to be updated at the end of the year and does not distinguishes between 1st and 2nd
        }
        else{
          # If prior ULCER then it is Exponential
          current_FAMPULCER_prob <- annual_p_weibull(microvascular_risk_equations$FAMPULCER,current_patient %>% select(risk_factors_microvascular),current_patient$YEAR)$p       
          current_patient$AMP1.EVENT <- rbinom(1,1,current_FAMPULCER_prob) #AMP.HIST to be updated at the end of the year
        }
      }
      else{ # if there is amputation history it can only be second and a 3rd one is not possible. 
        # Second amputation is exponential: AMP.HIST no need to be updated. But also a patient cannot have more than 2 amputations.
        if(current_AMP2_event == 0){
          current_SAMP_prob <- annual_p_weibull(microvascular_risk_equations$SAMP,current_patient %>% select(risk_factors_microvascular),current_patient$YEAR)$p       
          current_patient$AMP2.EVENT <- rbinom(1,1,current_SAMP_prob) 
          if(current_patient$AMP2.EVENT == 1){current_AMP2_event <- 1}
        }
      }
      
      # Renal failure is Exponential. This is assumed to happen only once; that's why the if condition below is used. 
      if(current_patient$RENAL.HIST == 0){
        current_RENALF_prob <- annual_p_weibull(microvascular_risk_equations$RENALF,current_patient %>% select(risk_factors_microvascular),current_patient$YEAR)$p       
        current_patient$RENAL.EVENT <- rbinom(1,1,current_RENALF_prob) #current_patient$RENAL.HIST updated after the year. Check here: I got really large values for the probability so have a second look
      }
      
      # There are four equations for death, depending on the events occurred in the current year and the history of previous events. 
      # Assumption (already mentioned above): occurrence of events only affects death probability in the year they occur.
      # The four equations are for:
      #                             1. Years with no history of previous events and no events in the current year
      #                             2. First year of events (so no previous history) excluding blindness or ulcer
      #                             3. Years with history of previous events but no events in the current year
      #                             4. Subsequent years (so there is previous history) of events excluding blindness or ulcer
      # Therefore, we need to define variables to determine what equation should be used.  
      
      # If any event happened in the current year, it should be captured with the following variable: 
      current_year_event <- sum(current_patient$CHF.EVENT, current_patient$IHD.EVENT, 
                                current_patient$MI.EVENT, # could be 1st or 2nd, no distinction
                                current_patient$STROKE.EVENT, # could be 1st or 2nd, no distinction
                                current_patient$BLIND.EVENT, current_patient$ULCER.EVENT, current_patient$AMP1.EVENT, 
                                current_patient$AMP2.EVENT, current_patient$RENAL.EVENT)
      
      # If any event except blindness and ulcer happened in the current year, it should be captured with the following variable:
      current_year_event_no_blind_no_ulcer <- sum(current_patient$CHF.EVENT, current_patient$IHD.EVENT, current_patient$MI.EVENT, 
                                                  current_patient$STROKE.EVENT, current_patient$AMP1.EVENT, current_patient$AMP2.EVENT, 
                                                  current_patient$RENAL.EVENT)
      
      # If current_year_event - current_year_event_no_blind_no_ulcer > 0 it means that either blindness or ulcer occurred in the current year.
      
      # .HIST variables are not updated for the current year. But the ones from the previous year are captured in this variable:
      current_hist  <- sum(current_patient$CHF.HIST, current_patient$IHD.HIST, current_patient$MI.HIST, 
                           current_patient$STROKE.HIST, current_patient$BLIND.HIST, current_patient$ULCER.HIST,
                           current_patient$AMP.HIST, current_patient$RENAL.HIST)
      
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
      
      # Sampling "dead" status
      current_DEATH_event <- rbinom(1,1,current_DEATH_prob) 
      current_patient$dead <- current_DEATH_event
      
      ### INFORMAL CARE AND PRODUCTIVITY COSTS: Added 29/08/2020
      
      # INFORMAL CARE
      # Probability of receiving at least weekly informal care for a whole year is calculated using a Bernoulli distribution.
      current_inf_care_prob <- annual_p_bernoulli(informal_care_equations$informal_care_coef, current_patient %>% select(risk_factors_informal))$p
      current_patient$INF.CARE <- rbinom(1,1,current_inf_care_prob) #INF.CARE = yes/no
      # If INF.CARE = yes, then calculate hours per day, and then total per year: Gamma distribution for now but other options are possible.
      # Make parameters of the gamma distribution input parameters of the function
      if(current_patient$INF.CARE == 1){current_patient$INF.CARE <- rgamma(1, 0.991532, 1/2.433387)*365.25}
      
      # PRODUCTIVITY LOSS -- main assumptions: 
      # 1. Short-term productivity loss costs: only for EMPLOYED patients based on sick days --> Assume fix 14 working days based on data. This could be changed.
      # 2. Permanent (one-off) productivity loss costs: only employed patients can get unemployed and incur productivity costs the year after they get unemployed. Friction method.
      # 3. After patients become unemployed, they are assumed to continue like that in the simulation --> No more productivity loss costs
      
      # If currently employed, then apply first short-term productivity costs
      if(current_patient$EMPLOYED == 1){
        # Short-term productivity costs (sick days)
        # Calculate number of worked hours per week: assume min 1, max 40, different for male and females
        # Gamma distribution assumed but could be something else. Make parameters of the gamma distribution input parameters of the function
        current_worked_hours_week <- if_else(current_patient$FEMALE == 1, min(max(1,rgamma(1, 4.59, 1/4.99)),40), min(40,max(1,rgamma(1, 13.1, 1/2.59)))) 
        current_worked_hours_day <- current_worked_hours_week/5
        # Assume 14 working days per year lost due to sickness as mentioned above and use it to calculate short-term prod. loss costs.
        sick_hours_year <- current_worked_hours_day*14 
        # Cost of an hour sick is different for male and females. Make input parameters of the function.
        cost_hour_sick <- if_else(current_patient$FEMALE == 1, 33.74, 40.46)
        # Multiply number of sick hours by the costs of an hour
        current_patient$PROD.LOSS <- sick_hours_year*cost_hour_sick
        
        # Permanent (one-off) prod. loss costs
        # Calculate first the probability of losing job (Bernoulli distribution).
        current_jobless_prob <- annual_p_bernoulli(informal_care_equations$prod_costs_coef, current_patient %>% select(risk_factors_prod))$p
        current_jobless <- rbinom(1,1,current_jobless_prob) # For the time being just show 0/1 but costs have to be calculated
        # If jobless, then apply permanent cost (one-off) cost: friction method as proportion of a maximum of 85 days. MAke 85 input parameter of the model.
        if(current_jobless == 1){ 
          # Add permanent cost (one-off) cost to cost of sick days
          current_patient$PROD.LOSS <- current_patient$PROD.LOSS + cost_hour_sick*85*current_worked_hours_day/8
          # Update employed status
          current_patient$EMPLOYED <- 0
        }
      }
      
      ####################################
      # Update patient characteristics   #
      ####################################
      
      # We first copy all the previous characteristics
      current_patient_update <- current_patient # Check: Do we need this? Or can we update directly current_patient?
      
      # Update the history characteristics 
      current_patient_update$CHF.HIST    <- if_else(current_patient$CHF.EVENT + current_patient$CHF.HIST == 0, 0, 1)
      current_patient_update$IHD.HIST    <- if_else(current_patient$IHD.EVENT + current_patient$IHD.HIST == 0, 0, 1)
      current_patient_update$MI.HIST     <- if_else(current_patient$MI.EVENT + current_patient$MI.HIST == 0,0,1) 
      current_patient_update$STROKE.HIST <- if_else(current_patient$STROKE.EVENT + current_patient$STROKE.HIST == 0, 0, 1)
      current_patient_update$BLIND.HIST  <- if_else(current_patient$BLIND.EVENT + current_patient$BLIND.HIST == 0, 0, 1)
      current_patient_update$ULCER.HIST  <- if_else(current_patient$ULCER.EVENT + current_patient$ULCER.HIST == 0, 0, 1)
      current_patient_update$AMP.HIST    <- if_else(current_patient$AMP1.EVENT + current_patient$AMP2.EVENT + current_patient$AMP.HIST == 0, 0, 1)
      current_patient_update$RENAL.HIST  <- if_else(current_patient$RENAL.EVENT + current_patient$RENAL.HIST == 0, 0, 1)
      
      # Update risk factors  
      
      #################################################################################################################
      # QUESTION: What remains unclear is whether all  risk factors change with time... and how...                    #
      # Decide what factors are assumed to be stable and which ones to change with time...                            #
      # For those that change, we need equations for predicting annual change based on other patient characteristics. #
      # Text below from MH2020 diabetes challenge:                                                                    #
      # It is important in each simulation all other factors are kept constant between simulations                    #
      # and limit variation to the intervention costs as per instructions in the steps below.                         #
      # This includes assumptions around biomarker evolution; i.e. HbA1c and systolic blood pressure                  #
      # should be kept constant over time and not allowed to change over time (i.e., drift).                          # 
      #################################################################################################################
      
      # BMI
      current_patient_update$BMI <- current_patient$BMI #- treatment_effect_BMI_input # treatment effect currently removed  
      # Based on the above BMI, we should update BMI1 and BMI3
      current_patient_update$BMI1 <- if_else(current_patient_update$BMI < 18.5, 1, 0)
      current_patient_update$BMI3 <- if_else(current_patient_update$BMI >= 25, 1, 0)
      
      current_patient_update$CURR.AGE <- current_patient$CURR.AGE + 1 
      current_patient_update$HbA1c <- current_patient$HbA1c  
      current_patient_update$LDL <- current_patient$LDL - treatment_effect_LDL_input #11/05/2020
      current_patient_update$SBP <- current_patient$SBP - treatment_effect_SBP_input #11/05/2020
      current_patient_update$YEAR <- current_patient$YEAR + 1
      current_patient_update$SDURATION <- current_patient$SDURATION + 1
      
      # Added 29/08/2020: age scaled parameters
      current_patient_update$CURR.AGE.SCALE.INF  <- (current_patient$CURR.AGE - 73.75225)/9.938705 #hardcoded values from data
      current_patient_update$CURR.AGE.SCALE.PROD <- (current_patient$CURR.AGE - 70.34077)/9.578899
      
      if(current_patient_update$CURR.AGE >= 65){
        current_patient_update$EMPLOYED <- 0
        current_patient_update$PROD.LOSS <- 0}
      
      # Force death at 100 years
      if(current_patient_update$CURR.AGE >= 99){current_patient_update$dead <- 1}
      
      # Note: Force death if updated continuous risk factors take unfeasible value (e.g. too high HbA1c)?
      # Consider logical bounds for update continuous attributes, e.g. > 0
      
      # When all characteristics are updated, we add these to the patient history
      simulation_patients_history <- bind_rows(simulation_patients_history,current_patient_update[history_characteristics])
      
      # And update current patient and go up to while loop
      current_patient <- current_patient_update
      
      # All the _event and .EVENT variables have to be reset to 0 because for the next year it only counts .HIST
      # Note "event_vars" defined in aux_functions.R
      current_patient[event_vars] <- 0
      if(current_patient$EMPLOYED == 0){current_patient$PROD.LOSS <- 0} # needed here or somewhere else? It seems to work though
    } #end while loop and move to another patient
    
    # Update patient index
    patient_index <- patient_index + 1
  } #end for loop in number of patients
  
  ###################################################
  ########## MAIN PART II: Calculate Costs ##########
  ###################################################
  
  # All cost inputs were sourced from the MH2020 diabetes challenge. We need NL costs.
  cost_discount_factor <- ((1+cost_disc_rate_input)^(simulation_patients_history$SDURATION - 1))
  
  #Ischemic heart disease/Angina
  fatal_IHD_cost      <- cost_inputs$IHD.FATAL*simulation_patients_history$IHD.EVENT*simulation_patients_history$dead
  nonfatal_IHD_cost   <- cost_inputs$IHD.NONFATAL*simulation_patients_history$IHD.EVENT*(1-simulation_patients_history$dead)
  subsequent_IHD_cost <- cost_inputs$IHD.SUB*simulation_patients_history$IHD.HIST 
  # Note: this is not completely correct now since subsequent costs are added to the first year too! Check this happens for all complications?
  simulation_patients_history$IHD.COST <- (fatal_IHD_cost + nonfatal_IHD_cost + subsequent_IHD_cost)/cost_discount_factor
  
  # Myocardial infarction
  fatal_MI_cost      <- cost_inputs$MI.FATAL*simulation_patients_history$MI.EVENT*simulation_patients_history$dead
  nonfatal_MI_cost   <- cost_inputs$MI.NONFATAL*simulation_patients_history$MI.EVENT*(1-simulation_patients_history$dead)
  subsequent_MI_cost <- cost_inputs$MI.SUB*simulation_patients_history$MI.HIST
  simulation_patients_history$MI.COST <- (fatal_MI_cost + nonfatal_MI_cost + subsequent_MI_cost)/cost_discount_factor
  
  # Heart failure
  fatal_CHF_cost      <- cost_inputs$CHF.FATAL*simulation_patients_history$CHF.EVENT*simulation_patients_history$dead
  nonfatal_CHF_cost   <- cost_inputs$CHF.NONFATAL*simulation_patients_history$CHF.EVENT*(1-simulation_patients_history$dead)
  subsequent_CHF_cost <- cost_inputs$CHF.SUB*simulation_patients_history$CHF.HIST
  simulation_patients_history$CHF.COST <- (fatal_CHF_cost + nonfatal_CHF_cost + subsequent_CHF_cost)/cost_discount_factor
  
  #Stroke
  fatal_STROKE_cost      <- cost_inputs$STROKE.FATAL*simulation_patients_history$STROKE.EVENT*simulation_patients_history$dead
  nonfatal_STROKE_cost   <- cost_inputs$STROKE.NONFATAL*simulation_patients_history$STROKE.EVENT*(1-simulation_patients_history$dead)
  subsequent_STROKE_cost <- cost_inputs$STROKE.SUB*simulation_patients_history$STROKE.HIST
  simulation_patients_history$STROKE.COST <- (fatal_STROKE_cost + nonfatal_STROKE_cost + subsequent_STROKE_cost)/cost_discount_factor
  
  # Amputation 	0	15,153	5,328	Alva et al. 2015 [1]
  # Note here we have in the model AMP1.EVENT and AMP2.EVENT but in terms of costs, there is no distinction
  # Question: in case of 2nd amp, should we consider the sum of the event + history of first?
  # Not sure whether costs for amputation are calculated ok. It should be consistent with other events when a 2nd event can occur.
  fatal_AMP_cost      <- cost_inputs$AMP.FATAL*simulation_patients_history$AMP1.EVENT*simulation_patients_history$dead + cost_inputs$AMP.FATAL*simulation_patients_history$AMP2.EVENT*simulation_patients_history$dead
  nonfatal_AMP_cost   <- cost_inputs$AMP.NONFATAL*simulation_patients_history$AMP1.EVENT*(1-simulation_patients_history$dead) + cost_inputs$AMP.NONFATAL*simulation_patients_history$AMP2.EVENT*(1-simulation_patients_history$dead)
  subsequent_AMP_cost <- cost_inputs$AMP.SUB*simulation_patients_history$AMP.HIST
  simulation_patients_history$AMP.COST <- (fatal_AMP_cost + nonfatal_AMP_cost + subsequent_AMP_cost)/cost_discount_factor
  
  #Blindness
  fatal_BLIND_cost      <- cost_inputs$BLIND.FATAL*simulation_patients_history$BLIND.EVENT*simulation_patients_history$dead
  nonfatal_BLIND_cost   <- cost_inputs$BLIND.NONFATAL*simulation_patients_history$BLIND.EVENT*(1-simulation_patients_history$dead)
  subsequent_BLIND_cost <- cost_inputs$BLIND.SUB*simulation_patients_history$BLIND.HIST
  simulation_patients_history$BLIND.COST <- (fatal_BLIND_cost + nonfatal_BLIND_cost + subsequent_BLIND_cost)/cost_discount_factor
  
  # Ulcer
  fatal_ULCER_cost      <- cost_inputs$ULCER.FATAL*simulation_patients_history$ULCER.EVENT*simulation_patients_history$dead
  nonfatal_ULCER_cost   <- cost_inputs$ULCER.NONFATAL*simulation_patients_history$ULCER.EVENT*(1-simulation_patients_history$dead)
  subsequent_ULCER_cost <- cost_inputs$ULCER.SUB*simulation_patients_history$ULCER.HIST
  simulation_patients_history$ULCER.COST <- (fatal_ULCER_cost + nonfatal_ULCER_cost + subsequent_ULCER_cost)/cost_discount_factor
  
  # Cost in the absence of complications 	1,990	Alva et al. 2015 [1]
  # Not sure about this but at this moment applied only in the year of death
  simulation_patients_history$NOCOMP.COST <- (cost_inputs$NOCOMP*simulation_patients_history$dead)/cost_discount_factor
  
  # Intervention costs: Not sure whether these are annual costs -->> should be but I find them quite low
  # For LOLA, this might be done differently, so set to 0.
  simulation_patients_history$TX.COST <- ((tx_cost_input*(1-simulation_patients_history$dead) + tx_cost_input/2*simulation_patients_history$dead))/cost_discount_factor
  
  # QUESTION: There is no half-cycle correction in these models -->> what are we assuming for the last year alive?
  # For example, in terms of treatment costs, do we assume a full year of costs? This assumption may apply to other costs too.
  # As it is now, we are assuming half-year costs for the year of death.
  
  
  ### Societal costs: Added 29/08/2020
  
  # Informal care and productivity loss: We calculated above INF.CARE & PROD.LOSS = 0/1 but costs have to be calculated here: 
  
  # Informal care cost per hour assumed = 14.95. Make this an input parameter in the model.
  simulation_patients_history$INF.CARE.COST <- (14.95*simulation_patients_history$INF.CARE*(1-simulation_patients_history$dead) + 14.95/2*simulation_patients_history$INF.CARE*simulation_patients_history$dead)/cost_discount_factor
  simulation_patients_history$PROD.LOSS.COST <- (simulation_patients_history$PROD.LOSS*(1-simulation_patients_history$dead) + 1/2*simulation_patients_history$PROD.LOSS*simulation_patients_history$dead)/cost_discount_factor
  
  # Future costs are dependent on age, gender and alive status: add explanation later
  future_cost_matrix <- inner_join(simulation_patients_history[c("CURR.AGE","FEMALE","dead")],future_medical_cost_inputs, by = 'CURR.AGE')
  future_cost_matrix <- inner_join(future_cost_matrix,future_nonmedical_cost_inputs, by = 'CURR.AGE')
  
  simulation_patients_history$FUTURE.COST.MEDICAL <- ((1-future_cost_matrix$FEMALE)*(1 - future_cost_matrix$dead)*future_cost_matrix[,6] + 
                                                        (1-future_cost_matrix$FEMALE)*future_cost_matrix$dead*future_cost_matrix[,4] + 
                                                        future_cost_matrix$FEMALE*(1 - future_cost_matrix$dead)*future_cost_matrix[,7] + 
                                                        future_cost_matrix$FEMALE*future_cost_matrix$dead*future_cost_matrix[,5])/cost_discount_factor 
  
  simulation_patients_history$FUTURE.COST.NONMEDICAL <- future_cost_matrix[,8]/cost_discount_factor 
  
  # Total annual costs discounted
  simulation_patients_history$TOTAL.COST <- (simulation_patients_history$IHD.COST + simulation_patients_history$MI.COST + simulation_patients_history$CHF.COST + 
                                               simulation_patients_history$STROKE.COST + simulation_patients_history$AMP.COST + simulation_patients_history$BLIND.COST + 
                                               simulation_patients_history$ULCER.COST + simulation_patients_history$NOCOMP.COST + simulation_patients_history$TX.COST + 
                                               simulation_patients_history$INF.CARE.COST + simulation_patients_history$PROD.LOSS.COST + 
                                               simulation_patients_history$FUTURE.COST.MEDICAL + simulation_patients_history$FUTURE.COST.NONMEDICAL) 
  
  ####################################################
  ########## MAIN PART III: Calculate QALYs ##########
  ####################################################
  
  # All utilities were sourced from the Mount Hood 2020 diabetes challenge and converted to NL values as explained above.
  # Discounting is applied to total QALYs only at this moment.
  qaly_discount_factor <- ((1+qol_disc_rate_input)^(simulation_patients_history$SDURATION - 1))
  
  ######################################################################################################################
  # QUESTION: similar to costs, there are more events defined in MH2020 than those included in our model.              #
  # Furthermore, it's unclear how to apply these decrements in time. For example, suppose a patient starts             #
  # with no complications. The utility for the 1st year would be then 0.785. Suppose at year 2, this patient           #
  # experiences a stroke; then the utility for the second year would be 0.785 - 0.164. This is clear.                  #
  # What happens in year 3 if there is no complications? Is the value again 0.785, 0.785 - 0.164 or anything else?     #
  # In other words, does the history of events impact the QoL over time? To me, it makes sense to say yes,             #
  # just think of a patient that has experienced several complications, it seems unrealistic to assume                 #
  # that the QoL would be back to normal. But on the other hand, is it reasonable to assume that the                   #
  # disutility associated to the history of event is the same than the one associated to the event itself?             #
  # Probably not either...                                                                                             #
  ######################################################################################################################
  
  # The following event variables are at this moment included in the model
  
  # Coronary heart disease group:	
  simulation_patients_history$CHD.QALY <- pmin(simulation_patients_history$CHF.EVENT*qol_inputs_NL$CHF.EVENT, 
                                               simulation_patients_history$IHD.EVENT*qol_inputs_NL$IHD.EVENT, 
                                               simulation_patients_history$MI.EVENT*qol_inputs_NL$MI.EVENT)
  
  # Cerebrovascular disease: only stroke
  simulation_patients_history$STROKE.QALY <- simulation_patients_history$STROKE.EVENT*qol_inputs_NL$STROKE.EVENT
  
  # Neuropathy group: so far we have ULCER.EVENT, AMP1.EVENT and AMP2.EVENT
  simulation_patients_history$NEUROPATHY.QALY <- pmin(simulation_patients_history$ULCER.EVENT*qol_inputs_NL$ULCER.EVENT, 
                                                      simulation_patients_history$AMP1.EVENT*qol_inputs_NL$AMP.EVENT, 
                                                      simulation_patients_history$AMP2.EVENT*qol_inputs_NL$AMP.EVENT)
  
  # Retinopathy: so far we only have BLIND.EVENT 
  simulation_patients_history$BLIND.QALY <- simulation_patients_history$BLIND.EVENT*qol_inputs_NL$BLIND.EVENT
  
  # Nephropathy: so far we only have RENAL.EVENT: assumed disutility from haemodialysis but not sure
  simulation_patients_history$RENAL.QALY <- simulation_patients_history$RENAL.EVENT*qol_inputs_NL$RENAL.EVENT 
  
  # Comorbidity: Excess BMI (each unit above 25 kg/m2)	-0.006	-0.008	-0.004
  simulation_patients_history$BMI.QALY <- (simulation_patients_history$BMI - 25)*qol_inputs_NL$BMI.HIGH
  # Note: 25 in formula above is HARDCODED!
  
  # Total annual utilities
  total_utils <- (qol_inputs_NL$BASELINE + simulation_patients_history$CHD.QALY + simulation_patients_history$STROKE.QALY
                  + simulation_patients_history$NEUROPATHY.QALY + simulation_patients_history$BLIND.QALY + simulation_patients_history$RENAL.QALY
                  + simulation_patients_history$BMI.QALY)
  
  # Annual discounted QALYs
  simulation_patients_history$QALY <- round((total_utils*(1-simulation_patients_history$dead) + (total_utils/2)*simulation_patients_history$dead)/qaly_discount_factor,4)
  
  ################################################################
  ########## MAIN PART IV: Calculate aggregated results ##########
  ################################################################
  
  # To be completed: 
  
  # Life expectancy
  mean_life_expectancy <- round(mean(simulation_patients_history[which(simulation_patients_history$dead==1),"SDURATION"]),4)
  
  # Event rates: note events calculated differently depending on how were defined: .EVENT or .HIST
  # Total number of events per patient during lifetime and rate per year
  cum_CHF <- aggregate(simulation_patients_history$CHF.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_CHF_rate <- round(sum(cum_CHF$x)/patient_size_input, 4) 
  
  cum_BLIND <- aggregate(simulation_patients_history$BLIND.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_BLIND_rate <- round(sum(cum_BLIND$x)/patient_size_input, 4)
  
  cum_ULCER <- aggregate(simulation_patients_history$ULCER.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_ULCER_rate <- round(sum(cum_ULCER$x)/patient_size_input, 4) 
  
  cum_AMP1 <- aggregate(simulation_patients_history$AMP1.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_AMP1_rate <- round(sum(cum_AMP1$x)/patient_size_input, 4) 
  
  cum_AMP2 <- aggregate(simulation_patients_history$AMP2.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_AMP2_rate <- round(sum(cum_AMP2$x)/patient_size_input, 4) 
  
  cum_MI <- aggregate(simulation_patients_history$MI.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_MI_rate <- round(sum(cum_MI$x)/patient_size_input, 4) 
  
  cum_RENAL <- aggregate(simulation_patients_history$RENAL.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_RENAL_rate <- round(sum(cum_RENAL$x)/patient_size_input, 4)#round(mean(cum_RENAL$x)/mean_life_expectancy, 4) # note this could be per patient or per time
  
  # Total costs per patient during lifetime and average per patient 
  total_costs_patient <- aggregate(simulation_patients_history$TOTAL.COST, list(Patient = simulation_patients_history$SIMID), sum)
  mean_total_costs <- round(sum(total_costs_patient$x)/patient_size_input, 2)
  
  # Breakdown costs: 
  
  ihd_costs_patient <- aggregate(simulation_patients_history$IHD.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_ihd_costs <- round(sum(ihd_costs_patient$x)/patient_size_input, 2)
  
  mi_costs_patient <- aggregate(simulation_patients_history$MI.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_mi_costs <- round(sum(mi_costs_patient$x)/patient_size_input, 2)
  
  chf_costs_patient <- aggregate(simulation_patients_history$CHF.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_chf_costs <- round(sum(chf_costs_patient$x)/patient_size_input, 2)
  
  stroke_costs_patient <- aggregate(simulation_patients_history$STROKE.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_stroke_costs <- round(sum(stroke_costs_patient$x)/patient_size_input, 2)
  
  amp_costs_patient <- aggregate(simulation_patients_history$AMP.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_amp_costs <- round(sum(amp_costs_patient$x)/patient_size_input, 2)
  
  blind_costs_patient <- aggregate(simulation_patients_history$BLIND.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_blind_costs <- round(sum(blind_costs_patient$x)/patient_size_input, 2)
  
  ulcer_costs_patient <- aggregate(simulation_patients_history$ULCER.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_ulcer_costs <- round(sum(ulcer_costs_patient$x)/patient_size_input, 2)
  
  mean_complication_costs <- mean_ihd_costs + mean_mi_costs + mean_chf_costs + mean_stroke_costs + mean_amp_costs + mean_blind_costs + mean_ulcer_costs
  
  # No complication costs
  
  nocomp_costs_patient <- aggregate(simulation_patients_history$NOCOMP.COST  , list(Patient = simulation_patients_history$SIMID), sum)
  mean_nocomp_costs <- round(sum(nocomp_costs_patient$x)/patient_size_input, 2)
  
  # To do, if relevant: simulation_patients_history$TX.COST 
  
  inf_care_costs_patient <- aggregate(simulation_patients_history$INF.CARE.COST, list(Patient = simulation_patients_history$SIMID), sum)
  mean_inf_care_costs <- round(sum(inf_care_costs_patient$x)/patient_size_input, 2)
  
  prod_loss_costs_patient <- aggregate(simulation_patients_history$PROD.LOSS.COST, list(Patient = simulation_patients_history$SIMID), sum)
  mean_prod_loss_costs <- round(sum(prod_loss_costs_patient$x)/patient_size_input, 2)
  
  future_medical_costs_patient <- aggregate(simulation_patients_history$FUTURE.COST.MEDICAL, list(Patient = simulation_patients_history$SIMID), sum)
  mean_future_medical_costs <- round(sum(future_medical_costs_patient$x)/patient_size_input, 2)
  
  future_nonmedical_costs_patient <- aggregate(simulation_patients_history$FUTURE.COST.NONMEDICAL, list(Patient = simulation_patients_history$SIMID), sum)
  mean_future_nonmedical_costs <- round(sum(future_nonmedical_costs_patient$x)/patient_size_input, 2)
  
  
  total_qalys_patient <- aggregate(simulation_patients_history$QALY, list(Patient = simulation_patients_history$SIMID), sum)
  mean_total_qalys <- round(sum(total_qalys_patient$x)/patient_size_input, 2)
  
  ### Return model outcomes: at this moment only clinical history is saved and returned
  ### Do not return clinical outcomes for the moment, except life expectancy, since these need to be validated.
  return(list(simulation_patients_history=simulation_patients_history, 
              mean_life_expectancy = mean_life_expectancy,
              #mean_CHF_rate = mean_CHF_rate, mean_BLIND_rate = mean_BLIND_rate, mean_ULCER_rate = mean_ULCER_rate,
              #mean_AMP1_rate = mean_AMP1_rate, mean_AMP2_rate = mean_AMP2_rate, mean_MI_rate = mean_MI_rate, mean_RENAL_rate = mean_RENAL_rate,
              mean_total_costs = mean_total_costs, mean_total_qalys = mean_total_qalys,
              mean_inf_care_costs = mean_inf_care_costs, mean_prod_loss_costs = mean_prod_loss_costs,
              mean_ihd_costs = mean_ihd_costs, mean_mi_costs = mean_mi_costs, mean_chf_costs = mean_chf_costs,
              mean_stroke_costs = mean_stroke_costs, mean_amp_costs = mean_amp_costs, mean_blind_costs = mean_blind_costs,
              mean_ulcer_costs = mean_ulcer_costs, mean_complication_costs = mean_complication_costs,
              mean_nocomp_costs = mean_nocomp_costs, mean_future_medical_costs = mean_future_medical_costs, 
              mean_future_nonmedical_costs = mean_future_nonmedical_costs))
  
} #end SMDMII_model_simulation function

# To run the model call the model function with the appropriate inputs in the correct order. 
sim_results_female <- SMDMII_model_simulation(5,     #patient_size_input: run 500 for LOLA
                                              1,       #female_input, 1 = female
                                              0,       #tx_cost_input
                                              0,       #treatment_effect_SBP_input from MH2020
                                              0,       #treatment_effect_LDL_input from MH2020
                                              0,       #treatment_effect_BMI_input from MH2020
                                              0.040,   #cost_disc_rate_input
                                              0.015,   #qol_disc_rate_input
                                              0,       #run_PSA_input, 0 == no PSA
                                              77       #seed_input
)        

sim_results_male <- SMDMII_model_simulation(5,     #patient_size_input: run 500 for LOLA
                                            0,       #female_input, 1 = female
                                            0,       #tx_cost_input
                                            0,       #treatment_effect_SBP_input from MH2020
                                            0,       #treatment_effect_LDL_input from MH2020
                                            0,       #treatment_effect_BMI_input from MH2020
                                            0.040,   #cost_disc_rate_input
                                            0.015,   #qol_disc_rate_input
                                            0,       #run_PSA_input, 0 == no PSA
                                            77       #seed_input
) 

# Results tables
sim_results_female_table <- matrix(c(sim_results_female$mean_complication_costs,sim_results_female$mean_nocomp_costs,
                                     sim_results_female$mean_inf_care_costs, sim_results_female$mean_prod_loss_costs, 
                                     sim_results_female$mean_future_medical_costs, sim_results_female$mean_future_nonmedical_costs,
                                     sim_results_female$mean_total_costs, sim_results_female$mean_total_qalys), nrow = 1)

colnames(sim_results_female_table) <- c("Complication costs", "No complication costs", "Informal care costs", "Productivity costs", 
                                        "Future medical costs", "Future non-medical costs", "Total costs", "Total QALYs")
rownames(sim_results_female_table) <- "Intervention"
sim_results_female_table

sim_results_male_table <- matrix(c(sim_results_male$mean_complication_costs,sim_results_male$mean_nocomp_costs,
                                   sim_results_male$mean_inf_care_costs, sim_results_male$mean_prod_loss_costs, 
                                   sim_results_male$mean_future_medical_costs, sim_results_male$mean_future_nonmedical_costs,
                                   sim_results_male$mean_total_costs, sim_results_male$mean_total_qalys), nrow = 1)

colnames(sim_results_male_table) <- colnames(sim_results_female_table)
rownames(sim_results_male_table) <- rownames(sim_results_female_table)
sim_results_male_table

# Variable defined to keep track of simulation time (delete afterwards)
end <- Sys.time()
end - init

####################################
########## TO DO LIST ##############
####################################
# Complication costs --> NL: For the moment just converted to EURO
# Patient age at baseline: 66 years in MH2020. What's more relevant for us? Younger patients? 55 YEARS IN INPUT FILE LOLA
# How to add costs and QALYs for the first year: Gimon to add

# FOR THE FUTURE:
# Implement treatment effects.
# Implement PSA.
# Validation.
# Improve code efficiency
# AGGREGATED RESULTS (in progress --> discuss what we want to show: e.g. distinguish between first and second events?)  
# BASELINE POPULATION (even though not needed for MH2020)
# HOW TO UPDATE RISK FACTORS: probably not needed for MH2020
