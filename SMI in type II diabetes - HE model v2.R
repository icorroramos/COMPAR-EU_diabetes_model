######################################################################################################
########## Simulation model for self-management interventions (SMI) in type II diabetes ##############
######################################################################################################

# Original work Copyright (C) 2020 Isaac Corro Ramos & ... 
# Institute for Medical Technology Assessment (iMTA) of the Erasmus University Rotterdam.

# This is the R code of the health economic simulation model for SMI interventions in type II diabetes.
# This model is developed as part of the COMPAR-EU project (https://self-management.eu/)

####################################
########## TO DO LIST ##############
####################################

# Implement societal costs
# Implement treatment effects.
# Validation
# Implement PSA.


########################################
########## PART I - SETUP ##############
########################################

# For packages needed, etc. see: 
source("setup.R")

# For functions used in the simulation (e.g. UKPDS, informal care and productivity costs) see:
source("aux_functions.R")

##############################################
########## PART II - SIMULATION ##############
##############################################


# --> Work in progress...

init <- Sys.time()

# The simulation consists of one main function called "SMDMII_model_simulation". 
# This function is used to 1) simulate patients' clinical history, 2) calculate costs and 3) calculate QALYs. 

# Read patient characteristics: 
baseline_characteristics_MH2020 <- read.csv("input/baseline_characteristics_MH2020.csv", sep=",")

# These are only direct costs, based on Mount Hood challenge. Please change them to NL.
cost_inputs_MH2020 <- read.csv("input/cost_inputs_MH2020.csv", sep=",")

# Also, these utilities were used in Mount Hood challenge. For Dutch utilities see utilities.R. Please change the input file.
qol_inputs_MH2020 <- read.csv("input/qol_inputs_MH2020.csv", sep=",")

# NOTE: all comments below might not be 100% accurate. This is still work in progress and did not have time to check everything properly.

# At this moment the input parameters of the "SMDMII_model_simulation" function are the following:
# 1. patient_size_input: number of patients included in the simulation.
# 2. Gender binary indicator variable (1 == female) --> MH2020 separate analysis for males and females
# 3. Treatment effect parameters: 

# 11/05/2020: MH2020 Challenge
# Three separate interventions represented by a permanent decrement in common risk factors:
# (A) 10mm Hg reduction in Systolic Blood Pressure;
# (B) 0.5 mmol/l reduction in LDL Cholesterol (for models with LDL parameter can assume a reduction in total of 0.5 mmol/l and no change in HDL).
# (C) 1-unit reduction in BMI
# Reductions from these interventions should only be applied to post-baseline cycles and baseline values should remain unchanged.

# 4. discount rates for costs and effects 
# 5. run_PSA_input: 1 == runs the model in probabilistic mode, 0 == deterministic.
# 6. seed_input: A random seed that it is used to ensure consistency in the model results. 

SMDMII_model_simulation <- function(patient_size_input,
                                    female_input,
                                    tx_cost_input,
                                    treatment_effect_SBP_input,
                                    treatment_effect_LDL_input,
                                    treatment_effect_BMI_input,
                                    cost_disc_rate_input,
                                    qol_disc_rate_input,
                                    run_PSA_input,
                                    seed_input){
  
  # At this moment in the simulation "simulation_baseline_patients" is simply an empty data frame that will be used to store patient
  # characteristics at baseline.
  simulation_baseline_patients <- data.frame(matrix(vector(), 0, length(history_characteristics), dimnames=list(c(), history_characteristics)),stringsAsFactors=F)
  # Note: history_characteristics <- c("SIMID", sort(risk_factors_simulation), "SDURATION","dead") is defined in "aux_functions.R".
  
  ##################################################
  ########## MAIN PART I: simulate events ##########
  ##################################################
  
  ### Regression coefficients for PSA -- once per PSA --> To be done
  if(run_PSA_input == 1){} # end if regression coef. PSA
  
  ### IMPORTANT: Set a random seed to be able to replicate the results.
  ### This first seed is used to draw the same pool of patients when the function is called multiple times.
  ### Decide later if this has to go inside or outside the function. It depends on whether this seed is going to be 
  ### always fixed (inside) or if it can change (outside - and make it random or user defined)
  
  set.seed(seed_input) 
  
  ### Select the sample used in the simulation: sample with replacement from the data set.
  ### At this moment "baseline_characteristics_MH2020" consists of only two patients with exactly the same 
  ### characteristics, except for gender.
  simulation_baseline_patients[1:patient_size_input,colnames(baseline_characteristics_MH2020)] <- baseline_characteristics_MH2020[which(baseline_characteristics_MH2020$FEMALE==female_input),]
  
  ### Because at baseline there are less variables than in the simulation, we have to define first some of these new variables:
  
  # Age of diabetes diagnosis
  simulation_baseline_patients$AGE.DIAG <- simulation_baseline_patients$CURR.AGE - simulation_baseline_patients$YEAR 
  
  # Added 29/08/2020: age scaled parameters and informal care and productivity loss indicators
  simulation_baseline_patients$CURR.AGE.SCALE.INF  <- (simulation_baseline_patients$CURR.AGE - 73.75225)/9.938705 #hard coded values from data
  simulation_baseline_patients$CURR.AGE.SCALE.PROD <- (simulation_baseline_patients$CURR.AGE - 70.34077)/9.578899
  simulation_baseline_patients$CURR.AGE.2 <- (simulation_baseline_patients$CURR.AGE)^2
  simulation_baseline_patients$INF.CARE <- 0
  
  # Changed baseline age in Excel (from original 66 to 55) to have meaningful employment/productivity
  # Assume retirement age 65 for now: this should also be an input in the model!
  # If age >= 65 (retirement age), then we assume no productivity loss
  ifelse(simulation_baseline_patients$CURR.AGE >= 65, simulation_baseline_patients$EMPLOYED <- 0,
         {baseline_employed_prob <- apply(simulation_baseline_patients %>% select(risk_factors_employment), 1, function(x) annual_p_bernoulli(employment_equations$employment_coef,x)$p)
         simulation_baseline_patients$EMPLOYED <- unlist(lapply(baseline_employed_prob, function(x) rbinom(1,1,x))) #EMPLOYED = yes/no
         })
  simulation_baseline_patients$PROD.LOSS <- 0 # We set this to 0 at baseline and further updated according to employment status
  
  # BMI categories < 18.5 and >= 25 are hard-coded. If the definition of the categories changes, these have to change too.
  simulation_baseline_patients$BMI1 <- if_else(simulation_baseline_patients$BMI < 18.5, 1, 0)
  simulation_baseline_patients$BMI3 <- if_else(simulation_baseline_patients$BMI >= 25, 1, 0)
  
  # Initialize all .EVENT variables
  simulation_baseline_patients[event_vars] <- 0 
  
  # Mind the units for all continuous variables since some of them were re-scaled in the UKPDS equations. 
  # See description in aux_functions.R
  # Double-check the more/less variables!
  simulation_baseline_patients$SBP <- simulation_baseline_patients$SBP/100
  simulation_baseline_patients$HDL <- simulation_baseline_patients$HDL*10
  simulation_baseline_patients$LDL <- simulation_baseline_patients$LDL*10
  simulation_baseline_patients$LDL35more <- if_else(simulation_baseline_patients$LDL >= 35, simulation_baseline_patients$LDL, 0)
  simulation_baseline_patients$eGFR <- simulation_baseline_patients$eGFR/10 
  simulation_baseline_patients$eGFR60more <- if_else(simulation_baseline_patients$eGFR >= 6, simulation_baseline_patients$eGFR, 0)
  simulation_baseline_patients$eGFR60less <- if_else(simulation_baseline_patients$eGFR <  6, simulation_baseline_patients$eGFR, 0)
  simulation_baseline_patients$HEART.R <- simulation_baseline_patients$HEART.R/10 
  
  ### Create the simulation patient history table (for now is just empty)
  simulation_patients_history <- simulation_baseline_patients[FALSE,c(history_characteristics)]
  
  ### Choose the 1st patient
  patient_index <- 1
  
  ### Begin the loop on the simulation size (i.e. the number of patients you want to simulate)
  for(patient_index in 1:patient_size_input){
    
    # Print the patient index to know how advanced is the simulation.
    if(run_PSA_input == 0){print(patient_index)}
    
    # Pick the current patient from those selected at baseline and set simulation ID.
    current_patient <- simulation_baseline_patients[patient_index,]
    current_patient$SIMID <- patient_index
    
    # Initialize tracking variables for second events. If a patient has a second MI, a second stroke or 
    # a second amputation, it is not  possible to have a 3rd. These variables keep track of this.
    current_SMI_event <- 0
    current_SSTROKE_event <- 0
    current_AMP2_event <- 0
    
    # Save the characteristics to be used in the simulation history 
    # QUESTION: for the moment everything is saved (those that are stable too, but we could keep only those changing) 
    simulation_patients_history <- bind_rows(simulation_patients_history,current_patient[history_characteristics])
    
    ### Start the "timed" simulation (while loop = clock)
    factor_for_seed <- 100 #test = 1 to reproduce results. I normally used 100 but can be anything large enough
    
    while(current_patient$dead==0){
      
      ### These seeds will be used to draw the annual event probabilities while the patient is alive 
      if(run_PSA_input == 0){
        set.seed(factor_for_seed*current_patient$SIMID + current_patient$SDURATION)
      }else{
        set.seed(factor_for_seed*seed_input*current_patient$SIMID + current_patient$SDURATION)
      }
      
      #####################################
      # Sample annual event probabilities #
      #####################################
      
      # Assumptions: 
      #              1. Multiple events can occur in one year.
      #              2. Occurrence of events only affects death probability in the year they occur.
      #              3. Occurrence of events affects other events probability in the year after they occur (through history variables).
      
      # QUESTION: HOWEVER, IN THE DESCRIPTION OF THE MODEL (UKPDS PAPER) IT SAYS "RANDOMLY ORDER AND RUN EVENT EQUATIONS" -->>
      # ORDER WOULD IMPLY SEQUENTIAL RUN OF EQUATIONS AND UPDATES COULD BE BASED ON CURRENT YEAR. 
      # NO FURTHER EXPLANATION TO THIS IS FOUND 
      
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
      # The variables MI.EVENT and MI.HIST do not distinguish between first and second. However, the model assumption
      # is that no more than 2 MI events are possible. 
      if(current_patient$MI.HIST == 0){
        
        # If no history of MI, then it is first. Different for males and females
        if(current_patient$FEMALE == 1){
          # First MI for female Weibull
          current_FMIFEMALE_prob  <- annual_p_weibull(macrovascular_risk_equations$FMIFEMALE,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
          current_patient$MI.EVENT <- rbinom(1,1,current_FMIFEMALE_prob) # Update current_patient$MI.HIST after the year simulation is finished
          
        }else{
          # First MI for male Exponential
          current_FMIMALE_prob <- annual_p_weibull(macrovascular_risk_equations$FMIMALE,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
          current_patient$MI.EVENT <- rbinom(1,1,current_FMIMALE_prob) # Update current_patient$MI.HIST after the year simulation is finished
        } # end if/else for gender
        
      }else{
        # Second MI is Exponential, regardless the gender. MI.HIST will remain == 1.
        # If a patient has a second MI, it is not  possible to have a 3rd. 
        # The variable "current_SMI_event" keeps track of this. 
        if(current_SMI_event == 0){
          current_SMI_prob <- annual_p_weibull(macrovascular_risk_equations$SMI,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
          current_patient$MI.EVENT <- rbinom(1,1,current_SMI_prob)
          if(current_patient$MI.EVENT == 1){current_SMI_event <- 1}
        }
        }
      
      
      # STROKE could be first or second. 
      # The variables STROKE.EVENT and STROKE.HIST do not distinguish between first and second. However, the model assumption
      # is that no more than 2 STROKE events are possible. 
      if(current_patient$STROKE.HIST == 0){
        # If no history of STROKE, then it is first and Weibull. 
        current_FSTROKE_prob <- annual_p_weibull(macrovascular_risk_equations$FSTROKE,current_patient %>% select(risk_factors_macrovascular),current_patient$YEAR)$p
        current_patient$STROKE.EVENT <- rbinom(1,1,current_FSTROKE_prob) #Update current_patient$STROKE.HIST after the year
        
      }else{ 
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
        }else{
          # If prior ULCER then it is Exponential
          current_FAMPULCER_prob <- annual_p_weibull(microvascular_risk_equations$FAMPULCER,current_patient %>% select(risk_factors_microvascular),current_patient$YEAR)$p       
          current_patient$AMP1.EVENT <- rbinom(1,1,current_FAMPULCER_prob) #AMP.HIST to be updated at the end of the year
          }
        }else{ # if there is amputation history it can only be second and a 3rd one is not possible. 
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
      current_year_event <- sum(current_patient$CHF.EVENT, 
                                current_patient$IHD.EVENT, 
                                current_patient$MI.EVENT, # could be 1st or 2nd, no distinction
                                current_patient$STROKE.EVENT, # could be 1st or 2nd, no distinction
                                current_patient$BLIND.EVENT, 
                                current_patient$ULCER.EVENT, 
                                current_patient$AMP1.EVENT, 
                                current_patient$AMP2.EVENT,
                                current_patient$RENAL.EVENT)

      # If any event except blindness and ulcer happened in the current year, it should be captured with the following variable:
      current_year_event_no_blind_no_ulcer <- sum(current_patient$CHF.EVENT, 
                                                  current_patient$IHD.EVENT, 
                                                  current_patient$MI.EVENT, 
                                                  current_patient$STROKE.EVENT, 
                                                  current_patient$AMP1.EVENT, 
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
      
      
      ### INFORMAL CARE AND PRODUCTIVITY COSTS ###
      
      # Added 29/08/2020
      current_inf_care_prob <- annual_p_bernoulli(informal_care_equations$informal_care_coef, current_patient %>% select(risk_factors_informal))$p
      current_patient$INF.CARE <- rbinom(1,1,current_inf_care_prob) #INF.CARE = yes/no
      # If INF.CARE = yes, then calculate hours per day, and then total per year
      if(current_patient$INF.CARE == 1){
        current_patient$INF.CARE <- rgamma(1, 0.991532, 1/2.433387)*365.25 # Make parameters of the gamma distribution input parameters of the function
        }
      
      # The assumption is: only employed patients can get unemployed and incur productivity costs. After patients become unemployed,
      # they are assumed to continue like that in the simulation
      
      if(current_patient$EMPLOYED == 1){
        current_prod_loss_prob <- annual_p_bernoulli(informal_care_equations$prod_costs_coef, current_patient %>% select(risk_factors_prod))$p
        current_patient$PROD.LOSS <- rbinom(1,1,current_prod_loss_prob) # For the time being just show 0/1 but costs have to be calculated
        if(current_patient$PROD.LOSS == 1){
          current_patient$PROD.LOSS <- if_else(current_patient$FEMALE == 1, rgamma(1, 4.59, 1/4.99)*52, rgamma(1, 13.1, 1/2.59)*52) # Make parameters of the gamma distribution input parameters of the function
        }
        }
      
      ####################################
      # Update patient characteristics   #
      ####################################
      
      # We first copy all the previous characteristics
      current_patient_update <- current_patient # Do we need this? Or can we update directly current_patient?
      
      # Update first the history characteristics 
      current_patient_update$CHF.HIST    <- if_else(current_patient$CHF.EVENT + current_patient$CHF.HIST == 0, 0, 1)
      current_patient_update$IHD.HIST    <- if_else(current_patient$IHD.EVENT + current_patient$IHD.HIST == 0, 0, 1)
      current_patient_update$MI.HIST     <- if_else(current_patient$MI.EVENT + current_patient$MI.HIST == 0,0,1) 
      current_patient_update$STROKE.HIST <- if_else(current_patient$STROKE.EVENT + current_patient$STROKE.HIST == 0, 0, 1)
      current_patient_update$BLIND.HIST  <- if_else(current_patient$BLIND.EVENT + current_patient$BLIND.HIST == 0, 0, 1)
      current_patient_update$ULCER.HIST  <- if_else(current_patient$ULCER.EVENT + current_patient$ULCER.HIST == 0, 0, 1)
      current_patient_update$AMP.HIST    <- if_else(current_patient$AMP1.EVENT + current_patient$AMP2.EVENT + current_patient$AMP.HIST == 0, 0, 1)
      current_patient_update$RENAL.HIST  <- if_else(current_patient$RENAL.EVENT + current_patient$RENAL.HIST == 0, 0, 1)
      
      
      # Update risk factors  
      
      ### QUESTION: What remains unclear is whether all  risk factors change with time... and how...
      ### Decide what factors are assumed to be stable and which ones to change with time...
      ### For those that change, we need equations for predicting annual change based on 
      ### other patient characteristics. -->> WORK IN PROGRESS
      
      #########
      
      # Text below from MH2020 diabetes challenge
      # It is important in each simulation all other factors are kept constant between simulations 
      # and limit variation to the intervention costs as per instructions in the steps below. 
      # This includes assumptions around biomarker evolution; i.e. HbA1c and systolic blood pressure 
      # should be kept constant over time and not allowed to change over time (i.e., drift).
      
      #########
      
      
      current_patient_update$BMI <- current_patient$BMI #- treatment_effect_BMI_input
      
      # based on the above BMI, we should update BMI1 and BMI3
      current_patient_update$BMI1 <- if_else(current_patient_update$BMI < 18.5, 1, 0)
      current_patient_update$BMI3 <- if_else(current_patient_update$BMI >= 25, 1, 0)
      
      # MIND this because maybe we only update BMI after the first cycle.
      # We have to see how this is done for all treatment effects
      # If it changes only after the first cycle, we need to consider and index or time variable (SDURATION) to
      # update the treatment effect only in the first cycle and left it unchanged in the remaining of the simulation
      
      ## NOW the effect is applied to all cycles, which I think it's wrong
      
      current_patient_update$CURR.AGE <- current_patient$CURR.AGE + 1 
      current_patient_update$HbA1c <- current_patient$HbA1c  
      current_patient_update$LDL <- current_patient$LDL - treatment_effect_LDL_input #11/05/2020
      current_patient_update$SBP <- current_patient$SBP - treatment_effect_SBP_input #11/05/2020
      current_patient_update$YEAR <- current_patient$YEAR + 1
      current_patient_update$SDURATION <- current_patient$SDURATION + 1
      
      # Added 29/08/2020: age scaled parameters
      current_patient_update$CURR.AGE.SCALE.INF  <- (current_patient$CURR.AGE - 73.75225)/9.938705 #hardcoded values from data
      current_patient_update$CURR.AGE.SCALE.PROD <- (current_patient$CURR.AGE - 70.34077)/9.578899
      
      if(current_patient_update$PROD.LOSS != 0){current_patient_update$EMPLOYED <- 0}
      
      if(current_patient_update$CURR.AGE >= 65){
        current_patient_update$EMPLOYED <- 0
        current_patient_update$PROD.LOSS <- 0
        }
      
      
      # Force death after 100 years?
      if(current_patient_update$CURR.AGE > 100){current_patient_update$dead <- 1}
        
      # Force death if updated continuous risk factors take unfeasible value (e.g. too high hba1c)?
      # Consider logical bounds for update continuous attributes, e.g. > 0
      
      ### When all characteristics have been updated, we add these to the patient history
      simulation_patients_history <- bind_rows(simulation_patients_history,current_patient_update[history_characteristics])
      
      ### And update current patient and go up to while loop
      current_patient <- current_patient_update
      
      ### All the _event and .EVENT variables have to be reset to 0 because for the next year it only counts .HIST
      ### Note "event_vars" defined in aux_functions.R
      current_patient[event_vars] <- 0
      
    } #end while loop and move to another patient
    
    ### Update patient index
    patient_index <- patient_index + 1
    
  } #end for loop in number of patients
  
  ###################################################
  ########## MAIN PART II: Calculate Costs ##########
  ###################################################
  
  # All costs below are sourced from the MH2020 diabetes challenge. Discounting only applied to total costs only at this moment.
  cost_discount_factor <- ((1+cost_disc_rate_input)^(simulation_patients_history$SDURATION - 1))
  
  
  #Ischemic heart disease/Angina
  fatal_IHD_cost      <- cost_inputs_MH2020$IHD.FATAL*simulation_patients_history$IHD.EVENT*simulation_patients_history$dead
  nonfatal_IHD_cost   <- cost_inputs_MH2020$IHD.NONFATAL*simulation_patients_history$IHD.EVENT*(1-simulation_patients_history$dead)
  subsequent_IHD_cost <- cost_inputs_MH2020$IHD.SUB*simulation_patients_history$IHD.HIST # this is not correct now since it is added to the first year too!
  simulation_patients_history$IHD.COST <- fatal_IHD_cost + nonfatal_IHD_cost + subsequent_IHD_cost
  
  # Myocardial infarction
  fatal_MI_cost      <- cost_inputs_MH2020$MI.FATAL*simulation_patients_history$MI.EVENT*simulation_patients_history$dead
  nonfatal_MI_cost   <- cost_inputs_MH2020$MI.NONFATAL*simulation_patients_history$MI.EVENT*(1-simulation_patients_history$dead)
  subsequent_MI_cost <- cost_inputs_MH2020$MI.SUB*simulation_patients_history$MI.HIST
  simulation_patients_history$MI.COST <- fatal_MI_cost + nonfatal_MI_cost + subsequent_MI_cost
  
  # Heart failure
  fatal_CHF_cost      <- cost_inputs_MH2020$CHF.FATAL*simulation_patients_history$CHF.EVENT*simulation_patients_history$dead
  nonfatal_CHF_cost   <- cost_inputs_MH2020$CHF.NONFATAL*simulation_patients_history$CHF.EVENT*(1-simulation_patients_history$dead)
  subsequent_CHF_cost <- cost_inputs_MH2020$CHF.SUB*simulation_patients_history$CHF.HIST
  simulation_patients_history$CHF.COST <- fatal_CHF_cost + nonfatal_CHF_cost + subsequent_CHF_cost
  
  #Stroke
  fatal_STROKE_cost      <- cost_inputs_MH2020$STROKE.FATAL*simulation_patients_history$STROKE.EVENT*simulation_patients_history$dead
  nonfatal_STROKE_cost   <- cost_inputs_MH2020$STROKE.NONFATAL*simulation_patients_history$STROKE.EVENT*(1-simulation_patients_history$dead)
  subsequent_STROKE_cost <- cost_inputs_MH2020$STROKE.SUB*simulation_patients_history$STROKE.HIST
  simulation_patients_history$STROKE.COST <- fatal_STROKE_cost + nonfatal_STROKE_cost + subsequent_STROKE_cost
  
  # Amputation 	0	15,153	5,328	Alva et al. 2015 [1]
  # Note here we have in the model AMP1.EVENT and AMP2.EVENT but in terms of costs, there is no distinction
  # Question: in case of 2nd amp, should we consider the sum of the event + history of first?
  fatal_AMP_cost      <- cost_inputs_MH2020$AMP.FATAL*simulation_patients_history$AMP1.EVENT*simulation_patients_history$dead + cost_inputs_MH2020$AMP.FATAL*simulation_patients_history$AMP2.EVENT*simulation_patients_history$dead
  nonfatal_AMP_cost   <- cost_inputs_MH2020$AMP.NONFATAL*simulation_patients_history$AMP1.EVENT*(1-simulation_patients_history$dead) + cost_inputs_MH2020$AMP.NONFATAL*simulation_patients_history$AMP2.EVENT*(1-simulation_patients_history$dead)
  subsequent_AMP_cost <- cost_inputs_MH2020$AMP.SUB*simulation_patients_history$AMP.HIST
  simulation_patients_history$AMP.COST <- fatal_AMP_cost + nonfatal_AMP_cost + subsequent_AMP_cost
  
  #I'm not sure if costs for amputation are calculated ok. It should be consistent with other events when a 2nd event can occur.
  
  
  #Blindness
  fatal_BLIND_cost      <- cost_inputs_MH2020$BLIND.FATAL*simulation_patients_history$BLIND.EVENT*simulation_patients_history$dead
  nonfatal_BLIND_cost   <- cost_inputs_MH2020$BLIND.NONFATAL*simulation_patients_history$BLIND.EVENT*(1-simulation_patients_history$dead)
  subsequent_BLIND_cost <- cost_inputs_MH2020$BLIND.SUB*simulation_patients_history$BLIND.HIST
  simulation_patients_history$BLIND.COST <- fatal_BLIND_cost + nonfatal_BLIND_cost + subsequent_BLIND_cost
  
  # Ulcer
  fatal_ULCER_cost      <- cost_inputs_MH2020$ULCER.FATAL*simulation_patients_history$ULCER.EVENT*simulation_patients_history$dead
  nonfatal_ULCER_cost   <- cost_inputs_MH2020$ULCER.NONFATAL*simulation_patients_history$ULCER.EVENT*(1-simulation_patients_history$dead)
  subsequent_ULCER_cost <- cost_inputs_MH2020$ULCER.SUB*simulation_patients_history$ULCER.HIST
  simulation_patients_history$ULCER.COST <- fatal_ULCER_cost + nonfatal_ULCER_cost + subsequent_ULCER_cost
  
  # Cost in the absence of complications 	1,990			Alva et al. 2015 [1]
  # Not sure about this but at this moment applied only in the year of death
  simulation_patients_history$NOCOMP.COST <- cost_inputs_MH2020$NOCOMP*simulation_patients_history$dead
  
  # Intervention	Mean cost (£)
  # Blood glucose intervention 1: 
  #   (0.5%-point reduction in HbA1c) 	 60
  # Blood glucose intervention 2:
  #   (0.5%-point reduction in HbA1c)	120
  # Weight intervention 1: 
  #   (1-unit reduction in BMI (kg/m2))	55
  # Weight intervention 2:
  #   (1-unit reduction in BMI (kg/m2))	100
  
  # QUESTION: Not sure whether these are annual costs -->> should be but I find them quite low
  simulation_patients_history$TX.COST <- (tx_cost_input*(1-simulation_patients_history$dead) + tx_cost_input/2*simulation_patients_history$dead)
  
  # QUESTION: There is no half-cycle correction in these models -->> what are we assuming for the last year alive?
  # For example, in terms of treatment costs, do we assume a full year of costs? This assumption may apply to other costs too.
  # As it is now, we are assuming half-year costs for the year of death.
  
  
  ### Societal costs: Added 29/08/2020
  
  # Informal care and productivity loss: We calculated above INF.CARE & PROD.LOSS = 0/1 but costs have to be calculated here: 
  
  # informal care cost per hour assumed = 15. Make this an input parameter in the model
  simulation_patients_history$INF.CARE.COST <- 15*simulation_patients_history$INF.CARE*(1-simulation_patients_history$dead) + 15/2*simulation_patients_history$INF.CARE*simulation_patients_history$dead
  simulation_patients_history$PROD.LOSS.COST <- 15*simulation_patients_history$PROD.LOSS*(1-simulation_patients_history$dead) + 15/2*simulation_patients_history$PROD.LOSS*simulation_patients_history$dead
  
  # Need to create file for informal care and productivity costs as inputs: for now I used dummy numerical values to check that it works
   
  
  # Also these have to be added to total_costs below. Decide whether we want to have all together or distinguish between direct and societal.
  
  # Total annual costs
  total_costs <- (simulation_patients_history$IHD.COST + simulation_patients_history$MI.COST + 
                  simulation_patients_history$CHF.COST + simulation_patients_history$STROKE.COST +
                  simulation_patients_history$AMP.COST + simulation_patients_history$BLIND.COST + 
                  simulation_patients_history$ULCER.COST + simulation_patients_history$NOCOMP.COST + 
                  simulation_patients_history$TX.COST + simulation_patients_history$INF.CARE.COST + simulation_patients_history$PROD.LOSS.COST) 
  
  # Annual discounted costs
  simulation_patients_history$TOTAL.COST <- round((total_costs*(1-simulation_patients_history$dead) + (total_costs/2)*simulation_patients_history$dead)/cost_discount_factor,2)
  
  ####################################################
  ########## MAIN PART III: Calculate QALYs ##########
  ####################################################
  
  # All utilities below are sourced from the Mount Hood 2020 diabetes challenge.
  qaly_discount_factor <- ((1+qol_disc_rate_input)^(simulation_patients_history$SDURATION - 1))
  # Applied to total QALYs only at this moment.
  
  # Based on the 2018 Mt. Hood challenge, two suggestions were made:
  #   
  #   1)	The additive quality-of-life (QoL) model is recommended when populating the 
  #       health utility values into the simulation model:  
  #       If a subject has experienced two different complications belonging to 2 different categories 
  #       of disease (e.g., stroke [in the category of cerebrovascular disease] and myocardial infarction
  #       [in the category of coronary heart disease]), the health utility value will be reduced 
  #       by 0.219 which is the sum of individual decrement for these 2 complications (i.e., 0.164+0.055).
  #       However, if a subject has experienced several complications within the same category 
  #       of disease (e.g., myocardial infarction [in the category of coronary heart disease] and 
  #       congestive heart failure [in the category of coronary heart disease]), the health utility 
  #       value will be reduced by 0.108 (the decrement for heart failure) which is the largest 
  #       decrement of these two complications. 
  # 
  #   2)	The utility decrement and its 95% confidence interval for renal transplant was assumed to be 
  #       half of those for hemodialysis.
  
  # Utility values by categories of diseases/complications
  #
  # Disease category	Complication level Utility/Disutility Values	Lower 95% CI	Upper 95% CI
  #
  # Baseline utility value:	
  #                         T2DM without complications	0.785	0.681	0.889
  #
  # Acute metabolic disorder:	
  #                           Minor hypoglycemia event	-0.014	-0.004	-0.004
  #                           Major hypoglycemia event	-0.047	-0.012	-0.012
  #
  # Comorbidity:	
  #             Excess BMI (each unit above 25 kg/m2)	-0.006	-0.008	-0.004
  #
  # Retinopathy:
  #               Cataract	-0.016	-0.031	-0.001
  #               Moderate non-proliferative background diabetic retinopathy	-0.040	-0.066	-0.014
  #               Moderate macular edema	-0.040	-0.066	-0.014
  #               Vision-threatening diabetic retinopathy	-0.070	-0.099	-0.041
  #               Severe vision loss	-0.074	-0.124	-0.025
  #
  # Nephropathy:	
  #               Proteinuria	-0.048	-0.091	-0.005
  #               Renal transplant	-0.082	-0.137	-0.027
  #               Hemodialysis	-0.164	-0.274	-0.054
  #               Peritoneal dialysis	-0.204	-0.342	-0.066
  #
  # Neuropathy:	
  #             Peripheral vascular disease	-0.061	-0.090	-0.032
  #             Neuropathy	-0.084	-0.111	-0.057
  #             Active ulcer	-0.170	-0.207	-0.133
  #             Amputation event	-0.280	-0.389	-0.170
  #
  # Cerebrovascular disease:
  #                           Stroke	-0.164	-0.222	-0.105
  #
  # Coronary heart disease:	
  #                         Myocardial infarction	-0.055	-0.067	-0.042
  #                         Ischemic heart disease	-0.090	-0.126	-0.054
  #                         Heart failure	-0.108	-0.169	-0.048
  
  # QUESTION: similar to costs, there are more events defined in MH2020 than those included in our model.
  # Furthermore, it's unclear how to apply these decrements in time. For example, suppose a patient starts 
  # with no complications. The utility for the 1st year would be then 0.785. Suppose at year 2, this patient
  # experiences a stroke; then the utility for the second year would be 0.785 - 0.164. This is clear.
  # What happens in year 3 if there is no complications? Is the value again 0.785, 0.785 - 0.164 or anything else?
  # In other words, does the history of events impact the QoL over time? To me, it makes sense to say yes,
  # just think of a patient that has experienced several complications, it seems unrealistic to assume
  # that the QoL would be back to normal. But on the other hand, is it reasonable to assume that the 
  # disutility associated to the history of event is the same than the one associated to the event itself?
  # Probably not either...
  
  # The following event variables are at this moment included in the model
  
  # Coronary heart disease group:	
  simulation_patients_history$CHD.QALY <- pmin(simulation_patients_history$CHF.EVENT*qol_inputs_MH2020$CHF.EVENT, 
                                               simulation_patients_history$IHD.EVENT*qol_inputs_MH2020$IHD.EVENT, 
                                               simulation_patients_history$MI.EVENT*qol_inputs_MH2020$MI.EVENT)
  
  # Cerebrovascular disease: it's only stroke
  simulation_patients_history$STROKE.QALY <- simulation_patients_history$STROKE.EVENT*qol_inputs_MH2020$STROKE.EVENT
  
  # Neuropathy group: so far we have ULCER.EVENT, AMP1.EVENT and AMP2.EVENT
  simulation_patients_history$NEUROPATHY.QALY <- pmin(simulation_patients_history$ULCER.EVENT*qol_inputs_MH2020$ULCER.EVENT, 
                                                      simulation_patients_history$AMP1.EVENT*qol_inputs_MH2020$AMP.EVENT, 
                                                      simulation_patients_history$AMP2.EVENT*qol_inputs_MH2020$AMP.EVENT)
  
  # Retinopathy: so far we only have BLIND.EVENT 
  simulation_patients_history$BLIND.QALY <- simulation_patients_history$BLIND.EVENT*qol_inputs_MH2020$BLIND.EVENT
  
  # Nephropathy: so far we only have RENAL.EVENT: assumed disutility from haemodialysis but not sure
  simulation_patients_history$RENAL.QALY <- simulation_patients_history$RENAL.EVENT*qol_inputs_MH2020$RENAL.EVENT 
  
  # Comorbidity: Excess BMI (each unit above 25 kg/m2)	-0.006	-0.008	-0.004
  simulation_patients_history$BMI.QALY <- (simulation_patients_history$BMI - 25)*qol_inputs_MH2020$BMI.HIGH
  # Note: 25 in formula above is HARDCODED!
  
  
  # Total annual utilities
  total_utils <- (qol_inputs_MH2020$BASELINE  
                  + simulation_patients_history$CHD.QALY 
                  + simulation_patients_history$STROKE.QALY
                  + simulation_patients_history$NEUROPATHY.QALY
                  + simulation_patients_history$BLIND.QALY
                  + simulation_patients_history$RENAL.QALY
                  + simulation_patients_history$BMI.QALY)
  
  # Annual discounted QALYs
  simulation_patients_history$QALY <- round((total_utils*(1-simulation_patients_history$dead) + (total_utils/2)*simulation_patients_history$dead)/qaly_discount_factor,4)
  
  
  ################################################################
  ########## MAIN PART IV: Calculate aggregated results ##########
  ################################################################
  
  # To be done: at this moment, only clinical outcomes (e.g. event rates and life expectancy can be calculated).
  
  ### Life expectancy
  mean_life_expectancy <- round(mean(simulation_patients_history[which(simulation_patients_history$dead==1),"SDURATION"]),4)
  
  ### Events calculated differently depending on how were defined: .EVENT or .HIST
  
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
  
  total_qalys_patient <- aggregate(simulation_patients_history$QALY, list(Patient = simulation_patients_history$SIMID), sum)
  mean_total_qalys <- round(sum(total_qalys_patient$x)/patient_size_input, 2)
  
  
  ### Return model outcomes: at this moment only clinical history is saved and returned
  return(list(simulation_patients_history=simulation_patients_history,
              mean_life_expectancy = mean_life_expectancy,
              mean_CHF_rate = mean_CHF_rate,
              mean_BLIND_rate = mean_BLIND_rate,
              mean_ULCER_rate = mean_ULCER_rate,
              mean_AMP1_rate = mean_AMP1_rate,
              mean_AMP2_rate = mean_AMP2_rate,
              mean_MI_rate = mean_MI_rate,
              mean_RENAL_rate = mean_RENAL_rate,
              mean_total_costs = mean_total_costs,
              mean_total_qalys = mean_total_qalys))
  
} #end SMDMII_model_simulation function

# To run the model simply call the model function with the appropriate inputs in the correct order. For example,
# the line below will run the model for 500 patients, deterministic, without treatment effects and with a 
# random seed equal to 77.


sim_results_female <- SMDMII_model_simulation(5,      #patient_size_input
                                              1,       #female_input, 1 = female
                                              60,      #tx_cost_input
                                              10/100,  #treatment_effect_SBP_input from MH2020
                                              0.5/10,  #treatment_effect_LDL_input from MH2020
                                              1,       #treatment_effect_BMI_input from MH2020
                                              0.035,   #cost_disc_rate_input
                                              0.035,   #qol_disc_rate_input
                                              0,       #run_PSA_input, 0 == no PSA
                                              77       #seed_input
                                              )        


# sim_results_male <- SMDMII_model_simulation(5,      #patient_size_input
#                                             0,       #female_input, 1 = female
#                                             60,      #tx_cost_input
#                                             10/100,  #treatment_effect_SBP_input from MH2020
#                                             0.5/10,  #treatment_effect_LDL_input from MH2020
#                                             1,       #treatment_effect_BMI_input from MH2020
#                                             0.035,   #cost_disc_rate_input
#                                             0.035,   #qol_disc_rate_input
#                                             0,       #run_PSA_input, 0 == no PSA       
#                                             77       #seed_input
#                                             )        
# 
# sim_results_male$mean_life_expectancy
# sim_results_female$mean_life_expectancy
# 
# sim_results_male$mean_CHF_rate
# sim_results_female$mean_CHF_rate
# 
# sim_results_male$mean_BLIND_rate
# sim_results_female$mean_BLIND_rate
# 
# sim_results_male$mean_ULCER_rate
# sim_results_female$mean_ULCER_rate
# 
# sim_results_male$mean_AMP1_rate
# sim_results_female$mean_AMP1_rate
# 
# sim_results_male$mean_AMP2_rate
# sim_results_female$mean_AMP2_rate
# 
# sim_results_male$mean_MI_rate
# sim_results_female$mean_MI_rate # Not sure why I did not consider 2nd MI or 2nd Stroke but I did 2n amputation
# 
# sim_results_male$mean_RENAL_rate
# sim_results_female$mean_RENAL_rate
# 
# 
# 
# #View(sim_results_male$simulation_patients_history)
# #View(sim_results_female$simulation_patients_history)
# 
# ### To discuss:
# 
# ### MODEL UPDATES: 1. ASSUMPTIONS
# ###                2. AGGREGATED RESULTS (in progress --> discuss what we want to show: e.g. distinguish between first and second events?)  
# ###                3. BASELINE POPULATION (even though not needed for MH2020)
# 
# ### NEXT STEPS:
# ###                     1. Implement societal costs
# ###                     2. MH2020: continue implementing treatment effects
# ###                     3. HOW TO UPDATE RISK FACTORS: probably needed for MH2020
# 
# end <- Sys.time()
# end - init
