#_____________________________________________________________________________________________________
#
########## Simulation model for self-management interventions (SMI) in type II diabetes ##############
#_____________________________________________________________________________________________________

# Original work Copyright (C) 2020 Isaac Corro Ramos, Gimon de Graaf & the iMTA COMPAR-EU team. 
# Institute for Medical Technology Assessment (iMTA) of the Erasmus University Rotterdam.

# This is the R code of the health economic simulation model for SMI interventions in type II diabetes.
# This model is developed as part of the COMPAR-EU project (https://self-management.eu/)


# PART I - SETUP ----------------------------------------------------------


# For packages needed, etc. see: 
source("setup.R")

# For functions, global lists, etc. used in the simulation (e.g. UKPDS, informal care and productivity costs) see:
source("R/aux_functions.R")

# The following objects should be present in the environment where the simulation function is called
# baseline_characteristics
# male_cost_inputs
# female_cost_inputs
# future_medical_cost_inputs
# future_nonmedical_cost_inputs
# qol_inputs
# qol_events_inputs

# TODO: annotate the objects above on their requirements of format/ content

# PART II - SIMULATION ----------------------------------------------------

# The simulation consists of one main function called "SMDMII_model_simulation". 
# This function is used to 1) simulate patients' clinical history, 2) calculate costs and 3) calculate QALYs. 

SMDMII_model_simulation <- function(patient_size_input, # numeric value > 0, patient_size_input: number of patients included in the simulation.
                                    female_input, # 1 = females, 0 = male
                                    tx_cost_input, # treatment costs -->> Need to implement
                                    treatment_effect_HbA1c_input, # vector including c(effect, start effect, end effect, start decline)
                                    treatment_effect_HDL_input, 
                                    treatment_effect_LDL_input,
                                    treatment_effect_BMI_input, 
                                    cost_disc_rate_input, # discount rates for costs and effects 
                                    qol_disc_rate_input, 
                                    retirement_age_input, # retirement age
                                    run_PSA_input, # run_PSA_input: 1 == runs the model in probabilistic mode, 0 == deterministic. Not implemented at the moment
                                    seed_input){ # A random seed that it is used to ensure consistency in the model results. 
  
  # At this moment "simulation_baseline_patients" is simply an empty data frame that will be used to store patient characteristics at baseline.
  simulation_baseline_patients <- data.frame(matrix(vector(), 0, length(history_characteristics), dimnames=list(c(), history_characteristics)),stringsAsFactors=F)
  # Note: the object "history_characteristics is defined in"the file "aux_functions.R".
  

  ########## MAIN PART I: simulate events ##########

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
  simulation_baseline_patients$CURR.AGE.SCALE.INF  <- (simulation_baseline_patients$CURR.AGE - 72.5474088)/10.4626624 
  # Hard-coded values from data: Primary analysis central countries assumed for UK
  simulation_baseline_patients$CURR.AGE.SCALE.PROD <- (simulation_baseline_patients$CURR.AGE - 60.2737989)/6.1177269 
  # Hard-coded values from data: Primary analysis central countries assumed for UK
  
  simulation_baseline_patients$CURR.AGE.2 <- (simulation_baseline_patients$CURR.AGE)^2
  simulation_baseline_patients$INF.CARE   <- 0
  
  # If age at baseline  >= retirement age, then we assume patient not employed and no productivity loss
  # Otherwise, age at baseline < retirement age, we calculate the probability of being employed (at baseline)
  ifelse(simulation_baseline_patients$CURR.AGE >= retirement_age_input, simulation_baseline_patients$EMPLOYED <- 0,
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
  
  # eGFR
  simulation_baseline_patients$eGFR       <- simulation_baseline_patients$eGFR/10 
  simulation_baseline_patients$eGFR60more <- if_else(simulation_baseline_patients$eGFR >= 6, simulation_baseline_patients$eGFR, 0)
  simulation_baseline_patients$eGFR60less <- if_else(simulation_baseline_patients$eGFR <  6, simulation_baseline_patients$eGFR, 0)
  
  # HDL
  simulation_baseline_patients$HDL <- simulation_baseline_patients$HDL*10
  
  # Heart rate
  simulation_baseline_patients$HEART.R <- simulation_baseline_patients$HEART.R/10 
  
  # LDL 
  simulation_baseline_patients$LDL       <- simulation_baseline_patients$LDL*10
  simulation_baseline_patients$LDL35more <- if_else(simulation_baseline_patients$LDL >= 35, simulation_baseline_patients$LDL, 0)
  
  # SBP
  simulation_baseline_patients$SBP <- simulation_baseline_patients$SBP/10
  
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
      
      # Sample annual event probabilities #
      
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
        #print(current_FSTROKE_prob)
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
        current_ULCER_prob  <- annual_p_logistic(microvascular_risk_equations$ULCER,current_patient %>% select(risk_factors_microvascular))$p #typo corrected
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
      
      # Note .HIST variables are not updated for the current year. But the ones from the previous year are captured in this variable:
      current_hist  <- sum(current_patient$CHF.HIST, current_patient$IHD.HIST, current_patient$MI.HIST, 
                           current_patient$STROKE.HIST, current_patient$BLIND.HIST, current_patient$ULCER.HIST,
                           current_patient$AMP.HIST, current_patient$RENAL.HIST)
      
      # The four equations for death are then the following: 
      
      # 1. If no history of previous events and no events in the current year, then gompertz distribution
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
      #print(current_DEATH_prob)
      current_patient$dead <- current_DEATH_event
      
      ### INFORMAL CARE AND PRODUCTIVITY COSTS: Added 29/08/2020
      
      # INFORMAL CARE
      # Probability of receiving at least weekly informal care for a whole year is calculated using a Bernoulli distribution.
      current_inf_care_prob <- annual_p_bernoulli(informal_care_equations$informal_care_coef, current_patient %>% select(risk_factors_informal))$p
      current_patient$INF.CARE <- rbinom(1,1,current_inf_care_prob) #INF.CARE = yes/no
      
      # If INF.CARE = yes, then calculate hours per day, and then total per year: Gamma distribution for now but other options are possible.
      # Make parameters of the gamma distribution input parameters of the function
      # Gamma distribution updated 10/12/2020: primary analysis central countries
      if(current_patient$INF.CARE == 1){current_patient$INF.CARE <- rgamma(1, 1.753152, 1/1.022368)*365.25}
      # Note here: Irene suggested we might better use the mean/median
      
      # PRODUCTIVITY LOSS -- main assumptions: 
      # 1. Short-term productivity loss costs: only for EMPLOYED patients based on sick days --> Assume fix 14 working days based on data. This could be changed.
      # 2. Permanent (one-off) productivity loss costs: only employed patients can get unemployed and incur productivity costs the year after they get unemployed. Friction method.
      # 3. After patients become unemployed, they are assumed to continue like that in the simulation --> No more productivity loss costs
      
      # If currently employed, then apply first short-term productivity costs
      if(current_patient$EMPLOYED == 1){
        # Short-term productivity costs (sick days)
        # Calculate number of worked hours per week: assume mean value, different for male and females. See e-mail 10/12/2020.
        # This could be updated if we are able to fit a distribution.
        current_worked_hours_week <- if_else(current_patient$FEMALE == 1, 26.5096273, 33.8016304) 
        current_worked_hours_day <- current_worked_hours_week/5
        
        # Assumption 11/12/2020: working days per year lost due to sickness --> Apply median value for central European countries.
        # The median was chosen as opposed to the mean due to the skewness of the distribution (mean seems to be very large).
        # There is no probability distribution fit to the data. Therefore, values are fixed. This could be updated if we manage to fit.
        # The model distinguishes between days lost due to 1) diabetes and 2) diabetes with heart attack/stroke. The latter is 
        # applied at the year of the event only. 
        # Values are 14 days for diabetes and 27.5 for diabetes + CHF/stroke
        
        # The definition of this variable in SHARE is: 'a heart attack including myocardial infarction or coronary thrombosis or 
        # any other heart problem including congestive heart failure'. 
        # So it could either be an MI event, or onset of heart failure. I therefore think these sick days should apply to both patients 
        # with congestive heart failure and patients with a myocardial infarction, but I think also to patients with ischemic heart disease. 
        
        sick_hours_year <- if_else(current_patient$STROKE.EVENT+current_patient$MI.EVENT+current_patient$CHF.EVENT+current_patient$IHD.EVENT>0, 
                                   current_worked_hours_day*27.5, 
                                   current_worked_hours_day*14)   
        
        # Cost of an hour sick: based the productivity costs per hour on Eurostat: for the UK it is 28.50 Euros per hour (not specified by gender) in 2019.
        # https://ec.europa.eu/eurostat/databrowser/view/tps00173/default/table?lang=en.
        # Inflated to 2020 is 28.93 -->> https://www.officialdata.org/europe/inflation/2019?amount=28.50
        # And converted into GBP is 26.6369 -->> https://www.xe.com/currencyconverter/convert/?Amount=28.93&From=EUR&To=GBP
        cost_hour_sick <- 26.6369
        # Multiply number of sick hours by the costs of an hour
        current_patient$PROD.LOSS <- sick_hours_year*cost_hour_sick
        
        # Permanent (one-off) prod. loss costs
        # Calculate first the probability of losing job (Bernoulli distribution).
        current_jobless_prob <- annual_p_bernoulli(informal_care_equations$prod_costs_coef, current_patient %>% select(risk_factors_prod))$p
        current_jobless <- rbinom(1,1,current_jobless_prob) 
        
        # If jobless, then apply permanent cost (one-off) cost: friction method as proportion of a maximum of 2.7 months (re-calculated as days).
        # To do: make 2.7 months an input parameter of the model. We can assume a friction period in the UK of 2.7 months based on this publication:
        #http://pure-oai.bham.ac.uk/ws/files/40800976/Kigozi_et_al_2017_Health_Economics.pdf
        
        if(current_jobless == 1){ 
          # Add permanent cost (one-off) cost to cost of sick days
          current_patient$PROD.LOSS <- current_patient$PROD.LOSS + cost_hour_sick*(2.7*365.25/12)*current_worked_hours_day/8
          # Update employed status
          current_patient$EMPLOYED <- 0
        }
      }
      
      # Update patient characteristics   #
      
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
      ### 23/11/2020: we are distinguishing between risk factors that will be updated in every loop iteration (1 year)
      ### and those depending on treatment effect. 
      ### For treatment effect parameters we are going to use two additional parameters: one to indicate the total duration
      ### of the treatment effect and another one to indicate when a decline on the effect will start.
      ### For example, if the total duration is 10 years, we will assume that after 10 years there is no more tx effect.
      ### If the decline starts at year 5, then we will assume full effect for years 1 to 5, and decline from 6 to 10.
      ### We need to decide how the decline will be modeled. We may consider several options. 
      ### Parameters that remain at their baseline value through the simulation will not be shown here.
      
      ### No Tx effect parameters
      current_patient_update$SDURATION <- current_patient$SDURATION + 1
      
      current_patient_update$CURR.AGE <- current_patient$CURR.AGE + 1 
      current_patient_update$YEAR <- current_patient$YEAR + 1 # Note YEAR is duration of diabetes
      
      # Added 29/08/2020: age scaled parameters
      current_patient_update$CURR.AGE.SCALE.INF  <- (current_patient$CURR.AGE - 73.75225)/9.938705 #hardcoded values from data
      current_patient_update$CURR.AGE.SCALE.PROD <- (current_patient$CURR.AGE - 70.34077)/9.578899 #hardcoded values from data
      
      if(current_patient_update$CURR.AGE >= retirement_age_input){
        current_patient_update$EMPLOYED <- 0
        current_patient_update$PROD.LOSS <- 0
        }
      
      # Force death at 99 years: DECIDE THIS TOO
      if(current_patient_update$CURR.AGE >= 99){current_patient_update$dead <- 1}
      
      ### Tx effect parameters
      
      # HbA1c
      
      # Note: unique(simulation_baseline_patients$HbA1c) has the HbA1c baseline value. Note we have at baseline only 2 patients
      # with the same characteristics except one is male and the other female. 
      # We need "unique" because simulation_baseline_patients is already a matrix with the size of the simulation 
      
      if(current_patient_update$SDURATION == treatment_effect_HbA1c_input[2]){
        current_patient_update$HbA1c <- min(max(5, current_patient$HbA1c + treatment_effect_HbA1c_input[1]),18) #23/11/2020
        }
      
      if(current_patient_update$SDURATION %in% treatment_effect_HbA1c_input[4]:treatment_effect_HbA1c_input[3]){
        current_patient_update$HbA1c <- min(current_patient$HbA1c - (treatment_effect_HbA1c_input[1])/(treatment_effect_HbA1c_input[3]-treatment_effect_HbA1c_input[4]),unique(simulation_baseline_patients$HbA1c)) #23/11/2020
        }
      
      
      # BMI
      current_patient_update$BMI <- current_patient$BMI #- treatment_effect_BMI_input # treatment effect currently removed  
      # Based on the above BMI, we should update BMI1 and BMI3
      current_patient_update$BMI1 <- if_else(current_patient_update$BMI < 18.5, 1, 0)
      current_patient_update$BMI3 <- if_else(current_patient_update$BMI >= 25, 1, 0)
      
      # LDL: same approach as for HbA1c14/12/2020
      # QUESTION: natural bounds for LDL???
      
      if(current_patient_update$SDURATION == treatment_effect_LDL_input[2]){
        current_patient_update$LDL <- min(max(0, current_patient$LDL + treatment_effect_LDL_input[1]),50)
      }
      if(current_patient_update$SDURATION %in% treatment_effect_LDL_input[4]:treatment_effect_LDL_input[3]){
        current_patient_update$LDL <- min(current_patient$LDL - (treatment_effect_LDL_input[1])/(treatment_effect_LDL_input[3]-treatment_effect_LDL_input[4]),unique(simulation_baseline_patients$LDL)) 
      }
      
      # HDL: same approach as for HbA1c14/12/2020
      # QUESTION: natural bounds for HDL???
      
      if(current_patient_update$SDURATION == treatment_effect_HDL_input[2]){
        current_patient_update$HDL <- min(max(0, current_patient$HDL + treatment_effect_HDL_input[1]),50)
      }
      if(current_patient_update$SDURATION %in% treatment_effect_HDL_input[4]:treatment_effect_HDL_input[3]){
        current_patient_update$HDL <- max(current_patient$HDL-(treatment_effect_HDL_input[1])/(treatment_effect_HDL_input[3]-treatment_effect_HDL_input[4]),unique(simulation_baseline_patients$HDL)) 
        
        # CAREFUL HERE: I put max instead of min above because the treatment effect is positive and not negative.
        # try to find a general solution, maybe using absolute value
        
      }
      
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
  
  
  ########## MAIN PART II: Calculate Costs ##########
  
  # All cost inputs were sourced from the MH2020 diabetes challenge. We need NL costs.
  cost_discount_factor <- (1+cost_disc_rate_input)^(simulation_patients_history$SDURATION)
  
  # Complication costs are dependent on age, gender and alive status
  complication_cost_inputs <- data.frame(ifelse(unique(simulation_patients_history$FEMALE) == 0, male_cost_inputs,female_cost_inputs) )
  
  if(unique(simulation_patients_history$FEMALE) == 0){complication_cost_inputs <- male_cost_inputs}else{complication_cost_inputs <- female_cost_inputs}
  
  complication_cost_matrix <- inner_join(simulation_patients_history[c("CURR.AGE","FEMALE","dead", 
                                                                       "IHD.EVENT", "IHD.HIST", "MI.EVENT", "MI.HIST", "CHF.EVENT", "CHF.HIST",
                                                                       "STROKE.EVENT", "STROKE.HIST", "AMP1.EVENT", "AMP2.EVENT", "AMP.HIST",
                                                                       "BLIND.EVENT","BLIND.HIST", "ULCER.EVENT", "ULCER.HIST")], 
                                         complication_cost_inputs, by = 'CURR.AGE')
  
  #Ischemic heart disease/Angina
  fatal_IHD_cost      <- (complication_cost_matrix$IHD.FATAL*complication_cost_matrix$IHD.EVENT*complication_cost_matrix$dead)/cost_discount_factor 
  nonfatal_IHD_cost   <- (complication_cost_matrix$IHD.NONFATAL*complication_cost_matrix$IHD.EVENT*(1-complication_cost_matrix$dead))/cost_discount_factor 
  subsequent_IHD_cost <- (complication_cost_matrix$IHD.SUB*complication_cost_matrix$IHD.HIST*(1-complication_cost_matrix$IHD.EVENT))/cost_discount_factor 
  simulation_patients_history$IHD.COST <- (fatal_IHD_cost + nonfatal_IHD_cost + subsequent_IHD_cost)
  
  # Myocardial infarction
  fatal_MI_cost      <- (complication_cost_matrix$MI.FATAL*complication_cost_matrix$MI.EVENT*complication_cost_matrix$dead)/cost_discount_factor
  nonfatal_MI_cost   <- (complication_cost_matrix$MI.NONFATAL*complication_cost_matrix$MI.EVENT*(1-complication_cost_matrix$dead))/cost_discount_factor
  subsequent_MI_cost <- (complication_cost_matrix$MI.SUB*complication_cost_matrix$MI.HIST*(1-complication_cost_matrix$MI.EVENT))/cost_discount_factor
  simulation_patients_history$MI.COST <- (fatal_MI_cost + nonfatal_MI_cost + subsequent_MI_cost)
  
  # Heart failure
  fatal_CHF_cost      <- (complication_cost_matrix$CHF.FATAL*complication_cost_matrix$CHF.EVENT*complication_cost_matrix$dead)/cost_discount_factor
  nonfatal_CHF_cost   <- (complication_cost_matrix$CHF.NONFATAL*complication_cost_matrix$CHF.EVENT*(1-complication_cost_matrix$dead))/cost_discount_factor
  subsequent_CHF_cost <- (complication_cost_matrix$CHF.SUB*complication_cost_matrix$CHF.HIST*(1-complication_cost_matrix$CHF.EVENT))/cost_discount_factor
  simulation_patients_history$CHF.COST <- (fatal_CHF_cost + nonfatal_CHF_cost + subsequent_CHF_cost)
  
  #Stroke
  fatal_STROKE_cost      <- (complication_cost_matrix$STROKE.FATAL*complication_cost_matrix$STROKE.EVENT*complication_cost_matrix$dead)/cost_discount_factor
  nonfatal_STROKE_cost   <- (complication_cost_matrix$STROKE.NONFATAL*complication_cost_matrix$STROKE.EVENT*(1-complication_cost_matrix$dead))/cost_discount_factor
  subsequent_STROKE_cost <- (complication_cost_matrix$STROKE.SUB*complication_cost_matrix$STROKE.HIST*(1-complication_cost_matrix$STROKE.EVENT))/cost_discount_factor
  simulation_patients_history$STROKE.COST <- (fatal_STROKE_cost + nonfatal_STROKE_cost + subsequent_STROKE_cost)
  
  # Amputation 	0	15,153	5,328	Alva et al. 2015 [1]
  # Note here we have in the model AMP1.EVENT and AMP2.EVENT but in terms of costs, there is no distinction
  # Question: in case of 2nd amp, should we consider the sum of the event + history of first?
  # Not sure whether costs for amputation are calculated ok. It should be consistent with other events when a 2nd event can occur.
  
  fatal_AMP_cost      <- (complication_cost_matrix$AMP.FATAL*complication_cost_matrix$AMP1.EVENT*complication_cost_matrix$dead + complication_cost_matrix$AMP.FATAL*complication_cost_matrix$AMP2.EVENT*complication_cost_matrix$dead)/cost_discount_factor
  nonfatal_AMP_cost   <- (complication_cost_matrix$AMP.NONFATAL*complication_cost_matrix$AMP1.EVENT*(1-complication_cost_matrix$dead) + complication_cost_matrix$AMP.NONFATAL*complication_cost_matrix$AMP2.EVENT*(1-complication_cost_matrix$dead))/cost_discount_factor
  subsequent_AMP_cost <- (complication_cost_matrix$AMP.SUB*complication_cost_matrix$AMP.HIST*(1-complication_cost_matrix$AMP1.EVENT)*(1-complication_cost_matrix$AMP2.EVENT))/cost_discount_factor
  simulation_patients_history$AMP.COST <- (fatal_AMP_cost + nonfatal_AMP_cost + subsequent_AMP_cost)

  #Blindness
  fatal_BLIND_cost      <- (complication_cost_matrix$BLIND.FATAL*complication_cost_matrix$BLIND.EVENT*complication_cost_matrix$dead)/cost_discount_factor
  nonfatal_BLIND_cost   <- (complication_cost_matrix$BLIND.NONFATAL*complication_cost_matrix$BLIND.EVENT*(1-complication_cost_matrix$dead))/cost_discount_factor
  subsequent_BLIND_cost <- (complication_cost_matrix$BLIND.SUB*complication_cost_matrix$BLIND.HIST*(1-complication_cost_matrix$BLIND.EVENT))/cost_discount_factor
  simulation_patients_history$BLIND.COST <- (fatal_BLIND_cost + nonfatal_BLIND_cost + subsequent_BLIND_cost)
  
  # Ulcer
  fatal_ULCER_cost      <- (complication_cost_matrix$ULCER.FATAL*complication_cost_matrix$ULCER.EVENT*complication_cost_matrix$dead)/cost_discount_factor
  nonfatal_ULCER_cost   <- (complication_cost_matrix$ULCER.NONFATAL*complication_cost_matrix$ULCER.EVENT*(1-complication_cost_matrix$dead))/cost_discount_factor
  subsequent_ULCER_cost <- (complication_cost_matrix$ULCER.SUB*complication_cost_matrix$ULCER.HIST*(1-complication_cost_matrix$ULCER.EVENT))/cost_discount_factor
  simulation_patients_history$ULCER.COST <- (fatal_ULCER_cost + nonfatal_ULCER_cost + subsequent_ULCER_cost)
 
  # Cost in the absence of complications 	1,990	Alva et al. 2015 [1]
  simulation_patients_history$NOCOMP.COST <- (complication_cost_matrix$NOCOMP*(1-complication_cost_matrix$IHD.EVENT)*(1-complication_cost_matrix$MI.EVENT)*
                                              (1-complication_cost_matrix$CHF.EVENT)*(1-complication_cost_matrix$STROKE.EVENT)*(1-complication_cost_matrix$AMP1.EVENT)*
                                                (1-complication_cost_matrix$AMP2.EVENT)*(1-complication_cost_matrix$BLIND.EVENT)*(1-complication_cost_matrix$ULCER.EVENT))/cost_discount_factor
  
  # Intervention costs: only applied the first year of simulation. Then equal to 0.
  simulation_patients_history$TX.COST <- if_else(simulation_patients_history$SDURATION == 0, (tx_cost_input)/cost_discount_factor,0)
  
  ### Societal costs: Added 29/08/2020
  
  # Informal care and productivity loss: We calculated above INF.CARE & PROD.LOSS = 0/1 but costs have to be calculated here: 
  
  # Informal care cost per hour assumed = ?20.26. Make this an input parameter in the model.
  # Calculation: From Weatherly et al. 2014 "Valuing Informal Care for Economic Evaluation". Table 1 estimate from Wilson et al. 2009 = ?13.11 (2004/UK)
  # Inflated to 2020 costs using https://www.officialdata.org/uk/inflation/2004?amount=13.11
  simulation_patients_history$INF.CARE.COST <- (20.26*simulation_patients_history$INF.CARE*(1-simulation_patients_history$dead) + 20.26/2*simulation_patients_history$INF.CARE*simulation_patients_history$dead)/cost_discount_factor
  
  simulation_patients_history$PROD.LOSS.COST <- (simulation_patients_history$PROD.LOSS*(1-simulation_patients_history$dead) + 1/2*simulation_patients_history$PROD.LOSS*simulation_patients_history$dead)/cost_discount_factor
  
  # Future costs are dependent on age, gender and alive status: add explanation later
  future_cost_matrix <- inner_join(simulation_patients_history[c("CURR.AGE","FEMALE","dead")],
                                   future_medical_cost_inputs, by = 'CURR.AGE')
  future_cost_matrix <- inner_join(future_cost_matrix,future_nonmedical_cost_inputs, by = 'CURR.AGE')
  
  simulation_patients_history$FUTURE.COST.MEDICAL <- ((1-future_cost_matrix$FEMALE)*(1-future_cost_matrix$dead)*future_cost_matrix[,6] + 
                                                        (1-future_cost_matrix$FEMALE)*future_cost_matrix$dead*future_cost_matrix[,4] + 
                                                        future_cost_matrix$FEMALE*(1-future_cost_matrix$dead)*future_cost_matrix[,7] + 
                                                        future_cost_matrix$FEMALE*future_cost_matrix$dead*future_cost_matrix[,5])/cost_discount_factor 
  
  simulation_patients_history$FUTURE.COST.NONMEDICAL <- future_cost_matrix[,8]/cost_discount_factor 
  
  # Total annual costs discounted
  simulation_patients_history$TOTAL.COST <- (simulation_patients_history$IHD.COST + 
                                             simulation_patients_history$MI.COST + 
                                             simulation_patients_history$CHF.COST + 
                                             simulation_patients_history$STROKE.COST + 
                                             simulation_patients_history$AMP.COST + 
                                             simulation_patients_history$BLIND.COST + 
                                             simulation_patients_history$ULCER.COST + 
                                             simulation_patients_history$NOCOMP.COST + 
                                             simulation_patients_history$TX.COST + 
                                             simulation_patients_history$INF.CARE.COST + 
                                             simulation_patients_history$PROD.LOSS.COST + 
                                             simulation_patients_history$FUTURE.COST.MEDICAL + 
                                             simulation_patients_history$FUTURE.COST.NONMEDICAL) 
  
  
  ########## MAIN PART III: Calculate QALYs ##########
  
  # Discounting is applied to total QALYs only at this moment.
  qaly_discount_factor <- ((1+qol_disc_rate_input)^(simulation_patients_history$SDURATION)) # Note if SDURATION starts at 0 then we don't need to add -1
  
  # Age and gender dependent
  qol_matrix <- inner_join(simulation_patients_history[c("CURR.AGE","FEMALE","CHF.EVENT","IHD.EVENT","MI.EVENT", "STROKE.EVENT",
                                                         "ULCER.EVENT", "AMP1.EVENT","AMP2.EVENT", "BLIND.EVENT", "RENAL.EVENT","BMI")],qol_inputs, by = 'CURR.AGE')
  
  # The following event utility decrements are at this moment included in the model (sourced from Beaudet at al.)
  
  # CHF.EVENT	IHD.EVENT	MI.EVENT	STROKE.EVENT	ULCER.EVENT	AMP.EVENT	BLIND.EVENT	RENAL.EVENT	BMI.HIGH
    # -0.089	   -0.074	    -0.045	  -0.135	     -0.14	      -0.231	   -0.061	     -0.135	   -0.005
   
  # Coronary heart disease group:	
  qol_matrix$CHF.EVENT <- qol_matrix$CHF.EVENT*qol_events_inputs$CHF
  qol_matrix$IHD.EVENT <- qol_matrix$IHD.EVENT*qol_events_inputs$IHD
  qol_matrix$MI.EVENT  <- qol_matrix$MI.EVENT*qol_events_inputs$MI
  
  chd_qaly_male   <- (1-qol_matrix$FEMALE)*pmin(qol_matrix$CHF.EVENT,qol_matrix$IHD.EVENT,qol_matrix$MI.EVENT)
  chd_qaly_female <- qol_matrix$FEMALE*pmin(qol_matrix$CHF.EVENT,qol_matrix$IHD.EVENT,qol_matrix$MI.EVENT)
  
  simulation_patients_history$CHD.QALY <- (chd_qaly_male + chd_qaly_female)/qaly_discount_factor
  
  # Cerebrovascular disease: only stroke
  qol_matrix$STROKE.EVENT <- qol_matrix$STROKE.EVENT*qol_events_inputs$STROKE
  stroke_qaly_male   <- (1-qol_matrix$FEMALE)*qol_matrix$STROKE.EVENT
  stroke_qaly_female <- qol_matrix$FEMALE*qol_matrix$STROKE.EVENT
  
  simulation_patients_history$STROKE.QALY <- (stroke_qaly_male + stroke_qaly_female)/qaly_discount_factor
  
  # Neuropathy group: so far we have ULCER.EVENT, AMP1.EVENT and AMP2.EVENT
  qol_matrix$ULCER.EVENT <- qol_matrix$ULCER.EVENT*qol_events_inputs$ULCER
  qol_matrix$AMP1.EVENT  <- qol_matrix$AMP1.EVENT*qol_events_inputs$AMP1
  qol_matrix$AMP2.EVENT  <- qol_matrix$AMP2.EVENT*qol_events_inputs$AMP2
  neuro_qaly_male   <- (1-qol_matrix$FEMALE)*pmin(qol_matrix$ULCER.EVENT,qol_matrix$AMP1.EVENT,qol_matrix$AMP2.EVENT)
  neuro_qaly_female <- qol_matrix$FEMALE*pmin(qol_matrix$ULCER.EVENT,qol_matrix$AMP1.EVENT,qol_matrix$AMP2.EVENT)
  simulation_patients_history$NEUROPATHY.QALY <- (neuro_qaly_male + neuro_qaly_female)/qaly_discount_factor
  
  # Retinopathy: so far we only have BLIND.EVENT 
  qol_matrix$BLIND.EVENT <- qol_matrix$BLIND.EVENT*qol_events_inputs$BLIND
  blind_qaly_male   <- (1-qol_matrix$FEMALE)*qol_matrix$BLIND.EVENT
  blind_qaly_female <- qol_matrix$FEMALE*qol_matrix$BLIND.EVENT
  simulation_patients_history$BLIND.QALY <- (blind_qaly_male + blind_qaly_female)/qaly_discount_factor
  
  
  # Nephropathy: so far we only have RENAL.EVENT: assumed disutility from haemodialysis but not sure
  qol_matrix$RENAL.EVENT <- qol_matrix$RENAL.EVENT*qol_events_inputs$RENAL
  renal_qaly_male   <- (1-qol_matrix$FEMALE)*qol_matrix$RENAL.EVENT
  renal_qaly_female <- qol_matrix$FEMALE*qol_matrix$RENAL.EVENT
  simulation_patients_history$RENAL.QALY <- (renal_qaly_male + renal_qaly_female)/qaly_discount_factor
  
  # Comorbidity: Excess BMI (each unit above 25 kg/m2)	-0.006	-0.008	-0.004
  qol_matrix$BMI <- max(qol_matrix$BMI-25,0)*(-0.005)
  bmi_qaly_male   <- (1-qol_matrix$FEMALE)*qol_matrix$BMI
  bmi_qaly_female <- qol_matrix$FEMALE*qol_matrix$BMI
  simulation_patients_history$BMI.QALY <- (bmi_qaly_male + bmi_qaly_female)/qaly_discount_factor
  # Note: 25 in formula above is HARDCODED!
  
  # Total annual utilities
  total_utils <- ((1-qol_matrix$FEMALE)*(qol_matrix$BASELINE.MALE) + qol_matrix$FEMALE*qol_matrix$BASELINE.FEMALE +
                    simulation_patients_history$CHD.QALY + simulation_patients_history$STROKE.QALY + 
                    simulation_patients_history$NEUROPATHY.QALY + simulation_patients_history$BLIND.QALY + 
                    simulation_patients_history$RENAL.QALY + simulation_patients_history$BMI.QALY)
                    
  # Annual discounted QALYs
  simulation_patients_history$QALY <- (total_utils*(1-simulation_patients_history$dead) + (total_utils/2)*simulation_patients_history$dead)
  
  #View(simulation_patients_history)
  
  ########## MAIN PART IV: Calculate aggregated results ##########
  
  # To be completed: add clinical outcomes - life expectancy and events
  
  # Life expectancy
  mean_life_expectancy <- mean(simulation_patients_history[which(simulation_patients_history$dead==1),"SDURATION"])
  
  # Event rates: note events calculated differently depending on how were defined: .EVENT or .HIST
  # Total number of events per patient during lifetime and rate per year
  cum_CHF <- aggregate(simulation_patients_history$CHF.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_CHF_rate <- sum(cum_CHF$x)/patient_size_input
  
  cum_BLIND <- aggregate(simulation_patients_history$BLIND.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_BLIND_rate <- sum(cum_BLIND$x)/patient_size_input
  
  cum_ULCER <- aggregate(simulation_patients_history$ULCER.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_ULCER_rate <- sum(cum_ULCER$x)/patient_size_input
  
  cum_AMP1 <- aggregate(simulation_patients_history$AMP1.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_AMP1_rate <- sum(cum_AMP1$x)/patient_size_input
  
  cum_AMP2 <- aggregate(simulation_patients_history$AMP2.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_AMP2_rate <- sum(cum_AMP2$x)/patient_size_input
  
  cum_MI <- aggregate(simulation_patients_history$MI.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_MI_rate <- sum(cum_MI$x)/patient_size_input
  
  cum_STROKE <- aggregate(simulation_patients_history$STROKE.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_STROKE_rate <- sum(cum_STROKE$x)/patient_size_input
  
  cum_RENAL <- aggregate(simulation_patients_history$RENAL.EVENT, list(Patient = simulation_patients_history$SIMID), sum)
  mean_RENAL_rate <- sum(cum_RENAL$x)/patient_size_input #mean(cum_RENAL$x)/mean_life_expectancy # note this could be per patient or per time
  
  # Total costs per patient during lifetime and average per patient 
  total_costs_patient <- aggregate(simulation_patients_history$TOTAL.COST, list(Patient = simulation_patients_history$SIMID), sum)
  mean_total_costs <- sum(total_costs_patient$x)/patient_size_input
  
  # Breakdown costs: 
  tx_costs_patient <- aggregate(simulation_patients_history$TX.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_tx_costs <- sum(tx_costs_patient$x)/patient_size_input
  
  ihd_costs_patient <- aggregate(simulation_patients_history$IHD.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_ihd_costs <- sum(ihd_costs_patient$x)/patient_size_input
  
  mi_costs_patient <- aggregate(simulation_patients_history$MI.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_mi_costs <- sum(mi_costs_patient$x)/patient_size_input
  
  chf_costs_patient <- aggregate(simulation_patients_history$CHF.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_chf_costs <- sum(chf_costs_patient$x)/patient_size_input
  
  stroke_costs_patient <- aggregate(simulation_patients_history$STROKE.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_stroke_costs <- sum(stroke_costs_patient$x)/patient_size_input
  
  amp_costs_patient <- aggregate(simulation_patients_history$AMP.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_amp_costs <- sum(amp_costs_patient$x)/patient_size_input
  
  blind_costs_patient <- aggregate(simulation_patients_history$BLIND.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_blind_costs <- sum(blind_costs_patient$x)/patient_size_input
  
  ulcer_costs_patient <- aggregate(simulation_patients_history$ULCER.COST , list(Patient = simulation_patients_history$SIMID), sum)
  mean_ulcer_costs <- sum(ulcer_costs_patient$x)/patient_size_input
  
  mean_complication_costs <- mean_ihd_costs + mean_mi_costs + mean_chf_costs + mean_stroke_costs + mean_amp_costs + mean_blind_costs + mean_ulcer_costs
  
  # No complication costs
  
  nocomp_costs_patient <- aggregate(simulation_patients_history$NOCOMP.COST  , list(Patient = simulation_patients_history$SIMID), sum)
  mean_nocomp_costs <- sum(nocomp_costs_patient$x)/patient_size_input
  
  # Informal care and productivity
  inf_care_costs_patient <- aggregate(simulation_patients_history$INF.CARE.COST, list(Patient = simulation_patients_history$SIMID), sum)
  mean_inf_care_costs <- sum(inf_care_costs_patient$x)/patient_size_input
  
  prod_loss_costs_patient <- aggregate(simulation_patients_history$PROD.LOSS.COST, list(Patient = simulation_patients_history$SIMID), sum)
  mean_prod_loss_costs <- sum(prod_loss_costs_patient$x)/patient_size_input
  
  future_medical_costs_patient <- aggregate(simulation_patients_history$FUTURE.COST.MEDICAL, list(Patient = simulation_patients_history$SIMID), sum)
  mean_future_medical_costs <- sum(future_medical_costs_patient$x)/patient_size_input
  
  future_nonmedical_costs_patient <- aggregate(simulation_patients_history$FUTURE.COST.NONMEDICAL, list(Patient = simulation_patients_history$SIMID), sum)
  mean_future_nonmedical_costs <- sum(future_nonmedical_costs_patient$x)/patient_size_input
  
  # QALYs
  total_qalys_patient <- aggregate(simulation_patients_history$QALY, list(Patient = simulation_patients_history$SIMID), sum)
  mean_total_qalys <- sum(total_qalys_patient$x)/patient_size_input
  
  
  ### Return model outcomes: 
  return(list(simulation_patients_history=simulation_patients_history, 
              mean_life_expectancy = mean_life_expectancy,
              mean_CHF_rate = mean_CHF_rate, mean_BLIND_rate = mean_BLIND_rate, mean_ULCER_rate = mean_ULCER_rate, mean_STROKE_rate = mean_STROKE_rate,
              mean_AMP1_rate = mean_AMP1_rate, mean_AMP2_rate = mean_AMP2_rate, mean_MI_rate = mean_MI_rate, mean_RENAL_rate = mean_RENAL_rate,
              mean_total_costs = mean_total_costs, 
              mean_total_qalys = mean_total_qalys,
              mean_inf_care_costs = mean_inf_care_costs, 
              mean_prod_loss_costs = mean_prod_loss_costs,
              mean_ihd_costs = mean_ihd_costs, 
              mean_mi_costs = mean_mi_costs, 
              mean_chf_costs = mean_chf_costs,
              mean_stroke_costs = mean_stroke_costs, 
              mean_amp_costs = mean_amp_costs, 
              mean_blind_costs = mean_blind_costs,
              mean_ulcer_costs = mean_ulcer_costs, 
              mean_complication_costs = mean_complication_costs,
              mean_nocomp_costs = mean_nocomp_costs, 
              mean_tx_costs = mean_tx_costs,
              mean_future_medical_costs = mean_future_medical_costs, 
              mean_future_nonmedical_costs = mean_future_nonmedical_costs)
         ) # end return parameters
  
} #end SMDMII_model_simulation function
