##################
# Set-up options #
##################

# Check and install the latest R version
# installr::updateR()
# update.packages() 

# Remove all objects from the simulation
rm(list = ls())

# Setting option for decimals
options(scipen = 3)

# Install and load required packages
# Sys.getenv() # check this for the Home path where packages are saved
.libPaths(Sys.getenv()[["R_LIBS_USER"]])

pkgs <- c("lattice", "MASS", "plyr", "survival", "tidyverse", "dplyr") # package names
#install.packages(pkgs)
inst <- lapply(pkgs, library, character.only = TRUE) # load them

# Load model function: 
source("R/SMI in type II diabetes - HE model v3.R")

# Load aux. functions, input parameters (from Excel), global lists, etc. 
source("R/aux_functions.R")

# Variable defined to keep track of simulation time (delete afterwards)
init <- Sys.time()

############################
# User adjustable settings #
############################

# Discount rates: please indicate the desired discount rates for costs and effects. Default: 0.035
discount_cost_input <- 0.035
discount_util_input <- 0.035	

# Please select running mode: 0 = deterministic, 1 = PSA. Default: 0
psa_input <- 0
 
# Please select number of PSA runs (only works if psa_input <- 1, otherwise will be ignored ). Default: 500 --> NOTE: not completely implemented
n_psa_input <- 5

# Please select the number of patients in simulation (default 1000 in deterministic run)
npats_input   <- 1000  

# Please indicate the name of the treatment to be identified
tx_label <- "Rank_1"

# Treatment effect inputs
treateff_start   <- 1 # Cycle in which treatment effect starts
treateff_end     <- 4 # Cycle in which treatment effect ends
treateff_decline <- 2 # Cycle in which treatment effect starts to decline (linearly)

treateff_hba1c <- -0.8232 # Treatment effect on HbA1c (in absolute %-points HbA1c)
treateff_hdl   <- 0 # Treatment effect on HDL-cholesterol (absolute effect, which unit??)
treateff_ldl   <- 0 # Treatment effect on LDL-cholesterol (absolute effect, which unit??)

# Tx effects are vectors: the current assumption is that the same start, end and decline is assumed for all effect modifiers
treatment_effect_HbA1c_input <- c(treateff_hba1c, treateff_start, treateff_end, treateff_decline)
treatment_effect_HDL_input   <- c(treateff_hdl, treateff_start, treateff_end, treateff_decline)
treatment_effect_LDL_input   <- c(treateff_ldl, treateff_start, treateff_end, treateff_decline)

treatment_effect_BMI_input <- 0 # This was taken from MH2020 but it si currently removed from the model function. Left here in case we want to include it again

# Cost inputs
tx_cost_input <- 0 # Total treatment cost --> Not sure if here or in Excel.     
retirement_age_input <- 65


# Quality of life inputs


# Gender: note at this moment the model distinguishes between males and females, This must be chosen here
female_input <- 1 #1 = female, 0 = male

# Set random seed for replication purposes
seed_input <- 77 # A random seed that it is used to ensure consistency in the model results. 

# NOTE: Input parameters (probabilities, costs and utilities) can be changed in the corresponding 
# csv files included in the folder "input". These can be changed to run for example scenario analyses
# without modifying the R code.
 
# # Directories to save results and plots (TEST folder created for users)
results_dir <- ("output/TEST/")
dir.create(results_dir)

graphics_dir <- ("graphics/TEST/") # not used at this moment
dir.create(graphics_dir)

# The sim.vars object collects all parameters that define the simulation into one object
# Then, it is saved with the output of the simulation. This way, if we  have multiple output files, we always have the information 
# on the relevant input parameters that were used to produce the output.
sim.vars <- list(npats_input, tx_cost_input, mget(apropos('treateff.')))

###########
# Results #
###########

# Get model results.
# Define sim_results as a function. This allows calling it multiple times for the PSA.
# We need to indicate which parameters change in the PSA by adding e.g. [x] to the input parameters.
# We need first to have a vector of inputs defined in advance (probably in aux file).

sim_results <- function(x){
  SMDMII_model_simulation(npats_input,
                          female_input,
                          tx_cost_input,
                          treatment_effect_HbA1c_input, 
                          treatment_effect_HDL_input,
                          treatment_effect_LDL_input, 
                          treatment_effect_BMI_input,
                          discount_cost_input, 
                          discount_util_input, 
                          retirement_age_input, 
                          psa_input, #might not be needed inside the function
                          seed_input)
}

# Collect results depending on whether the simulation is deterministic or probabilistic 

if(psa_input == 0){
  
  # Run simulation
  sim_results <- sim_results()
  
  # Results tables: cost effectiveness and clinical results
  sim_CE_results_table <- matrix(c(sim_results$mean_complication_costs,
                                   sim_results$mean_nocomp_costs,
                                   sim_results$mean_tx_costs,
                                   sim_results$mean_inf_care_costs,
                                   sim_results$mean_prod_loss_costs,
                                   sim_results$mean_future_medical_costs,
                                   sim_results$mean_future_nonmedical_costs,
                                   sim_results$mean_total_costs,
                                   sim_results$mean_total_qalys), nrow = 1)
  
  colnames(sim_CE_results_table) <- c("Complication costs", "No complication costs", 
                                      "Tx costs","Informal care costs", "Productivity costs",
                                      "Future medical costs", "Future non-medical costs", "Total costs", "Total QALYs")
  rownames(sim_CE_results_table) <- tx_label
  
  
  sim_clinical_results_table <- matrix(c(sim_results$mean_life_expectancy,
                                         sim_results$mean_CHF_rate,
                                         sim_results$mean_MI_rate,
                                         sim_results$mean_BLIND_rate,
                                         sim_results$mean_ULCER_rate,
                                         sim_results$mean_AMP1_rate,
                                         sim_results$mean_AMP2_rate,
                                         sim_results$mean_RENAL_rate,
                                         sim_results$mean_STROKE_rate), nrow = 1)
  
  colnames(sim_clinical_results_table) <- c("Life expectancy", "CHF rate", "MI rate", "Blindness rate", "Ulcer rate",
                                            "1st amputation rate", "2nd amputation rate", "Renal failure rate", "Stroke rate")
  rownames(sim_clinical_results_table) <- tx_label
  
  
  # KM curves from simulated data (to be compared with UKPDS - validation):
  n_years <- 1:max(sim_results$simulation_patients_history$SDURATION)
  current_survival <- rep(1,length(event_vars)+1)
  KM_data <- sim_results$simulation_patients_history[,c(event_vars,"dead")]
  KM_data <- KM_data[FALSE,]
  for(i in n_years){
    current_survival <- current_survival*(1-colMeans(sim_results$simulation_patients_history[which(sim_results$simulation_patients_history$SDURATION == i),c(event_vars,"dead")]))
    KM_data[i,] <- current_survival
  }
  
  
  # Decide what to save
  export_csv(KM_data)
  export_csv(sim_clinical_results_table)
  export_csv(sim_CE_results_table)
  
  
  }else{ # PSA functionality (loops) seems to work - a bit slow in any case. Need to implement with proper inputs
    
    # run simulation
    sim_results <- lapply(1:n_psa_input, sim_results)
    
    
    # Results tables: cost effectiveness and clinical results
    
    complication_costs      <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_complication_costs))
    nocomp_costs            <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_nocomp_costs))
    tx_costs                <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_tx_costs))
    inf_care_costs          <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_inf_care_costs))
    prod_loss_costs         <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_prod_loss_costs))
    future_medical_costs    <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_future_medical_costs))
    future_nonmedical_costs <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_future_nonmedical_costs))
    total_costs             <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_total_costs))
    total_qalys             <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_total_qalys))
    
    # This is a summary but we also need results per PSA for the plots
    sim_CE_results_table <- matrix(c(mean(complication_costs), 
                                     mean(nocomp_costs),
                                     mean(tx_costs),
                                     mean(inf_care_costs),
                                     mean(prod_loss_costs),
                                     mean(future_medical_costs),
                                     mean(future_nonmedical_costs),
                                     mean(total_costs),
                                     mean(total_qalys)), nrow = 1)
    
    colnames(sim_CE_results_table) <- c("Complication costs", "No complication costs", 
                                        "Tx costs","Informal care costs", "Productivity costs",
                                        "Future medical costs", "Future non-medical costs", "Total costs", "Total QALYs")
    rownames(sim_CE_results_table) <- tx_label
    
    
    # Clinical results 
    life_expectancy <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_life_expectancy))
    CHF_rate        <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_CHF_rate))
    MI_rate         <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_MI_rate))
    BLIND_rate      <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_BLIND_rate))
    ULCER_rate      <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_ULCER_rate))
    AMP1_rate       <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_AMP1_rate))
    AMP2_rate       <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_AMP2_rate))
    RENAL_rate      <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_RENAL_rate))
    STROKE_rate     <- unlist(lapply(1:n_psa_input, function(z) sim_results[[z]]$mean_STROKE_rate))
    
    sim_clinical_results_table <- matrix(c(mean(life_expectancy),
                                           mean(CHF_rate),
                                           mean(MI_rate),
                                           mean(BLIND_rate),
                                           mean(ULCER_rate),
                                           mean(AMP1_rate),
                                           mean(AMP2_rate),
                                           mean(RENAL_rate),
                                           mean(STROKE_rate)), nrow = 1)

    colnames(sim_clinical_results_table) <- c("Life expectancy", "CHF rate", "MI rate", "Blindness rate", "Ulcer rate",
                                              "1st amputation rate", "2nd amputation rate", "Renal failure rate", "Stroke rate")
    rownames(sim_clinical_results_table) <- tx_label
    
    
    # Decide what to save
    export_csv(sim_clinical_results_table)
    export_csv(sim_CE_results_table)
    export_csv(total_costs)
    export_csv(total_qalys)
    
    # CE-plane and CEAC not produced here because this is a single arm simulation. 
    # PSA results will be saved and plots can be done in another script.
    
    }


# Variable defined to keep track of simulation time (delete afterwards)
end <- Sys.time()
end - init


# SINK RESULTS DATA -------------------------------------------------------

# Probably we need to explain the purpose of this Rdata file:
# Change the name automatically so that it is not overwritten?
# How to access it after it is saved?
# Anything else?

# save(sim.vars,
#      sim_CE_results_female_table,
#      sim_clinical_results_female_table,
#      sim_CE_results_table,
#      sim_clinical_results_table,
#      sim_CE_results_male_table,
#      sim_clinical_results_male_table,
#      sim_CE_results_male_table_comp,
#      sim_clinical_results_male_table_comp,
#      file = 'Sim_results.Rdata')

# Optional - save to csv's --> Probably delete - see above export_csv
# write.csv(sim_CE_results_female_table, 'sim_CE_female_int.csv', quote = FALSE)
# write.csv(sim_clinical_results_female_table, 'sim_clin_female_int.csv', quote = FALSE)
# write.csv(sim_CE_results_table, 'sim_CE_female_comp.csv', quote = FALSE)
# write.csv(sim_clinical_results_table, 'sim_clin_female_comp.csv', quote = FALSE)
