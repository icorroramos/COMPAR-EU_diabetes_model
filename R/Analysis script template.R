# Model function: 
source("R/SMI in type II diabetes - HE model v3.R")

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
 
# Please select number of PSA runs. Default: 500 --> NOTE: not yet implemented
n_psa_input <- 500

# Please select the number of patients in simulation (default 1000 in deterministic run)
npats_input   <- 10  

# Treatment effect inputs
treateff_start_input   <- 1 # Cycle in which treatment effect starts
treateff_end_input     <- 4 # Cycle in which treatment effect ends
treateff_decline_input <- 2 # Cycle in which treatment effect starts to decline (linearly)

treateff_hba1c_input <- 0 # Treatment effect on HbA1c (in absolute %-points HbA1c)
treateff_hdl_input   <- 0 # Treatment effect on HDL-cholesterol (absolute effect, which unit??)
treateff_ldl_input   <- 0 # Treatment effect on LDL-cholesterol (absolute effect, which unit??)

# Cost inputs

# Total treatment cost --> Not sure if here or in Excel.
tx_cost_input <- 0     


# Quality of life inputs


# Gender: note at this moment the model distinguishes between males and females, This must be chosen here
female_input <- 0 #1 = female, 0 = male

# NOTE: Input parameters (probabilities, costs and utilities) can be changed in the corresponding 
# csv files included in the folder "input". These can be changed to run for example scenario analyses
# without modifying the R code.
 
# # Directories to save results and plots (TEST folder created for users)
results_dir <- ("output/TEST/")
dir.create(results_dir)

graphics_dir <- ("graphics/TEST/")
dir.create(graphics_dir)

# The sim.vars object collects all parameters that define the simulation into one object
# Then, it is saved with the output of the simulation. This way, if we  have multiple output files, we always have the information 
# on the relevant input parameters that were used to produce the output.
sim.vars <- list(npats_input, tx_cost_input, mget(apropos('treateff.')))

# <----------





###########################################################



###########################################################






# Direct costs of diabetes-related complications for UK are age-gender dependent.
male_cost_inputs   <- read.csv("input/UK/Event cost male 2020.csv", sep=",")
female_cost_inputs <- read.csv("input/UK/Event cost female 2020.csv", sep=",")

# When a societal perspective is adopted, we also have future costs. This are obtained from the PAID online tool.
future_medical_cost_inputs    <- read.csv("input/UK/UKPAID__AllORUnrelated_Costs_2020.csv", sep=",") # Ingelin: 15/12/2020
future_nonmedical_cost_inputs <- read.csv("input/UK/UK_nonmedical_futurecosts_data.csv", sep=",") # UK costs updated 11/12/2020 - corrected version from Hamraz

# Utilities are UK-based and age/gender dependent
qol_inputs <- read.csv("input/UK/qol_inputs_UK.csv", sep=",")

# Utility decrements associated to diabetes-related events
qol_events_inputs <- read.csv("input/UK/qol_events_inputs_UK.csv", sep=",")


# TODO: also extract output processing scripts to here.

# Results 

sim_results <- SMDMII_model_simulation(npats_input,  #patient_size_input: run 500 for LOLA
                                       female_input,  #female_input, 1 = female
                                       tx_cost_input, #tx_cost_input --> Gimon
                                       rep(0,4), #treatment_effect_HbA1c_input --> from COMPAR + Assumption
                                       rep(0,4),  #treatment_effect_HDL_input
                                       rep(0,4),  #treatment_effect_LDL_input --> from COMPAR + Assumption
                                       0, #treatment_effect_BMI_input from MH2020
                                       discount_cost_input, 
                                       discount_util_input, 
                                       65, #retirement_age_input
                                       psa_input, 
                                       77 #seed_input
                                       )

# Results tables
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
rownames(sim_CE_results_table) <- "Comparator"
# View(sim_CE_results_table)

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

rownames(sim_clinical_results_table) <- "Name" # define above auto?

# View(sim_clinical_results_table)


# KM data (to be compraed with UKPDS - validation):
n_years <- 1:max(sim_results$simulation_patients_history$SDURATION)
current_survival <- rep(1,length(event_vars)+1)
KM_data <- sim_results$simulation_patients_history[,c(event_vars,"dead")]
KM_data <- KM_data[FALSE,]
for(i in n_years){
current_survival <- current_survival*(1-colMeans(sim_results$simulation_patients_history[which(sim_results$simulation_patients_history$SDURATION == i),c(event_vars,"dead")]))
KM_data[i,] <- current_survival
}

# View(tail(KM_data))


# Export tables into Excel using this function:
export_csv <- function(object_input){
  write.csv(object_input, paste(results_dir, substitute(object_input),"_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".csv", sep = ""))
}


clinical_results_table <- sim_clinical_results_table # Validation, delete afterwards
  
export_csv(KM_data)
export_csv(clinical_results_table)


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

# Optional - save to csv's
# write.csv(sim_CE_results_female_table, 'sim_CE_female_int.csv', quote = FALSE)
# write.csv(sim_clinical_results_female_table, 'sim_clin_female_int.csv', quote = FALSE)
# write.csv(sim_CE_results_table, 'sim_CE_female_comp.csv', quote = FALSE)
# write.csv(sim_clinical_results_table, 'sim_clin_female_comp.csv', quote = FALSE)
# 
# write.csv(sim_CE_results_male_table,'sim_CE_male_int.csv', quote = FALSE)
# write.csv(sim_clinical_results_male_table, 'sim_clin_male_int.csv', quote = FALSE)
# write.csv(sim_CE_results_male_table_comp, 'sim_CE_male_comp.csv', quote = FALSE)
# write.csv(sim_clinical_results_male_table_comp, 'sim_clin_male_comp.csv', quote = FALSE)