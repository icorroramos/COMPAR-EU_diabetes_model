# Clear all objects and clear out memory
rm(list = ls())
gc()
#Sys.sleep(90)
memory.size (max = F)

# ANALYSIS INPUT VARIABLES ------------------------------------------------
country.id <- 'UK' # Choose from: 'UK', 'NL', 'DE', 'ES', 'GR'


# Load model function: 
source("R/SMI in type II diabetes - HE model v3.R")

# Load aux. functions, input parameters (from Excel), global lists, etc. 
source("R/aux_functions.R")

# Load saved vector with 1000 RNG seed values
load('R/Psa_seed.Rdata')

# Please select running mode: 0 = deterministic, 1 = PSA. Default: 0
psa_input <- 0

# Please select number of PSA runs (only works if psa_input <- 1, otherwise will be ignored ). Default: 500 --> NOTE: not completely implemented
n_psa_input <- 1000

# Please select the number of patients in simulation (default 1000 in deterministic run)
npats_input  <- 150

# Please indicate the name of the treatment to be identified
tx_label <- "Usual care"

# Cost inputs
tx_cost_input <- 0 # Total treatment cost --> Not sure if here or in Excel.     
retirement_age_input <- 66

# Print time to be able to estimate finishing time
print(Sys.time())


##################
# Set-up options #
##################

# Variable defined to keep track of simulation time (delete afterwards)
init <- Sys.time()

# Get treatment effec draws
draws.hba1c <- rpert(n_psa_input, min = -1.4638, mode = -0.8232, max = -0.1825)
draws.bmi <- rpert(n_psa_input, min = -3.0261, mode = -1.7, max = -0.3739)
draws.sbp <- rpert(n_psa_input, min = -18.7532, mode = -10.03, max = -1.3068)


# Create empty objects to store PSA data in
sim.var.data <- matrix(nrow = n_psa_input, ncol = 11)
female.int.data <- list()
female.comp.data <- list()
male.int.data <- list()
male.comp.data <- list()

for (i in 1:n_psa_input){
  # Treatment effect inputs
  treateff_start   <- 1 # Cycle in which treatment effect starts
  treateff_end     <- 3 # Cycle in which treatment effect ends
  treateff_decline <- 2 # Cycle in which treatment effect starts to decline (linearly)
  
  treateff_hba1c <- draws.hba1c[i] # Treatment effect on HbA1c (in absolute %-points HbA1c)
  treateff_hdl   <- 0 * 0.02586 * 10 # TRANSFORMATION *10 FOR MODEL INPUT # Treatment effect on HDL-cholesterol (absolute effect in mmol/l)
  treateff_ldl   <- 0 * 0.02586 * 10 # TRANSFORMATION *10 FOR MODEL INPUT# Treatment effect on LDL-cholesterol (absolute effect in mmol/l)
  treateff_bmi <- draws.bmi[i] # Treatment effect on BMI (in absolute points)
  treateff_sbp <- draws.sbp[i] / 10 # TRANSFORMATION /10 FOR MODEL INPUT# Treatment effect on SBP (in absolute mmHg)
  
  
  # Tx effects are vectors: the current assumption is that the same start, end and decline is assumed for all effect modifiers
  treatment_effect_HbA1c_input <- c(treateff_hba1c, treateff_start, treateff_end, treateff_decline)
  treatment_effect_HDL_input <- c(treateff_hdl, treateff_start, treateff_end, treateff_decline)
  treatment_effect_LDL_input <- c(treateff_ldl, treateff_start, treateff_end, treateff_decline)
  treatment_effect_BMI_input <- c(treateff_bmi, treateff_start, treateff_end, treateff_decline)
  treatment_effect_SBP_input <- c(treateff_sbp, treateff_start, treateff_end, treateff_decline)
  
  
  # RESULT OUTPUT FOR OUTPUT SCRIPT -----------------------------------------
  sim.results.female <- SMDMII_model_simulation(patient_size_input = npats_input,
                                                female_input = 1,
                                                tx_cost_input = tx_cost_input,
                                                treatment_effect_HbA1c_input = treatment_effect_HbA1c_input, 
                                                treatment_effect_HDL_input = treatment_effect_HDL_input,
                                                treatment_effect_LDL_input = treatment_effect_LDL_input, 
                                                treatment_effect_BMI_input = treatment_effect_BMI_input,
                                                treatment_effect_SBP_input = treatment_effect_SBP_input,
                                                cost_disc_rate_input = discount_factors[1], 
                                                qol_disc_rate_input = discount_factors[2],
                                                retirement_age_input = retirement_age_input, 
                                                inf_care_hours_input = inf_care_hours_input, 
                                                worked_hours_input = worked_hours_input,
                                                working_days_lost_input = working_days_lost_input,
                                                cost_hour_sick_input = cost_hour_sick_input,
                                                friction_period_input = friction_period_input, 
                                                informal_care_coef_input = informal_care_coef_input,
                                                employment_coef_input = employment_coef_input,
                                                prod_costs_coef_input = prod_costs_coef_input, 
                                                inf_care_age_scale_input = inf_care_age_scale_input, 
                                                prod_loss_age_scale_input = prod_loss_age_scale_input,
                                                run_PSA_input = psa_input,
                                                seed_input = psa.seed[i])
  
  sim.results.male <- SMDMII_model_simulation(patient_size_input = npats_input,
                                              female_input = 0,
                                              tx_cost_input = tx_cost_input,
                                              treatment_effect_HbA1c_input = treatment_effect_HbA1c_input, 
                                              treatment_effect_HDL_input = treatment_effect_HDL_input,
                                              treatment_effect_LDL_input = treatment_effect_LDL_input, 
                                              treatment_effect_BMI_input = treatment_effect_BMI_input,
                                              treatment_effect_SBP_input = treatment_effect_SBP_input,
                                              cost_disc_rate_input = discount_factors[1], 
                                              qol_disc_rate_input = discount_factors[2],
                                              retirement_age_input = retirement_age_input, 
                                              inf_care_hours_input = inf_care_hours_input, 
                                              worked_hours_input = worked_hours_input,
                                              working_days_lost_input = working_days_lost_input,
                                              cost_hour_sick_input = cost_hour_sick_input,
                                              friction_period_input = friction_period_input, 
                                              informal_care_coef_input = informal_care_coef_input,
                                              employment_coef_input = employment_coef_input,
                                              prod_costs_coef_input = prod_costs_coef_input, 
                                              inf_care_age_scale_input = inf_care_age_scale_input, 
                                              prod_loss_age_scale_input = prod_loss_age_scale_input,
                                              run_PSA_input = psa_input,
                                              seed_input = psa.seed[i])
  
  
  # Store iteration results in data objects
  sim.var.data[i, ] <- c(psa.seed[i], npats_input, tx_cost_input, unlist(mget(apropos('treateff.'))))
  female.int.data[[i]] <- unlist(sim.results.female[-1])
  male.int.data[[i]] <- unlist(sim.results.male[-1])
}

colnames(sim.var.data) <- c('Seed', 'Npats', 'Tx_Cost', 'BMI', 'Decline', 'End', 'HbA1c', 'HDL', 'LDL', 'SBP', 'Start')
psa.results.female <- as.data.frame(do.call(rbind, female.int.data))
psa.results.male <- as.data.frame(do.call(rbind, male.int.data))

# Print simulation duration
end <- Sys.time()
print(end - init)


# Save simulation results
save(sim.var.data,
     psa.results.female,
     psa.results.male,
     file = 'output/PSA_Rank1.RData')