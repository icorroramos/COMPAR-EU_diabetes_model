##################
# Set-up options #
##################

# Check and install the latest R version
# installr::updateR()
# update.packages() 

# Setting option for decimals
options(scipen = 3)

# Install and load required packages
# Sys.getenv() # check this for the Home path where packages are saved
.libPaths(Sys.getenv()[["R_LIBS_USER"]])

pkgs <- c("lattice", "MASS", "plyr", "survival", "tidyverse", "dplyr", "beepr") # package names
#install.packages(pkgs)
inst <- lapply(pkgs, library, character.only = TRUE) # load them

# ANALYSIS INPUT VARIABLES ------------------------------------------------
country.id <- 'UK' # Choose from: 'UK', 'NL', 'DE', 'ES', 'GR'

# Please select running mode: 0 = deterministic, 1 = PSA. Default: 0
psa_input <- 0

# Please select number of PSA runs (only works if psa_input <- 1, otherwise will be ignored ). Default: 500 --> NOTE: not completely implemented
n_psa_input <- 5

# Please select the number of patients in simulation (default 1000 in deterministic run)
npats_input   <- 1000 

# Please indicate the name of the treatment to be identified
tx_label <- "Usual care"

# Cost inputs
tx_cost_input <- 0 # Total treatment cost --> Not sure if here or in Excel.     
retirement_age_input <- 66

# Set random seed for replication purposes
seed_input <- 958 # A random seed that it is used to ensure consistency in the model results.

# LOAD MODEL FUNCTION AND AUX FUNCTIONS/ DATA -----------------------------

# Load model function: 
source("R/SMI in type II diabetes - HE model v3.R")
#source("R/OLD SMI in type II diabetes - HE model v3.R")

# Load aux. functions, input parameters (from Excel), global lists, etc. 
source("R/aux_functions.R")

init <- Sys.time()

# Print time to be able to estimate finishing time
print(Sys.time())

# Treatment effect inputs
treateff_start   <- 1 # Cycle in which treatment effect starts
treateff_end     <- 3 # Cycle in which treatment effect ends
treateff_decline <- 2 # Cycle in which treatment effect starts to decline (linearly)

# Treatment effects
treateff_hba1c <- -0.8232 # Treatment effect on HbA1c (in absolute %-points HbA1c)
treateff_hdl   <-  0 * 0.02586 * 10 # TRANSFORMATION *10 FOR MODEL INPUT # Treatment effect on HDL-cholesterol (absolute effect in mmol/l)
treateff_ldl   <-  0 * 0.02586 * 10 # TRANSFORMATION *10 FOR MODEL INPUT# Treatment effect on LDL-cholesterol (absolute effect in mmol/l)
treateff_bmi   <- -1.7 # Treatment effect on BMI (in absolute points)
treateff_sbp   <- -10.03 / 10 # TRANSFORMATION /10 FOR MODEL INPUT# Treatment effect on SBP (in absolute mmHg)

# Tx effects are vectors: the current assumption is that the same start, end and decline is assumed for all effect modifiers
treatment_effect_HbA1c_input <- c(treateff_hba1c, treateff_start, treateff_end, treateff_decline)
treatment_effect_HDL_input   <- c(treateff_hdl, treateff_start, treateff_end, treateff_decline)
treatment_effect_LDL_input   <- c(treateff_ldl, treateff_start, treateff_end, treateff_decline)
treatment_effect_BMI_input   <- c(treateff_bmi, treateff_start, treateff_end, treateff_decline)
treatment_effect_SBP_input   <- c(treateff_sbp, treateff_start, treateff_end, treateff_decline)

# NOTE: Input parameters (probabilities, costs and utilities) can be changed in the corresponding 
# csv files included in the folder "input". These can be changed to run for example scenario analyses
# without modifying the R code.

# # Directories to save results and plots (TEST folder created for users)
# results_dir <- ("output/TEST/")
# dir.create(results_dir)
# 
# graphics_dir <- ("graphics/TEST/") # not used at this moment
# dir.create(graphics_dir)

# The sim.vars object collects all parameters that define the simulation into one object
# Then, it is saved with the output of the simulation. This way, if we  have multiple output files, we always have the information 
# on the relevant input parameters that were used to produce the output.
sim.vars <- list(seed_input, npats_input, tx_cost_input, mget(apropos('treateff.')))

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
                                              seed_input = seed_input)

sim.results.female.comp <- SMDMII_model_simulation(patient_size_input = npats_input,
                                                   female_input = 1,
                                                   tx_cost_input = tx_cost_input,
                                                   treatment_effect_HbA1c_input = rep(0,4),
                                                   treatment_effect_HDL_input = rep(0,4),
                                                   treatment_effect_LDL_input = rep(0,4),
                                                   treatment_effect_BMI_input = rep(0,4),
                                                   treatment_effect_SBP_input = rep(0,4),
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
                                                   seed_input = seed_input)

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
                                            seed_input = seed_input)

sim.results.male.comp <- SMDMII_model_simulation(patient_size_input = npats_input,
                                                 female_input = 0,
                                                 tx_cost_input = tx_cost_input,
                                                 treatment_effect_HbA1c_input = rep(0,4),
                                                 treatment_effect_HDL_input = rep(0,4),
                                                 treatment_effect_LDL_input = rep(0,4),
                                                 treatment_effect_BMI_input = rep(0,4),
                                                 treatment_effect_SBP_input = rep(0,4),
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
                                                 seed_input = seed_input)

# Print simulation duration
end <- Sys.time()
print(end - init)

# Save simulation results
save(sim.vars,
     sim.results.female,
     sim.results.female.comp,
     sim.results.male,
     sim.results.male.comp,
     file = 'output/techver/TEST_RUN_Test10.RData') #file = paste0('output/', comp, '_seed_', seed_input, '.RData')