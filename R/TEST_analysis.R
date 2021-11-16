##################
# Set-up options #
##################

# Check and install the latest R version
# installr::updateR()
# update.packages() 


# # Setting option for decimals
# options(scipen = 3)
# 
# # Install and load required packages
# # Sys.getenv() # check this for the Home path where packages are saved
# .libPaths(Sys.getenv()[["R_LIBS_USER"]])
# 
# pkgs <- c("lattice", "MASS", "plyr", "survival", "tidyverse", "dplyr", "beepr") # package names
# #install.packages(pkgs)
# inst <- lapply(pkgs, library, character.only = TRUE) # load them

# Load model function: 
source("SMI in type II diabetes - HE model v3.R")

# Load aux. functions, input parameters (from Excel), global lists, etc. 
source("aux_functions.R")


# ANALYSIS INPUT VARIABLES ------------------------------------------------
country.id <- 'DE' # Choose from: 'UK', 'NL', 'DE', 'ES', 'GR'

# Discount rates: please indicate the desired discount rates for costs and effects. Default: 0.035
discount_cost_input <- 0.035
discount_util_input <- 0.035	

# Please select running mode: 0 = deterministic, 1 = PSA. Default: 0
psa_input <- 0

# Please select number of PSA runs (only works if psa_input <- 1, otherwise will be ignored ). Default: 500 --> NOTE: not completely implemented
n_psa_input <- 5

# Please select the number of patients in simulation (default 1000 in deterministic run)
npats_input   <- 5000 # 1000  

# Please indicate the name of the treatment to be identified
tx_label <- "Usual care"

# Cost inputs
tx_cost_input <- 0 # Total treatment cost --> Not sure if here or in Excel.     
retirement_age_input <- 66

# Set random seed for replication purposes
seed_input <- 77 # A random seed that it is used to ensure consistency in the model results.


init <- Sys.time()

# Print time to be able to estimate finishing time
print(Sys.time())

# Treatment effect inputs
treateff_start   <- 1 # Cycle in which treatment effect starts
treateff_end     <- 3 # Cycle in which treatment effect ends
treateff_decline <- 2 # Cycle in which treatment effect starts to decline (linearly)

treateff_hba1c <- 0
treateff_hdl   <- 0 
treateff_ldl   <- 0 
treateff_bmi   <- 0
treateff_sbp   <- 0

# Tx effects are vectors: the current assumption is that the same start, end and decline is assumed for all effect modifiers
treatment_effect_HbA1c_input <- c(treateff_hba1c, treateff_start, treateff_end, treateff_decline)
treatment_effect_HDL_input   <- c(treateff_hdl, treateff_start, treateff_end, treateff_decline)
treatment_effect_LDL_input   <- c(treateff_ldl, treateff_start, treateff_end, treateff_decline)
treatment_effect_BMI_input   <- c(treateff_bmi, treateff_start, treateff_end, treateff_decline)
treatment_effect_SBP_input   <- c(treateff_sbp, treateff_start, treateff_end, treateff_decline)


# Country-specific parameters
inf_care_hours_input      <- 1.25318 # 15/11/2021: country-specific median = 1.25318 for UK, NL, DE and median = 2.26803 for ES, GR
worked_hours_input        <- c(28, 38) # 15//11/21: VECTOR(female, male) = (28, 38) for UK, NL, DE and (30, 37.5) for ES, GR
working_days_lost_input   <- c(14, 27.5) # 15/11/21: VECTOR(diabetes, diabetes & CHF/stroke) = (14, 27.5) for UK, NL, DE and (10, 20) for ES, GR
cost_hour_sick_input      <- 25.781 # 15/11/21: UK=25.781 GBP, NL = 36.8 EUR, DE = 36.6 EUR, ES = 22.8 EUR, GR = 16.9 EUR
friction_period_input     <- 82.18125 # 15/11/21: UK = 82.18125 days, NL= 85 days, DE = 69 days, ES = 75 days, GR = 98.6175 days
informal_care_coef_input  <- "informal_care_coef_UK_NL_DE" # 15/11/21: regression coefficients -->> two equations informal_care_coef_UK_NL_DE and informal_care_coef_ES_GR -->> see aux_functions
employment_coef_input     <- "employment_coef_UK_NL_DE" # 15/11/21: regression coefficients -->> two equations employment_coef_UK_NL_DE and  employment_coef_ES_GR -->> see aux_functions
prod_costs_coef_input     <- "prod_costs_coef_UK_NL_DE" # 15/11/2021: regression coefficients -->> two equations prod_costs_coef_UK_NL_DE and prod_costs_coef_ES_GR   -->> see aux_functions
inf_care_age_scale_input  <-  c(72.5474088, 10.4626624) # 15/11/2021: VECTOR(mean, sd) = (72.5474088, 10.4626624) for UK, NL, DE and VECTOR(mean, sd) = (76.0286769, 9.6290315) for ES, GR
prod_loss_age_scale_input <-  c(60.2737989, 60.2737989) # 15/11/2021: VECTOR(mean, sd) = (60.2737989, 60.2737989) for UK, NL, DE and VECTOR(mean, sd) = (62.5992071, 6.6265962) for ES, GR

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
sim.results.female <- SMDMII_model_simulation(npats_input,
                                              female_input = 1,
                                              tx_cost_input,
                                              treatment_effect_HbA1c_input, 
                                              treatment_effect_HDL_input,
                                              treatment_effect_LDL_input, 
                                              treatment_effect_BMI_input,
                                              treatment_effect_SBP_input,
                                              discount_cost_input, 
                                              discount_util_input, 
                                              retirement_age_input, 
                                              inf_care_hours_input, 
                                              worked_hours_input, 
                                              working_days_lost_input, 
                                              cost_hour_sick_input, 
                                              friction_period_input, 
                                              informal_care_coef_input, 
                                              employment_coef_input, 
                                              prod_costs_coef_input, 
                                              inf_care_age_scale_input, 
                                              prod_loss_age_scale_input, 
                                              psa_input,
                                              seed_input)




sim.results.female.comp <- SMDMII_model_simulation(npats_input,
                                                   female_input = 1,
                                                   tx_cost_input,
                                                   treatment_effect_HbA1c_input = rep(0,4), 
                                                   treatment_effect_HDL_input = rep(0,4),
                                                   treatment_effect_LDL_input = rep(0,4), 
                                                   treatment_effect_BMI_input = rep(0,4),
                                                   treatment_effect_SBP_input = rep(0,4),
                                                   discount_cost_input, 
                                                   discount_util_input, 
                                                   retirement_age_input, 
                                                   psa_input,
                                                   seed_input)

sim.results.male <- SMDMII_model_simulation(npats_input,
                                            female_input = 0,
                                            tx_cost_input,
                                            treatment_effect_HbA1c_input, 
                                            treatment_effect_HDL_input,
                                            treatment_effect_LDL_input, 
                                            treatment_effect_BMI_input,
                                            treatment_effect_SBP_input,
                                            discount_cost_input, 
                                            discount_util_input, 
                                            retirement_age_input, 
                                            psa_input,
                                            seed_input)

sim.results.male.comp <- SMDMII_model_simulation(npats_input,
                                                 female_input = 0,
                                                 tx_cost_input,
                                                 treatment_effect_HbA1c_input = rep(0,4), 
                                                 treatment_effect_HDL_input = rep(0,4),
                                                 treatment_effect_LDL_input = rep(0,4), 
                                                 treatment_effect_BMI_input = rep(0,4),
                                                 treatment_effect_SBP_input = rep(0,4),
                                                 discount_cost_input, 
                                                 discount_util_input, 
                                                 retirement_age_input, 
                                                 psa_input,
                                                 seed_input)

# Print simulation duration
end <- Sys.time()
print(end - init)


# Save simulation results
save(sim.vars,
     sim.results.female,
     sim.results.female.comp,
     sim.results.male,
     sim.results.male.comp,
     file = paste0('output/', comp, '_seed_', seed_input, '.RData'))

