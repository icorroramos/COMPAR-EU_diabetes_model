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

pkgs <- c("lattice", "MASS", "plyr", "survival", "tidyverse", "dplyr") # package names
#install.packages(pkgs)
inst <- lapply(pkgs, library, character.only = TRUE) # load them

# Load model function: 
source("R/SMI in type II diabetes - HE model v3.R")

# Load aux. functions, input parameters (from Excel), global lists, etc. 
source("R/aux_functions.R")

# Variable defined to keep track of simulation time (delete afterwards)
init <- Sys.time()


# User adjustable settings #

# Adjust target population
# baseline_characteristics <- read.csv("input/UK/baseline_characteristics_UK_rank_1_study_pop.csv", sep=",")

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
tx_label <- "Usual care"

# Treatment effect inputs
treateff_start   <- 1 # Cycle in which treatment effect starts
treateff_end     <- 3 # Cycle in which treatment effect ends
treateff_decline <- 2 # Cycle in which treatment effect starts to decline (linearly)

treateff_hba1c <- -0.5616 # Treatment effect on HbA1c (in absolute %-points HbA1c)
treateff_hdl   <- 0  * 10 # TRANSFORMATION *10 FOR MODEL INPUT # Treatment effect on HDL-cholesterol (absolute effect in mmol/l)
treateff_ldl   <- 0  * 10 # TRANSFORMATION *10 FOR MODEL INPUT# Treatment effect on LDL-cholesterol (absolute effect in mmol/l)
treateff_bmi <- 0 # Treatment effect on BMI (in absolute points)
treateff_sbp <- -10.5144 / 10 # TRANSFORMATION /10 FOR MODEL INPUT# Treatment effect on SBP (in absolute mmHg)

# # Treatment effects currently not in use
# treateff_QoL

# Tx effects are vectors: the current assumption is that the same start, end and decline is assumed for all effect modifiers
treatment_effect_HbA1c_input <- c(treateff_hba1c, treateff_start, treateff_end, treateff_decline)
treatment_effect_HDL_input <- c(treateff_hdl, treateff_start, treateff_end, treateff_decline)
treatment_effect_LDL_input <- c(treateff_ldl, treateff_start, treateff_end, treateff_decline)
treatment_effect_BMI_input <- c(treateff_bmi, treateff_start, treateff_end, treateff_decline)
treatment_effect_SBP_input <- c(treateff_sbp, treateff_start, treateff_end, treateff_decline)


# Cost inputs
tx_cost_input <- 0 # Total treatment cost --> Not sure if here or in Excel.     
retirement_age_input <- 66

# Set random seed for replication purposes
seed_input <- 77 # A random seed that it is used to ensure consistency in the model results. 

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


# Print simulation duration
end <- Sys.time()
print(end - init)


# Save simulation results
save(sim.vars,
     sim.results.female,
     sim.results.male,
     file = 'output/Rank5_basecase.RData')  #file = 'output/Rank1_basecase.RData') # file = 'output/Rank1_spec_tartetpop.RData') #



