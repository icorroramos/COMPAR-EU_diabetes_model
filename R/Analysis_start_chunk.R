# # Clear all objects and clear out memory
# rm(list = ls())
# gc()
# Sys.sleep(90)
# memory.size (max = F)


# ANALYSIS INPUT VARIABLES ------------------------------------------------
### Variables moving to SURF scripts for use on SURF-Sara Lisa
# # Set a run identifier
# run_id <- 1
# 
# # Set random seed for replication purposes
# seed_input <- 958 # A random seed that it is used to ensure consistency in the model results.
# 
# # # Please select the number of patients in simulation (default 1000 in deterministic run)
# npats_input   <- 25 # 1000
#
# # Country identifier
# country.id <- 'GR' # Choose from: 'UK', 'NL', 'DE', 'ES', 'GR'

#--------


# Duration of the effect from SMIs
treateff_start   <- 1 # Cycle in which treatment effect starts
treateff_end     <- 3 # Cycle in which treatment effect ends
treateff_decline <- 2 # Cycle in which treatment effect starts to decline (linearly)


# Please select running mode: 0 = deterministic, 1 = PSA. Default: 0
psa_input <- 0

# Please select number of PSA runs (only works if psa_input <- 1, otherwise will be ignored ). Default: 500 --> NOTE: not completely implemented
n_psa_input <- 5

# Please indicate the name of the treatment to be identified
tx_label <- "Usual care"

# Cost inputs
tx_cost_input <- 0 # Total treatment cost --> Not sure if here or in Excel.     


# LOAD MODEL FUNCTION AND INPUT DATA --------------------------------------

# Load model function: 
source("R/SMI in type II diabetes - HE model v3.R")

# Load aux. functions, input parameters (from Excel), global lists, etc. 
source("R/aux_functions.R")


# Print time to be able to estimate finishing time
print(Sys.time())
