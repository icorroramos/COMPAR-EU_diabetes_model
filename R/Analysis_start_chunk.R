# Clear all objects and clear out memory
rm(list = ls())
gc()
Sys.sleep(120)
memory.size (max = F)


# # Setting option for decimals
# options(scipen = 3)
# 
# # Install and load required packages
# # Sys.getenv() # check this for the Home path where packages are saved
# .libPaths(Sys.getenv()[["R_LIBS_USER"]])
# 
# pkgs <- c("lattice", "MASS", "plyr", "survival", "tidyverse", "dplyr", "mc2d") # package names
# #install.packages(pkgs)
# inst <- lapply(pkgs, library, character.only = TRUE) # load them

# Load model function: 
source("R/SMI in type II diabetes - HE model v3.R")

# Load aux. functions, input parameters (from Excel), global lists, etc. 
source("R/aux_functions.R")



# ANALYSIS INPUT VARIABLES ------------------------------------------------

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

# Print time to be able to estimate finishing time
print(Sys.time())
