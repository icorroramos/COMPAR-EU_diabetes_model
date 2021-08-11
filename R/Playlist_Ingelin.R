# Setting option for decimals
options(scipen = 3)

# Install and load required packages
# Sys.getenv() # check this for the Home path where packages are saved
.libPaths(Sys.getenv()[["R_LIBS_USER"]])

pkgs <- c("lattice", "MASS", "plyr", "survival", "tidyverse", "dplyr", "mc2d") # package names
#install.packages(pkgs)
inst <- lapply(pkgs, library, character.only = TRUE) # load them

source('R/PSA_Rank11.R')

source('R/PSA_Rank12.R')

source('R/PSA_Rank13.R')