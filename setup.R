##############
### Set-up ###
##############

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