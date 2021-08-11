# Setting option for decimals
options(scipen = 3)

# Install and load required packages
# Sys.getenv() # check this for the Home path where packages are saved
.libPaths(Sys.getenv()[["R_LIBS_USER"]])

pkgs <- c("lattice", "MASS", "plyr", "survival", "tidyverse", "dplyr", "mc2d") # package names
#install.packages(pkgs)
inst <- lapply(pkgs, library, character.only = TRUE) # load them


# source('R/PSA_all_SMI_vs_UC.R')
# 
# source('R/PSA_Rank1.R')
# 
# source('R/PSA_Rank2.R')

source('R/PSA_Rank3.R')

source('R/PSA_Rank4.R')

source('R/PSA_Rank5.R')

source('R/PSA_Rank6.R')

rm(list = ls())
gc()
Sys.sleep(120)
memory.size (max = F)

source('R/PSA_Rank7.R')

source('R/PSA_Rank8.R')

source('R/PSA_Rank9.R')

source('R/PSA_Rank10.R')

rm(list = ls())
gc()
Sys.sleep(120)
memory.size (max = F)

source('R/PSA_Rank11.R')

source('R/PSA_Rank12.R')

source('R/PSA_Rank13.R')