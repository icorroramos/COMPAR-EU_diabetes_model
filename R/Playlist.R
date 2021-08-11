# Setting option for decimals
options(scipen = 3)

# Install and load required packages
# Sys.getenv() # check this for the Home path where packages are saved
.libPaths(Sys.getenv()[["R_LIBS_USER"]])

pkgs <- c("lattice", "MASS", "plyr", "survival", "tidyverse", "dplyr", "mc2d") # package names
#install.packages(pkgs)
inst <- lapply(pkgs, library, character.only = TRUE) # load them

source('R/Analysis_all_SMI_vs_UC.R')

source('R/Analysis_Usual_Care.R')

source('R/Analysis_Rank1.R')

source('R/Analysis_Rank2.R')

source('R/Analysis_Rank3.R')

source('R/Analysis_Rank4.R')

source('R/Analysis_Rank5.R')

rm(list = ls())
gc()
Sys.sleep(120)
memory.size (max = F)

source('R/Analysis_Rank6.R')

source('R/Analysis_Rank7.R')

source('R/Analysis_Rank8.R')

source('R/Analysis_Rank9.R')

source('R/Analysis_Rank10.R')

source('R/Analysis_Rank11.R')

source('R/Analysis_Rank12.R')

rm(list = ls())
gc()
Sys.sleep(120)
memory.size (max = F)

source('R/Analysis_Rank13.R')

source('R/Analysis_Rank1_spec_targetpop.R')

source('R/Analysis_Rank2_spec_targetpop.R')

source('R/Analysis_Rank6_spec_targetpop.R')

source('R/Analysis_Rank6_shorteff.R')

