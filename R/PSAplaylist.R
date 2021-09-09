# Setting option for decimals
options(scipen = 3)

# Install and load required packages
# Sys.getenv() # check this for the Home path where packages are saved
.libPaths(Sys.getenv()[["R_LIBS_USER"]])

pkgs <- c("lattice", "MASS", "plyr", "survival", "tidyverse", "dplyr", "mc2d") # package names
#install.packages(pkgs)
inst <- lapply(pkgs, library, character.only = TRUE) # load them

load('R/Psa_seed.Rdata', verbose = TRUE)

for (i in 1:5) {
  gc()
  Sys.sleep(20)
  memory.size (max = F)
  
  comp <- 'UC_vs_any_SMI'
  
  treateff_hba1c <- -0.391 # Treatment effect on HbA1c (in absolute %-points HbA1c)
  treateff_hdl   <- 0.268 * 0.02586 * 10 # TRANSFORMATION *10 FOR MODEL INPUT # Treatment effect on HDL-cholesterol (absolute effect in mmol/l)
  treateff_ldl   <- -1.78 * 0.02586 * 10 # TRANSFORMATION *10 FOR MODEL INPUT# Treatment effect on LDL-cholesterol (absolute effect in mmol/l)
  treateff_bmi <- -0.283 # Treatment effect on BMI (in absolute points)
  treateff_sbp <- -2.04 / 10 # TRANSFORMATION /10 FOR MODEL INPUT# Treatment effect on SBP (in absolute mmHg)
  
  seed_input <- psa.seed[i]
  source('R/TEST_analysis.R')
}

rm(list = ls())
gc()
Sys.sleep(120)
memory.size (max = F)

load('R/Psa_seed.Rdata', verbose = TRUE)

for (i in 1:5) {
  gc()
  Sys.sleep(20)
  memory.size (max = F)
  
  comp <- 'Rank 8'
  
  treateff_hba1c <- -1.4248 # Treatment effect on HbA1c (in absolute %-points HbA1c)
  treateff_hdl <- 5.4857 * 0.02586 * 10 # TRANSFORMATION *10 FOR MODEL INPUT # Treatment effect on HDL-cholesterol (absolute effect in mmol/l)
  treateff_ldl <- -11.9313 * 0.02586 * 10 # TRANSFORMATION *10 FOR MODEL INPUT# Treatment effect on LDL-cholesterol (absolute effect in mmol/l)
  treateff_bmi <- 0 # Treatment effect on BMI (in absolute points)
  treateff_sbp <- 0  / 10 # TRANSFORMATION /10 FOR MODEL INPUT# Treatment effect on SBP (in absolute mmHg)
  
  seed_input <- psa.seed[i]
  source('R/TEST_analysis.R')
}




source('R/PSA_all_SMI_vs_UC.R')

source('R/PSA_Usual_care.R')

source('R/PSA_Rank1.R')

rm(list = ls())
gc()
Sys.sleep(120)
memory.size (max = F)


source('R/PSA_Rank2.R')

source('R/PSA_Rank3.R')

source('R/PSA_Rank4.R')

rm(list = ls())
gc()
Sys.sleep(120)
memory.size (max = F)

source('R/PSA_Rank5.R')

source('R/PSA_Rank6.R')

source('R/PSA_Rank7.R')

rm(list = ls())
gc()
Sys.sleep(120)
memory.size (max = F)

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