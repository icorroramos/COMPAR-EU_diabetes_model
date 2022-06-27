
# Setting option for decimals
options(scipen = 3)

pkgs <- c("lattice", "MASS", "plyr", "survival", "tidyverse", "dplyr", "mc2d") # package names
lapply(pkgs, library, character.only = TRUE) # load them


# SET PARAMETERS SPECIFIC TO CURRENT RUN ----------------------------------
country.id <- 'UK' # Choose from: 'UK', 'NL', 'DE', 'ES', 'GR'

# CALL ANALYSES SCRIPTS ---------------------------------------------------
npats_input   <- 2500
load('R/run_seeds.Rdata')

start.runs <- Sys.time()

for (i in 1:8) {
  run_id <- i
  seed_input <- run.seeds[i]
  source('R/Analysis_LookAHEAD.R')
  # source('R/Analysis_LookEXTREME.R')
}

end.runs <- Sys.time()
print(paste('Duration of all', i, 'runs for EuEHA scenario:'))
print(end.runs - start.runs)
