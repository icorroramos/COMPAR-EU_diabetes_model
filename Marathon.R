### SCRIPT TO RUN ALL COUNTRIES & ALL ANALYSES LOCALLY 

# Setting option for decimals
options(scipen = 3)


pkgs <- c("lattice", "MASS", "plyr", "survival", "tidyverse", "dplyr", "mc2d") # package names
lapply(pkgs, library, character.only = TRUE) # load them


countries.torun <- c('UK', 'NL', 'DE', 'ES', 'GR')

# CALL ANALYSES SCRIPTS ---------------------------------------------------
npats_input   <- 2
load('R/run_seeds.Rdata')

start.marathon <- Sys.time()

for (h in 1:length(countries.torun)) {
  start.country <- Sys.time()
  country.id <- countries.torun[h]
  print(paste('Sarting analysis for', country.id))
  
  # UC run
  start.runs <- Sys.time()
  for (i in 1:8) {
    run_id <- i
    seed_input <- run.seeds[i]
    source('R/Analysis_Usual_Care.R')
  }
  
  end.runs <- Sys.time()
  print(paste('Duration of all', i, 'runs for Usual Care:'))
  print(end.runs - start.runs)
  print('')
  
  # Any SMI vs UC
  start.runs <- Sys.time()
  for (i in 1:8) {
    run_id <- i
    seed_input <- run.seeds[i]
    source('R/Analysis_all_SMI_vs_UC.R')
  }
  
  end.runs <- Sys.time()
  print(paste('Duration of all', i, 'runs for Any SMI:'))
  print(end.runs - start.runs)
  print('')
  
  # All numbered interventions
  for (j in 1:13) {
    start.runs <- Sys.time()
    
    for (i in 1:8) {
      run_id <- i
      seed_input <- run.seeds[i]
      source(paste0('R/Analysis_Rank', j, '.R'))
    }
    
    end.runs <- Sys.time()
    print(paste('Duration of all', i, 'runs for Rank', j, ':'))
    print(end.runs - start.runs)
    print('')
  }
  
  end.country <-Sys.time()
  print(paste('Duration of all analyses for:', country.id))
  print(end.country - start.country)
  
  
  # Clear memory before going to next country
  rm(list=setdiff(ls(), c('countries.torun', 'npats_input', 'run.seeds', 'start.marathon', 'h')))
  gc()
  Sys.sleep(10) #FIXME: reset to 120
  memory.size (max = F)
}

end.marathon <-Sys.time()
print(paste('Duration of complete marathon'))
print(end.marathon - start.marathon)
