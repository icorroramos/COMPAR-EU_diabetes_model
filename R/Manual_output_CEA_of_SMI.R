# LOAD FUNCTIONS AND SET WD -----------------------------------------------
wd.dir <- getwd()

# Load stitcher function
source(paste0(wd.dir, '/R/Stitcher.R'))


## Override automatic working directory setting for manual selection of output to analyze
setwd("C:/Work/GitHub projects/COMPAR-EU_diabetes_model")
#setwd('\\\\mydocuments.eur.nl@SSL/DavWWWRoot/shared/groups/IMTA-600-H2020COMPAR-EU/600 Projectuitvoer/Disease models/Diabetes/Testing analysis settings')
name.loc2 <- '5000_pats_seed-77/' # Name of subfolder with second run results (for use with stitcher function)

# Fraction females in target pop
frac.fem.allcountries <- data.frame(
  UK = 26706/58171, # McGovern
  NL = 0.464, # van Herp 2017
  DE = 0.472, # DPV and DIVE databases
  ES = 0.471, # Luis A Vázquez(2014)
  GR = 0.495 # Liatis 2019
)
frac.fem <- as.numeric(frac.fem.allcountries[country.id])

wtp <- 20000

# MANUAL ANALYSIS ---------------------------------------------------------

# LOAD SIMULATION DATA ----------------------------------------------------

# load('All_SMI_vs_UC_basecase.RData', verbose = TRUE)
# load('Usual_care_outcomes.RData', verbose = TRUE)
load('Rank1_basecase.RData', verbose = TRUE)
# load('Rank2_basecase.RData', verbose = TRUE)
# load('Rank3_basecase.RData', verbose = TRUE)
# load('Rank4_basecase.RData', verbose = TRUE)
# load('Rank5_basecase.RData', verbose = TRUE)
# load('Rank6_basecase.RData', verbose = TRUE)
# load('Rank7_basecase.RData', verbose = TRUE)
# load('Rank8_basecase.RData', verbose = TRUE)
# load('Rank9_basecase.RData', verbose = TRUE)
# load('Rank10_basecase.RData', verbose = TRUE)
# load('Rank11_basecase.RData', verbose = TRUE)
# load('Rank12_basecase.RData', verbose = TRUE)
# load('Rank13_basecase.RData', verbose = TRUE)

# load('Rank1_spec_targetpop.RData', verbose = TRUE)
# load('Rank6_spec_targetpop.RData', verbose = TRUE)

# load('Rank6_shorteff.RData', verbose = TRUE)
# load('Rank8_shorter_treateff.RData', verbose = TRUE)
# load('Rank8_longer_treateff.RData', verbose = TRUE)

load('All_SMI_vs_UC_basecase_run_1.RData', verbose = TRUE)

load('Usual_care_outcomes_run_1.RData', verbose = TRUE)

load('Rank1_basecase_run_1.RData', verbose = TRUE)

# ANALYSIS OF MEANS -------------------------------------------------------

f.int.means <- as.data.frame(sim.results.female[-1])
f.uc.means <- as.data.frame(sim.results.female.comp[-1])
m.int.means <- as.data.frame(sim.results.male[-1])
m.uc.means <- as.data.frame(sim.results.male.comp[-1])

f.delta <- f.int.means - f.uc.means
m.delta <- m.int.means - m.uc.means

w.int.means <- frac.fem * f.int.means + (1 - frac.fem) * m.int.means
w.uc.means <- frac.fem * f.uc.means + (1 - frac.fem) * m.uc.means

delta.means <- w.int.means - w.uc.means

relevant.out <- delta.means %>%
  transmute(life_expecanty = round(mean_life_expectancy, 3),
            total_costs = round(mean_total_discounted_costs),
            total_qalys = round(mean_total_discounted_qalys, 3),
            headroom = round(wtp * mean_total_discounted_qalys - mean_total_discounted_costs, 2))

print(relevant.out)

compare.sex <- data.frame(female = t(f.delta), male = t(m.delta))
print(compare.sex)


# Code for old simulation function output
relevant.out.old <- delta.means %>%
  transmute(life_expecanty = round(mean_life_expectancy, 3),
            total_costs = round(mean_total_costs),
            total_qalys = round(mean_total_qalys, 3),
            headroom = round(wtp * mean_total_qalys - mean_total_costs, 2))

print(relevant.out.old)

compare.sex <- data.frame(female = t(f.delta), male = t(m.delta))
print(compare.sex)

# ANALYSIS OF INNER LOOP --------------------------------------------------
# FOR EXTENDED DATASETS
female.pp.uc <- sim.results.female.comp %>% 
  group_by(SIMID) %>% 
  summarise(life_expectancy = max(SDURATION), 
            total_qaly = sum(QALY),
            total_cost = sum(TOTAL.COST)
  )

male.pp.uc <- sim.results.male.comp %>% 
  group_by(SIMID) %>% 
  summarise(life_expectancy = max(SDURATION), 
            total_qaly = sum(QALY),
            total_cost = sum(TOTAL.COST)
  )

female.pp.int <- sim.results.female %>% 
  group_by(SIMID) %>% 
  summarise(life_expectancy = max(SDURATION), 
            total_qaly = sum(QALY),
            total_cost = sum(TOTAL.COST)
  )

male.pp.int <- sim.results.male %>% 
  group_by(SIMID) %>% 
  summarise(life_expectancy = max(SDURATION), 
            total_qaly = sum(QALY),
            total_cost = sum(TOTAL.COST)
  )

female.uc.conversion <- matrix(nrow = nrow(female.pp.uc), ncol = 3)
male.uc.conversion <- matrix(nrow = nrow(male.pp.int), ncol = 3)
female.int.conversion <- matrix(nrow = nrow(female.pp.int), ncol = 3)
male.int.conversion <- matrix(nrow = nrow(male.pp.int), ncol = 3)

for (i in 1:nrow(female.pp.uc)) {
  female.uc.conversion[i, ] <- colSums(female.pp.uc[1:i, -1] / i)
  male.uc.conversion[i, ] <- colSums(male.pp.uc[1:i, -1] / i)
  female.int.conversion[i, ] <- colSums(female.pp.int[1:i, -1] / i)
  male.int.conversion[i, ] <- colSums(male.pp.int[1:i, -1] / i)
}

delta.female.conversion <- female.int.conversion - female.uc.conversion
delta.female.conversion <- cbind(delta.female.conversion, delta.female.conversion[, 2] * wtp - delta.female.conversion[, 3])
delta.male.conversion <- male.int.conversion - male.uc.conversion
delta.male.conversion <- cbind(delta.male.conversion, delta.male.conversion[, 2] * wtp - delta.male.conversion[, 3])

delta.weighed.conversion <- frac.fem * delta.female.conversion + (1 - frac.fem) * delta.male.conversion

plot(x = 1:nrow(delta.female.conversion), y = delta.female.conversion[ ,2], type = 'l', xlab = '# Females', ylab = 'Incremental QALYs')
plot(x = 1:nrow(delta.female.conversion), y = delta.female.conversion[ ,4], type = 'l', xlab = '# Females', ylab = 'Incremental Headroom')

plot(x = 1:nrow(delta.weighed.conversion), y = delta.weighed.conversion[ ,4], type = 'l', xlab = '# Patients (M/F weighed', ylab = 'Incremental Headroom')


# CONVERGENCE CHECKS ON SEED VALUES ---------------------------------------

for (i in 1:5) {
  print(i)
  load(paste0('output/UC_vs_any_SMI_seed_', psa.seed[i], '.RData'))
  seed.converge <- calc.convergence()
  # Plots of incremental QALY (weighed)
  png(filename = paste0('output/Converge_Incr_QALY_UCvSMI_seed_', psa.seed[i], '-', i, '.png'), height = 15, width = 20, units = 'cm', res = 150)
  plot(x = 1:nrow(seed.converge), y = seed.converge[ ,2], type = 'l', xlab = '# Patients (M/F weighed)', ylab = 'Incremental QALY', main = paste('Any SMI vs UC - seed', psa.seed[i]))
  dev.off()
  # Plots of headroom
  png(filename = paste0('output/Converge_Headroom_UCvSMI_seed_', psa.seed[i], '-', i, '.png'), height = 15, width = 20, units = 'cm', res = 150)
  plot(x = 1:nrow(seed.converge), y = seed.converge[ ,4], type = 'l', xlab = '# Patients (M/F weighed)', ylab = 'Headroom', main = paste('Any SMI vs UC - seed', psa.seed[i]))
  dev.off()
}

for (i in 1:5) {
  print(i)
  load(paste0('output/Rank 8_seed_', psa.seed[i], '.RData'))
  seed.converge <- calc.convergence()
  # Plots of incremental QALY (weighed)
  png(filename = paste0('output/Converge_Incr_QALY_Rank8_seed_', psa.seed[i], '-', i, '.png'), height = 15, width = 20, units = 'cm', res = 150)
  plot(x = 1:nrow(seed.converge), y = seed.converge[ ,2], type = 'l', xlab = '# Patients (M/F weighed)', ylab = 'Incremental QALY', main = paste('Rank 8 - seed', psa.seed[i]))
  dev.off()
  # Plots of headroom
  png(filename = paste0('output/Converge_Headroom_Rank8_seed_', psa.seed[i], '-', i, '.png'), height = 15, width = 20, units = 'cm', res = 150)
  plot(x = 1:nrow(seed.converge), y = seed.converge[ ,4], type = 'l', xlab = '# Patients (M/F weighed)', ylab = 'Headroom', main = paste('Rank 8 - seed', psa.seed[i]))
  dev.off()
}





# DETAILED ANALYSIS -------------------------------------------------------
# diff.fem <- data.frame(
#   Usual_care = t(f.uc.means),
#   SMI = t(f.int.means),
#   Incr = t(f.int.means) - t(f.uc.means)
# )
# 
# diff.male <- data.frame(
#   Usual_care = t(m.uc.means),
#   SMI = t(m.int.means),
#   Incr = t(m.int.means) - t(m.uc.means)
# )
# 
# print(relevant.out)
# print(diff.fem)
# print(diff.male)
# 
# 
# write.csv(diff.fem, 'Deltas_female.csv')
# write.csv(diff.male, 'Deltas_male.csv')