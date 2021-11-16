### OUTPUT OF ANALYSES OF 13 RANKED INTERVENTIONS ###

# Load stitcher function
source('R/Stitcher.R')

# Select working directory to choose outputs from
# setwd('output')
# setwd('output/Deterministic_1000_pats')
# setwd('output/UK/Lisa_15000')
setwd('output/DE/')

name.loc2 <- '5000_pats_seed-77/' # Name of subfolder with second run results

# INPUT PARAMETERS FOR CALCULATIONS ---------------------------------------
# Fraction females in target pop
frac.fem <- 0.472 # UK: 26706/58171 

#Lambda (willingness to pay)
wtp <- 20000

# FUNCTIONS FOR AUTOMATED OUTPUT CALCULATIONS -----------------------------

calc.deltas <- function(comp.name) {
  f.int.means <- as.data.frame(sim.results.female[-1])
  f.uc.means <- as.data.frame(sim.results.female.comp[-1])
  m.int.means <- as.data.frame(sim.results.male[-1])
  m.uc.means <- as.data.frame(sim.results.male.comp[-1])
  
  w.int.means <- frac.fem * f.int.means + (1 - frac.fem) * m.int.means
  w.uc.means <- frac.fem * f.uc.means + (1 - frac.fem) * m.uc.means
  
  delta.means <- w.int.means - w.uc.means
  
  relevant.out <- delta.means %>% 
    transmute(comp = comp.name,
              added_days_life_expecanty = round(mean_life_expectancy * 365.25, 1),
              incremental_total_costs = round(mean_total_costs),
              incremental_total_qalys = round(mean_total_qalys, 3),
              headroom = round(wtp * mean_total_qalys - mean_total_costs))
}

calc.incidences <- function(comp.name, per.popsize) {

  f.int.means <- as.data.frame(sim.results.female[-1])
  f.uc.means <- as.data.frame(sim.results.female.comp[-1])
  m.int.means <- as.data.frame(sim.results.male[-1])
  m.uc.means <- as.data.frame(sim.results.male.comp[-1])
  
  w.int.means <- frac.fem * f.int.means + (1 - frac.fem) * m.int.means
  w.uc.means <- frac.fem * f.uc.means + (1 - frac.fem) * m.uc.means
  
  delta.means <- w.int.means - w.uc.means
  
  incidences <- delta.means %>% 
    transmute(intervention = comp.name,
              diff_CHD = mean_CHF_rate * per.popsize,
              diff_MI = mean_MI_rate * per.popsize,
              diff_stroke = mean_STROKE_rate * per.popsize,
              diff_blindness = mean_BLIND_rate * per.popsize,
              diff_renal_fail = mean_RENAL_rate * per.popsize,
              diff_amputation = (mean_AMP1_rate + mean_AMP2_rate) * per.popsize,
              diff_ulcer = mean_ULCER_rate * per.popsize
  )
  return(incidences)
}

calc.convergence <- function() {
  female.pp.uc <- sim.results.female.comp$simulation_patients_history %>% 
    group_by(SIMID) %>% 
    summarise(life_expectancy = max(SDURATION), 
              total_qaly = sum(QALY),
              total_cost = sum(TOTAL.COST),
              .groups = 'drop'
    )
  
  male.pp.uc <- sim.results.male.comp$simulation_patients_history %>% 
    group_by(SIMID) %>% 
    summarise(life_expectancy = max(SDURATION), 
              total_qaly = sum(QALY),
              total_cost = sum(TOTAL.COST),
              .groups = 'drop'
    )
  
  female.pp.int <- sim.results.female$simulation_patients_history %>% 
    group_by(SIMID) %>% 
    summarise(life_expectancy = max(SDURATION), 
              total_qaly = sum(QALY),
              total_cost = sum(TOTAL.COST),
              .groups = 'drop'
    )
  
  male.pp.int <- sim.results.male$simulation_patients_history %>% 
    group_by(SIMID) %>% 
    summarise(life_expectancy = max(SDURATION), 
              total_qaly = sum(QALY),
              total_cost = sum(TOTAL.COST),
              .groups = 'drop'
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
  return(delta.weighed.conversion)
}


# AUTOMATED PROCESSING OF SIM OUTPUT --------------------------------------
load('All_SMI_vs_UC_basecase.RData')
# stitcher('All_SMI_vs_UC_basecase', '', name.loc2)

output.deltas <- calc.deltas('Any SMI')
output.incidence <- calc.incidences('Any SMI', 1000)

for (i in 1:13) {
  print(i)
  load('Usual_care_outcomes.RData')
  load(paste0('Rank', i, '_basecase.RData'))
  # stitcher('Usual_care_outcomes', '', name.loc2)
  # stitcher(paste0('Rank', i, '_basecase'), '', name.loc2)
  output.rank <- calc.deltas(paste('Rank', i))
  output.deltas <- rbind(output.deltas, output.rank)
  incidence.rank <- calc.incidences(paste('Rank', i), 1000)
  output.incidence <- rbind(output.incidence, incidence.rank)
}

# # Remove all sim results to prevent use of base case comparator in alternate analyses
# rm(sim.results.female, sim.results.male, sim.results.female.comp, sim.results.male.comp)
# 
# for (i in c(1, 2)) {
#   print(i)
#   # load(paste0('Rank', i, '_spec_targetpop.RData'))
#   stitcher(paste0('Rank', i, '_spec_targetpop'), '', name.loc2)
#   output.rank <- calc.deltas(paste('Rank', i, 'alt pop'))
#   output.deltas <- rbind(output.deltas, output.rank)
#   incidence.rank <- calc.incidences(paste('Rank', i, 'alt pop'), 1000)
#   output.incidence <- rbind(output.incidence, incidence.rank)
# }

# # Add short effect rank 6 analysis
# rm(sim.results.female, sim.results.male, sim.results.female.comp, sim.results.male.comp)
# load('Rank6_shorteff.RData')
# output.deltas <- rbind(output.deltas, calc.deltas('Rank 6 short eff'))

print(output.deltas)
print(output.incidence)

# Sink out put to CSV
write.csv(output.deltas, 'Deltas_all_analyses.csv')
write.csv(output.incidence, 'Difference_complication_incidence.csv')


# AUTOMATED CONVERGENCE CHECK ---------------------------------------------
# Convergence Headroom - choose between files or use the stitcher

# load('All_SMI_vs_UC_basecase.RData')
stitcher('All_SMI_vs_UC_basecase', '', name.loc2)
allvuc.converge <- calc.convergence()
png(filename = 'All_vs_UC_converge.png', height = 15, width = 20, units = 'cm', res = 150)
plot(x = 1:nrow(allvuc.converge), y = allvuc.converge[ ,4], type = 'l', xlab = '# Patients (M/F weighed)', ylab = 'Incremental Headroom', main = 'Any SMI vs UC')
dev.off()

for (i in 1:13) {
  print(i)
  # load('Usual_care_outcomes.RData')
  # load(paste0('Rank', i, '_basecase.RData'))
  stitcher('Usual_care_outcomes', '', name.loc2)
  stitcher(paste0('Rank', i, '_basecase'), '', name.loc2)
  converge.rank <- calc.convergence()
  png(filename = paste0('Rank_', i, '_convergence.png'), height = 15, width = 20, units = 'cm', res = 150)
  plot(x = 1:nrow(converge.rank), y = converge.rank[ ,4], type = 'l', xlab = '# Patients (M/F weighed)', ylab = 'Incremental Headroom', main = paste('Rank', i))
  dev.off()
}

# Convergence QALYs
load('All_SMI_vs_UC_basecase.RData')
allvuc.converge <- calc.convergence()
png(filename = 'All_vs_UC_converge_QALY.png', height = 15, width = 20, units = 'cm', res = 150)
plot(x = 1:nrow(allvuc.converge), y = allvuc.converge[ ,2], type = 'l', xlab = '# Patients (M/F weighed)', ylab = 'Incremental QALY', main = 'Any SMI vs UC')
dev.off()

for (i in 1:13) {
  print(i)
  load('Usual_care_outcomes.RData')
  load(paste0('Rank', i, '_basecase.RData'))
  converge.rank <- calc.convergence()
  png(filename = paste0('Rank_', i, '_convergence_QALY.png'), height = 15, width = 20, units = 'cm', res = 150)
  plot(x = 1:nrow(converge.rank), y = converge.rank[ ,2], type = 'l', xlab = '# Patients (M/F weighed)', ylab = 'Incremental QALY', main = paste('Rank', i))
  dev.off()
}

# MANUAL ANALYSIS ---------------------------------------------------------

# Paste datasets with different seeds together
load('All_SMI_vs_UC_basecase.RData', verbose = TRUE)
res.names <- c('female', 'female.comp', 'male', 'male.comp')
for (i in 1:length(res.names)) {
  sim.data <- get(paste0('sim.results.', res.names[i]))
  pat.hist <- sim.data[[1]]
  pat.hist$SIMID <- pat.hist$SIMID + 5000
  assign(paste0('sim.results.', res.names[i], '.extended'), pat.hist)
}

load('All_SMI_vs_UC_basecase_5000.RData', verbose = TRUE)
for (i in 1:length(res.names)) {
  sim.data <- get(paste0('sim.results.', res.names[i]))
  first <- sim.data[[1]]
  second <- get(paste0('sim.results.', res.names[i], '.extended'))
  combo <- rbind(first, second)
  assign(paste0('sim.results.', res.names[i]), combo)
}


# LOAD SIMULATION DATA ----------------------------------------------------

# load('All_SMI_vs_UC_basecase.RData', verbose = TRUE)
# load('output/Usual_care_outcomes.RData', verbose = TRUE)
load('Rank1_basecase.RData', verbose = TRUE)
# load('Rank2_basecase.RData', verbose = TRUE)
# load('output/Rank3_basecase.RData', verbose = TRUE)
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

# load('All_SMI_vs_UC_seed77.RData', verbose = TRUE)
# load('All_SMI_vs_UC_seed1984.RData', verbose = TRUE)
# load('All_SMI_vs_UC_seed15.RData', verbose = TRUE)
load('All_SMI_vs_UC_seed265979.RData', verbose = TRUE)


# load('Usual_care_outcomes_5000.RData', verbose = TRUE)
# load('Usual_care_outcomes_4000_pats.RData', verbose = TRUE)
# load('Rank3_basecase_4000_pats.RData', verbose = TRUE)
# 
# load('Usual_care_outcomes_2000_pats_altseed.RData', verbose = TRUE)
# load('Rank3_basecase_2000_pats_altseed.RData', verbose = TRUE)

# ANALYSIS OF MEANS -------------------------------------------------------

f.int.means <- as.data.frame(sim.results.female[-1])
f.uc.means <- as.data.frame(sim.results.female.comp[-1])
m.int.means <- as.data.frame(sim.results.male[-1])
m.uc.means <- as.data.frame(sim.results.male.comp[-1])

w.int.means <- frac.fem * f.int.means + (1 - frac.fem) * m.int.means
w.uc.means <- frac.fem * f.uc.means + (1 - frac.fem) * m.uc.means

delta.means <- w.int.means - w.uc.means

relevant.out <- delta.means %>%
  transmute(life_expecanty = round(mean_life_expectancy, 3),
            total_costs = round(mean_total_costs),
            total_qalys = round(mean_total_qalys, 3),
            headroom = round(wtp * mean_total_qalys - mean_total_costs, 2))


print(relevant.out)



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