### OUTPUT OF ANALYSES OF 13 RANKED INTERVENTIONS ###


# DASHBOARD - INPUT PARAMTERS ---------------------------------------------
# Select country for output analysis based on identifier
# country.id <- 'UK'

# Number of runs to compile for complete output
n.runs <- 8


# LOAD FUNCTIONS AND SET WD -----------------------------------------------
# wd.dir <- getwd()

# Load stitcher function
source(paste0(wd.dir, '/R/Stitcher.R'))

# Select working directory for output based on country id
setwd(paste0(wd.dir, '/output/', country.id, '/'))

# Fraction females in target pop
frac.fem.allcountries <- data.frame(
  UK = 26706/58171, # McGovern
  NL = 0.464, # van Herp 2017
  DE = 0.472, # DPV and DIVE databases
  ES = 0.471, # Luis A Vázquez(2014)
  GR = 0.495 # Liatis 2019
)
frac.fem <- as.numeric(frac.fem.allcountries[country.id])

## Functions for automated output processing
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
              incremental_total_costs = round(mean_total_discounted_costs),
              incremental_total_qalys = round(mean_total_discounted_qalys, 3),
              incremental_undisc_qalys = round(mean_total_qalys, 3),
              headroom_20k = round(20000 * mean_total_discounted_qalys - mean_total_discounted_costs),
              headroom_50k = round(50000 * mean_total_discounted_qalys - mean_total_discounted_costs))
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
#load('All_SMI_vs_UC_basecase.RData')
# stitcher('All_SMI_vs_UC_basecase', '', name.loc2)
sewing_machine(n.runs, 'Usual_care_outcomes_run_')
uc.run.seeds <- sewing.seeds
sewing_machine(n.runs, 'All_SMI_vs_UC_basecase_run_')

output.deltas <- calc.deltas('Any SMI')
output.incidence <- calc.incidences('Any SMI', 1000)

for (i in 1:13) {
  print(i)
  sewing_machine(n.runs, paste0('Rank', i, '_basecase_run_'))
  if(!identical(sewing.seeds, uc.run.seeds)) stop(paste('Seeds for intervention', i, 'do not match those for usual care'))
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
#   load(paste0('Rank', i, '_spec_targetpop.RData'))
#   # stitcher(paste0('Rank', i, '_spec_targetpop'), '', name.loc2)
#   output.rank <- calc.deltas(paste('Rank', i, 'alt pop'))
#   output.deltas <- rbind(output.deltas, output.rank)
#   incidence.rank <- calc.incidences(paste('Rank', i, 'alt pop'), 1000)
#   output.incidence <- rbind(output.incidence, incidence.rank)
# }
# 
# # Add short effect rank 6 analysis
# rm(sim.results.female, sim.results.male, sim.results.female.comp, sim.results.male.comp)
# load('Rank6_shorteff.RData')
# output.deltas <- rbind(output.deltas, calc.deltas('Rank 6 short eff'))

print(output.deltas)
print(output.incidence)

# Sink out put to CSV
write.csv(output.deltas, paste0('Deltas_all_analyses_', country.id, '.csv'))
write.csv(output.incidence, paste0('Difference_complication_incidence_', country.id, '.csv'))


# AUTOMATED CONVERGENCE CHECK ---------------------------------------------
# Convergence Headroom - choose between files or use the stitcher

# load('All_SMI_vs_UC_basecase.RData')
# stitcher('All_SMI_vs_UC_basecase', '', name.loc2)
sewing_machine(n.runs, 'Usual_care_outcomes_run_')
sewing_machine(n.runs, 'All_SMI_vs_UC_basecase_run_')

allvuc.converge <- calc.convergence()
png(filename = 'All_vs_UC_converge.png', height = 15, width = 20, units = 'cm', res = 150)
plot(x = 1:nrow(allvuc.converge), y = allvuc.converge[ ,4], type = 'l', xlab = '# Patients (M/F weighed)', ylab = 'Incremental Headroom', main = 'Any SMI vs UC')
dev.off()

for (i in 1:13) {
  print(i)
  # load('Usual_care_outcomes.RData')
  # load(paste0('Rank', i, '_basecase.RData'))
  # stitcher('Usual_care_outcomes', '', name.loc2)
  # stitcher(paste0('Rank', i, '_basecase'), '', name.loc2)
  sewing_machine(n.runs, paste0('Rank', i, '_basecase_run_'))
  converge.rank <- calc.convergence()
  png(filename = paste0('Rank_', i, '_convergence.png'), height = 15, width = 20, units = 'cm', res = 150)
  plot(x = 1:nrow(converge.rank), y = converge.rank[ ,4], type = 'l', xlab = '# Patients (M/F weighed)', ylab = 'Incremental Headroom', main = paste('Rank', i))
  dev.off()
}

# Convergence QALYs
# load('All_SMI_vs_UC_basecase.RData')
sewing_machine(n.runs, 'Usual_care_outcomes_run_')
sewing_machine(n.runs, 'All_SMI_vs_UC_basecase_run_')

allvuc.converge <- calc.convergence()
png(filename = 'All_vs_UC_converge_QALY.png', height = 15, width = 20, units = 'cm', res = 150)
plot(x = 1:nrow(allvuc.converge), y = allvuc.converge[ ,2], type = 'l', xlab = '# Patients (M/F weighed)', ylab = 'Incremental QALY', main = 'Any SMI vs UC')
dev.off()

for (i in 1:13) {
  print(i)
  # load('Usual_care_outcomes.RData')
  # load(paste0('Rank', i, '_basecase.RData'))
  sewing_machine(n.runs, paste0('Rank', i, '_basecase_run_'))
  converge.rank <- calc.convergence()
  png(filename = paste0('Rank_', i, '_convergence_QALY.png'), height = 15, width = 20, units = 'cm', res = 150)
  plot(x = 1:nrow(converge.rank), y = converge.rank[ ,2], type = 'l', xlab = '# Patients (M/F weighed)', ylab = 'Incremental QALY', main = paste('Rank', i))
  dev.off()
}


