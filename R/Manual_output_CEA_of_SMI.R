# LOAD FUNCTIONS AND SET WD -----------------------------------------------
setwd('C:/Work/GitHub projects/COMPAR-EU_diabetes_model')

# DASHBOARD - INPUT PARAMTERS ---------------------------------------------
# Select country for output analysis based on identifier
country.id <- 'UK'

# Number of runs to compile for complete output
n.runs <- 8

# Choose which WTP threshold is used in the convergence plots
wtp.conversion <- 50000

# Cohortsize to calculate incidence on (i.e. incidence per x patients)
inc.cohortsize <- 1000

# LOAD FUNCTIONS AND SET WD -----------------------------------------------
# Load stitcher function
source('R/Stitcher.R')

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


# FUNCTIONS TO CALCULATE SCENARIOS ----------------------------------------

calc.short.horizon <- function(sim.output, analysis.horizon){
  full.data <- sim.output[[1]] %>% select(SIMID, SDURATION, discounted_TOTAL.COST, discounted_QALY)
  within.horizon <- full.data %>% filter(SDURATION <= analysis.horizon)
  patient.totals <- within.horizon %>% 
    group_by(SIMID) %>% 
    summarize(total_qalys = sum(discounted_QALY),
              total_costs = sum(discounted_TOTAL.COST))
  means <-patient.totals %>% summarise_all(mean) %>% select(-SIMID)
}

calc.additional.output <- function(sim.output) {
  patient.cost.data <- sim.output[[1]] %>% select(SIMID, SDURATION, IHD.COST, MI.COST, CHF.COST, STROKE.COST, AMP.COST, BLIND.COST, ULCER.COST, RENAL.COST, NOCOMP.COST, INF.CARE.COST, PROD.LOSS.COST, FUTURE.COST.MEDICAL, FUTURE.COST.NONMEDICAL, cost_discount_factor)
  discounted <- apply(patient.cost.data[, -c(1, 2, 16)], 2, '/', patient.cost.data[, 16])
  patient.cost.data.discounted <- cbind(patient.cost.data[, 1:2], discounted)
  patient.totals <- patient.cost.data.discounted %>% 
    group_by(SIMID) %>% 
    summarize(IHD_cost_disc = sum(IHD.COST),
              MI_cost_disc = sum(MI.COST), 
              CHF_cost_disc = sum(CHF.COST),
              stroke_cost_disc = sum(STROKE.COST),
              amput_cost_disc = sum(AMP.COST),
              blind_cost_disc = sum(BLIND.COST),
              ulcer_cost_disc = sum(ULCER.COST),
              renal_cost_disc = sum(RENAL.COST),
              uncomp_cost_disc = sum(NOCOMP.COST),
              informal_care_cost_disc = sum(INF.CARE.COST),
              production_losses_disc = sum(PROD.LOSS.COST),
              future_medical_cost_disc = sum(FUTURE.COST.MEDICAL),
              future_nonmedical_cost_disc = sum(FUTURE.COST.NONMEDICAL))
  mean.costs <- patient.totals %>% summarise_all(mean) %>% select(-SIMID)
  return(mean.costs)
}

calc.deltas.forD6.1 <- function(comp.name) {
  # Grab outcomes that are delivered by simulation function
  f.int.means <- as.data.frame(sim.results.female[-1])
  f.uc.means <- as.data.frame(sim.results.female.comp[-1])
  m.int.means <- as.data.frame(sim.results.male[-1])
  m.uc.means <- as.data.frame(sim.results.male.comp[-1])
  
  w.int.means <- frac.fem * f.int.means + (1 - frac.fem) * m.int.means
  w.uc.means <- frac.fem * f.uc.means + (1 - frac.fem) * m.uc.means
  
  delta.means <- w.int.means - w.uc.means
  
  # Get the additional output from the full patient histories via the calc.additional.output function
  f.int.add.cost <- calc.additional.output(sim.results.female)
  f.uc.add.cost <- calc.additional.output(sim.results.female.comp)
  m.int.add.cost <- calc.additional.output(sim.results.male)
  m.uc.add.cost <- calc.additional.output(sim.results.male.comp)
  
  w.int.add.means <- frac.fem * f.int.add.cost + (1 - frac.fem) * m.int.add.cost
  w.uc.add.means <- frac.fem * f.uc.add.cost + (1 - frac.fem) * m.uc.add.cost
  
  delta.add.means <- w.int.add.means - w.uc.add.means
  
  all.means <- cbind(delta.means, delta.add.means)
  
  relevant.out <- all.means %>% 
    transmute(comp = comp.name,
              life_expecanty = round(mean_life_expectancy, 3),
              incremental_total_qalys = round(mean_total_discounted_qalys, 3),
              related_medical_costs = round(sum(IHD_cost_disc, MI_cost_disc, CHF_cost_disc, stroke_cost_disc, amput_cost_disc, blind_cost_disc, ulcer_cost_disc, renal_cost_disc, uncomp_cost_disc)),
              unrelated_medical_costs = round(future_medical_cost_disc),
              informal_care_costs = round(informal_care_cost_disc),
              productivity_losses = round(production_losses_disc),
              future_nonmedical_cons = round(future_nonmedical_cost_disc),
              incremental_total_costs = round(mean_total_discounted_costs),
              headroom_20k = round(20000 * mean_total_discounted_qalys - mean_total_discounted_costs),
              headroom_50k = round(50000 * mean_total_discounted_qalys - mean_total_discounted_costs)
    )
}

# CALCULATE SCNARIO OUTCOMES ----------------------------------------------

# Results for scenario with extended treatment effect
sewing_machine(n.runs, 'output/UK/Usual_care_outcomes_run_') # Usual care runs
sewing_machine(n.runs, 'output/Scenarios/All_SMI_vs_UC_extended_effect_')
output.extended <- calc.deltas.forD6.1('Any SMI extended effect')

# Results for scenario with lifetime effect
sewing_machine(n.runs, 'output/UK/Usual_care_outcomes_run_') # Usual care runs
sewing_machine(n.runs, 'output/Scenarios/All_SMI_vs_UC_lifetime_effect_')
output.lifetime <- calc.deltas.forD6.1('Any SMI lifetime effect')


# Calculating short time horizons
sewing_machine(n.runs, paste0('output/', country.id, '/Usual_care_outcomes_run_'))
uc.run.seeds <- sewing.seeds
sewing_machine(n.runs, paste0('output/', country.id, '/All_SMI_vs_UC_basecase_run_'))


horizon.scenario <- 10

f.int.sh <- calc.short.horizon(sim.results.female, horizon.scenario)
f.uc.sh <- calc.short.horizon(sim.results.female.comp, horizon.scenario)
m.int.sh <- calc.short.horizon(sim.results.male, horizon.scenario)
m.uc.sh <- calc.short.horizon(sim.results.male.comp, horizon.scenario)

w.int.means.sh <- frac.fem * f.int.sh + (1 - frac.fem) * m.int.sh
w.uc.means.sh <- frac.fem * f.uc.sh + (1 - frac.fem) * m.uc.sh

delta.means.sh <- w.int.means.sh - w.uc.means.sh


sh.scenario.result <- delta.means.sh %>% 
  transmute(incremental_total_qalys = round(total_qalys, 3),
            incremental_total_costs = round(total_costs),
            headroom_20k = round(20000 * total_qalys - total_costs),
            headroom_50k = round(50000 * total_qalys - total_costs))

print(sh.scenario.result)



# SCENARIO FOR EUHEA PRESENTATION -----------------------------------------
sewing_machine(n.runs, 'output/UK/Usual_care_outcomes_run_') # Usual care runs
sewing_machine(n.runs, 'output/UK/LookAHEAD_run_')
output.lookahead <- calc.deltas.forD6.1('LookAHEAD trial effect')

# Extreme effect scenario
sewing_machine(n.runs, 'output/UK/Usual_care_outcomes_run_') # Usual care runs
sewing_machine(n.runs, 'output/UK/LookEXTREME_run_')
output.extreme <- calc.deltas.forD6.1('Extreme LookAHED trial effect')


# TEMP DUMP
data.out <- rbind(output.lookahead, output.extreme)
write.csv(data.out, 'Output_EUHEA_scenarios.csv', quote = FALSE)

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
female.pp.uc <- sim.results.female.comp$simulation_patients_history %>% 
  group_by(SIMID) %>% 
  summarise(life_expectancy = max(SDURATION), 
            total_qaly = sum(QALY),
            total_cost = sum(TOTAL.COST)
  )

male.pp.uc <- sim.results.male.comp$simulation_patients_history %>% 
  group_by(SIMID) %>% 
  summarise(life_expectancy = max(SDURATION), 
            total_qaly = sum(QALY),
            total_cost = sum(TOTAL.COST)
  )

female.pp.int <- sim.results.female$simulation_patients_history %>% 
  group_by(SIMID) %>% 
  summarise(life_expectancy = max(SDURATION), 
            total_qaly = sum(QALY),
            total_cost = sum(TOTAL.COST)
  )

male.pp.int <- sim.results.male$simulation_patients_history %>% 
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




# Old automated 
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