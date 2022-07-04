### OUTPUT OF ANALYSES OF 13 RANKED INTERVENTIONS ###
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

## Functions for automated output processing
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
              incremental_total_qalys = round(mean_total_discounted_qalys, 3),
              incremental_total_costs = round(mean_total_discounted_costs),
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
              total_qaly = sum(discounted_QALY),
              total_cost = sum(discounted_TOTAL.COST),
              .groups = 'drop'
    )
  
  male.pp.uc <- sim.results.male.comp$simulation_patients_history %>% 
    group_by(SIMID) %>% 
    summarise(life_expectancy = max(SDURATION), 
              total_qaly = sum(discounted_QALY),
              total_cost = sum(discounted_TOTAL.COST),
              .groups = 'drop'
    )
  
  female.pp.int <- sim.results.female$simulation_patients_history %>% 
    group_by(SIMID) %>% 
    summarise(life_expectancy = max(SDURATION), 
              total_qaly = sum(discounted_QALY),
              total_cost = sum(discounted_TOTAL.COST),
              .groups = 'drop'
    )
  
  male.pp.int <- sim.results.male$simulation_patients_history %>% 
    group_by(SIMID) %>% 
    summarise(life_expectancy = max(SDURATION), 
              total_qaly = sum(discounted_QALY),
              total_cost = sum(discounted_TOTAL.COST),
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
  delta.female.conversion <- cbind(delta.female.conversion, delta.female.conversion[, 2] * wtp.conversion - delta.female.conversion[, 3])
  delta.male.conversion <- male.int.conversion - male.uc.conversion
  delta.male.conversion <- cbind(delta.male.conversion, delta.male.conversion[, 2] * wtp.conversion - delta.male.conversion[, 3])
  
  delta.weighed.conversion <- frac.fem * delta.female.conversion + (1 - frac.fem) * delta.male.conversion
  return(delta.weighed.conversion)
}


# AUTOMATED PROCESSING OF SIM OUTPUT --------------------------------------
#load('All_SMI_vs_UC_basecase.RData')
# stitcher('All_SMI_vs_UC_basecase', '', name.loc2)
sewing_machine(n.runs, paste0('output/', country.id, '/Usual_care_outcomes_run_'))
uc.run.seeds <- sewing.seeds
sewing_machine(n.runs, paste0('output/', country.id, '/All_SMI_vs_UC_basecase_run_'))

output.deltas <- calc.deltas.forD6.1('Any SMI') #calc.deltas('Any SMI') 
output.incidence <- calc.incidences('Any SMI', inc.cohortsize)

# Calculate full output for usual care (first column in result tables D6.1)
f.uc.means <- as.data.frame(sim.results.female.comp[-1])
m.uc.means <- as.data.frame(sim.results.male.comp[-1])
w.uc.means <- frac.fem * f.uc.means + (1 - frac.fem) * m.uc.means
f.uc.add.cost <- calc.additional.output(sim.results.female.comp)
m.uc.add.cost <- calc.additional.output(sim.results.male.comp)
w.uc.add.means <- frac.fem * f.uc.add.cost + (1 - frac.fem) * m.uc.add.cost
all.uc.means <- cbind(w.uc.means, w.uc.add.means)
uc.data.out <-  all.uc.means %>% 
  transmute(life_expecanty = round(mean_life_expectancy, 3),
            total_discounted_qalys = round(mean_total_discounted_qalys, 3),
            total_undiscounted_qalys = round(mean_total_qalys, 3),
            related_medical_costs = round(sum(IHD_cost_disc, MI_cost_disc, CHF_cost_disc, stroke_cost_disc, amput_cost_disc, blind_cost_disc, ulcer_cost_disc, renal_cost_disc, uncomp_cost_disc)),
            unrelated_medical_costs = round(future_medical_cost_disc),
            informal_care_costs = round(informal_care_cost_disc),
            productivity_losses = round(production_losses_disc),
            future_nonmedical_cons = round(future_nonmedical_cost_disc),
            incremental_total_costs = round(mean_total_discounted_costs),
            CHF_incidence = round(mean_CHF_rate * inc.cohortsize, 3), 
            MI_incidence =  round(mean_MI_rate * inc.cohortsize, 3),
            Stroke_incidence =  round(mean_STROKE_rate * inc.cohortsize, 3),
            Blindness_incidence =  round(mean_BLIND_rate * inc.cohortsize, 3),
            ERSD_incidence =  round(mean_RENAL_rate * inc.cohortsize, 3),
            Amputation_incidence =  round(mean_AMP1_rate + mean_AMP2_rate * inc.cohortsize, 3),
            Ulcer_incidence =  round(mean_ULCER_rate * inc.cohortsize, 3)
  )


for (i in 1:13) {
  print(i)
  sewing_machine(n.runs, paste0('output/', country.id, '/Rank', i, '_basecase_run_'))
  if(!identical(sewing.seeds, uc.run.seeds)) stop(paste('Seeds for intervention', i, 'do not match those for usual care'))
  output.rank <- calc.deltas.forD6.1(paste('Rank', i))
  output.deltas <- rbind(output.deltas, output.rank)
  incidence.rank <- calc.incidences(paste('Rank', i), inc.cohortsize)
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
#   incidence.rank <- calc.incidences(paste('Rank', i, 'alt pop'), inc.cohortsize)
#   output.incidence <- rbind(output.incidence, incidence.rank)
# }
# 
# # Add short effect rank 6 analysis
# rm(sim.results.female, sim.results.male, sim.results.female.comp, sim.results.male.comp)
# load('Rank6_shorteff.RData')
# output.deltas <- rbind(output.deltas, calc.deltas('Rank 6 short eff'))

# print(output.deltas)
# print(output.incidence)


# Sink out put to CSV
write.csv(uc.data.out, paste0('output/Usual_care_outcomes_', country.id, '.csv'))
# write.csv(output.deltas, paste0('output/Deltas_all_analyses_', country.id, '.csv'))
# write.csv(output.incidence, paste0('output/Difference_complication_incidence_', country.id, '.csv'))



 
# AUTOMATED CONVERGENCE CHECK ---------------------------------------------
# Load intervention names for titles
int.names <- read.csv('R/Intervention_names.csv', header = FALSE) %>% pull(`V1`)

# Convergence Headroom - choose between files or use the stitcher

# load('All_SMI_vs_UC_basecase.RData')
# stitcher('All_SMI_vs_UC_basecase', '', name.loc2)
sewing_machine(n.runs, paste0('output/', country.id, '/Usual_care_outcomes_run_'))
sewing_machine(n.runs, paste0('output/', country.id, '/All_SMI_vs_UC_basecase_run_'))

allvuc.converge <- calc.convergence()
 # Convergence plot Headroom
png(filename = paste0('output/', country.id, '/All_vs_UC_converge_WTP_', wtp.conversion / 1000, 'K.png'), height = 15, width = 20, units = 'cm', res = 150)
plot(x = 1:nrow(allvuc.converge), y = allvuc.converge[ ,4], type = 'l', xlab = '# Patients (M/F weighed)', ylab = paste('Incremental Headroom (discounted, WTP = ', wtp.conversion/1000, 'K)'), main = 'Overall SMI effect')
dev.off()
  # Convergence plot QALYs
png(filename = paste0('output/', country.id, '/All_vs_UC_converge_QALY.png'), height = 15, width = 20, units = 'cm', res = 150)
plot(x = 1:nrow(allvuc.converge), y = allvuc.converge[ ,2], type = 'l', xlab = '# Patients (M/F weighed)', ylab = 'Incremental QALY (discounted)', main = 'Overall SMI effect')
dev.off()

for (i in 1:13) {
  print(i)
  sewing_machine(n.runs, paste0('output/', country.id, '/Rank', i, '_basecase_run_'))
  converge.rank <- calc.convergence()
  # Plot Headroom
  png(filename = paste0('output/', country.id, '/Rank_', i, '_convergence_WTP_', wtp.conversion / 1000, 'K.png'), height = 15, width = 20, units = 'cm', res = 150)
  plot(x = 1:nrow(converge.rank), y = converge.rank[ ,4], type = 'l', xlab = '# Patients (M/F weighed)', ylab = paste('Incremental Headroom (discounted, WTP = ', wtp.conversion/1000, 'K)'), main = int.names[i])
  dev.off()
  # Plot QALYs
  png(filename = paste0('output/', country.id, '/Rank_', i, '_convergence_QALY.png'), height = 15, width = 20, units = 'cm', res = 150)
  plot(x = 1:nrow(converge.rank), y = converge.rank[ ,2], type = 'l', xlab = '# Patients (M/F weighed)', ylab = paste('Incremental Headroom (discounted, WTP = ', wtp.conversion/1000, 'K)'), main = int.names[i])
  dev.off()
}


