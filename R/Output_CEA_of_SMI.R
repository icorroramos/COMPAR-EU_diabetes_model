### OUTPUT OF ANALYSES OF 13 RANKED INTERVENTIONS ###


# INPUT PARAMETERS FOR CALCULATIONS ---------------------------------------
# Fraction females in target pop
frac.fem <- 26706/58171

#Lambda (willingness to pay)
wtp <- 20000


# FUNCTION WITH CALCULATIONS ----------------------------------------------
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
              life_expecanty = round(mean_life_expectancy, 3),
              total_costs = round(mean_total_costs),
              total_qalys = round(mean_total_qalys, 3),
              headroom = round(wtp * mean_total_qalys - mean_total_costs, 2))
}


# AUTOMATED PROCESSING OF SIM OUTPUT --------------------------------------
load('output/All_SMI_vs_UC_basecase.RData')

output.deltas <- calc.deltas('Any SMI')

for (i in 1:13) {
  print(i)
  load('output/Usual_care_outcomes.RData')
  load(paste0('output/Rank', i, '_basecase.RData'))
  output.rank <- calc.deltas(paste('Rank', i))
  output.deltas <- rbind(output.deltas, output.rank)
}

# Remove all sim results to prevent use of base case comparator in alternate analyses
rm(sim.results.female, sim.results.male, sim.results.female.comp, sim.results.male.comp)

for (i in c(1, 2, 6)) {
  print(i)
  load(paste0('output/Rank', i, '_spec_targetpop.RData'))
  output.rank <- calc.deltas(paste('Rank', i, 'alt pop'))
  output.deltas <- rbind(output.deltas, output.rank)
}

# Add short effect rank 6 analysis
rm(sim.results.female, sim.results.male, sim.results.female.comp, sim.results.male.comp)
load('output/Rank6_shorteff.RData')
output.deltas <- rbind(output.deltas, calc.deltas('Rank 6 short eff'))

print(output.deltas)

# Sink out put to CSV
write.csv(output.deltas, 'output/Deltas_all_analyses.csv')




# MANUAL ANALYSIS ---------------------------------------------------------

# LOAD SIMULATION DATA ----------------------------------------------------
# load('output/All_SMI_vs_UC_basecase.RData', verbose = TRUE)
# load('output/Usual_care_outcomes.RData', verbose = TRUE)
 load('output/Rank1_basecase.RData', verbose = TRUE)
# load('output/Rank2_basecase.RData', verbose = TRUE)
# load('output/Rank3_basecase.RData', verbose = TRUE)
# load('output/Rank4_basecase.RData', verbose = TRUE)
# load('output/Rank5_basecase.RData', verbose = TRUE)
# load('output/Rank6_basecase.RData', verbose = TRUE)
# load('output/Rank7_basecase.RData', verbose = TRUE)
# load('output/Rank8_basecase.RData', verbose = TRUE)
# load('output/Rank9_basecase.RData', verbose = TRUE)
# load('output/Rank10_basecase.RData', verbose = TRUE)
# load('output/Rank11_basecase.RData', verbose = TRUE)
# load('output/Rank12_basecase.RData', verbose = TRUE)
# load('output/Rank13_basecase.RData', verbose = TRUE)

# load('output/Rank1_spec_targetpop.RData', verbose = TRUE)
# load('output/Rank6_spec_targetpop.RData', verbose = TRUE)

# load('output/Rank8_shorter_treateff.RData', verbose = TRUE)
# load('output/Rank8_longer_treateff.RData', verbose = TRUE)

# load('output/All_SMI_vs_UC_seed77.RData', verbose = TRUE)
# load('output/All_SMI_vs_UC_seed1984.RData', verbose = TRUE)
# load('output/All_SMI_vs_UC_seed15.RData', verbose = TRUE)
# load('output/All_SMI_vs_UC_seed265979.RData', verbose = TRUE)

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
# write.csv(diff.fem, 'output/Deltas_female.csv')
# write.csv(diff.male, 'output/Deltas_male.csv')