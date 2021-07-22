### CREATE OUTPUTS FOR INTERNAL DISCUSSION JULY 2021


# LOAD SIMULATION DATA ----------------------------------------------------
# load('output/All_SMI_vs_UC_basecase.RData', verbose = TRUE)
# load('output/Usual_care_outcomes.RData', verbose = TRUE)
# load('output/Rank1_basecase.RData', verbose = TRUE)
# load('output/Rank2_basecase.RData', verbose = TRUE)
# load('output/Rank3_basecase.RData', verbose = TRUE)
# load('output/Rank6_basecase.RData', verbose = TRUE)
# load('output/Rank8_basecase.RData', verbose = TRUE)
# load('output/Rank1_spec_tartetpop.RData', verbose = TRUE)
# load('output/Rank6_spec_targetpop.RData', verbose = TRUE)
# load('output/Rank8_shorter_treateff.RData', verbose = TRUE)
# load('output/Rank8_longer_treateff.RData', verbose = TRUE)

load('output/All_SMI_vs_UC_seed77.RData', verbose = TRUE)
# load('output/All_SMI_vs_UC_seed1984.RData', verbose = TRUE)
# load('output/All_SMI_vs_UC_seed15.RData', verbose = TRUE)
# load('output/All_SMI_vs_UC_seed265979.RData', verbose = TRUE)


# ADDITIONAL PARAMETERS ---------------------------------------------------

# Fraction females in target pop
frac.fem <- 1 - 0.566737715306181

#Lambda (willingness to pay)
wtp <- 20000


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
diff.fem <- data.frame(
  Usual_care = t(f.uc.means),
  SMI = t(f.int.means),
  Incr = t(f.int.means) - t(f.uc.means)
)

diff.male <- data.frame(
  Usual_care = t(m.uc.means),
  SMI = t(m.int.means),
  Incr = t(m.int.means) - t(m.uc.means)
)

print(relevant.out)
print(diff.fem)
print(diff.male)


write.csv(diff.fem, 'output/Deltas_female.csv')
write.csv(diff.male, 'output/Deltas_male.csv')