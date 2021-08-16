# INPUT PARAMETERS FOR CALCULATIONS ---------------------------------------
# Fraction females in target pop
frac.fem <- 26706/58171

#Lambda (willingness to pay)
wtp <- 20000

# load('output/PSA_Rank1.RData', verbose = TRUE)
# 
# load('output/PSA_Rank2.RData', verbose = TRUE)

load('output/PSA_All_SMI_vs_UC_basecase.RData', verbose = TRUE)
load('output/PSA_Rank11.RData', verbose = TRUE)

#PSA Analysis
w.int.means <- frac.fem * psa.results.female + (1 - frac.fem) * psa.results.male
w.uc.means <- frac.fem * psa.results.female.comp + (1 - frac.fem) * psa.results.male.comp

delta.means <- w.int.means - w.uc.means


psa.out <- data.frame(
  CI_life_expect = quantile(delta.means$mean_life_expectancy, c(0.05, 0.95)),
  CI_total_cost = quantile(delta.means$mean_total_costs, c(0.05, 0.95)),
  CI_qalys = quantile(delta.means$mean_total_qalys, c(0.05, 0.95)),
)