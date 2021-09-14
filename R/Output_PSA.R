# INPUT PARAMETERS FOR CALCULATIONS ---------------------------------------
# Fraction females in target pop
frac.fem <- 26706/58171

#Lambda (willingness to pay)
wtp <- 20000


# FUNCTIONS FOR AUTOMATED OUTPUT PROCESSING -------------------------------

calc.psaout <- function(comp.name) {
  w.int.means <- frac.fem * psa.results.female + (1 - frac.fem) * psa.results.male
  w.uc.means <- frac.fem * psa.results.female.comp + (1 - frac.fem) * psa.results.male.comp
  
  delta.means <- w.int.means - w.uc.means
  delta.means$headroom <- wtp * delta.means$mean_total_qalys - delta.means$mean_total_costs
  
  relevant.output <- data.frame(
    comp = comp.name,
    life_exp_mean = round(mean(delta.means$mean_life_expectancy), 3),
    life_exp_p2.5 = round(quantile(delta.means$mean_life_expectancy, 0.025), 3),
    life_exp_p97.5 = round(quantile(delta.means$mean_life_expectancy, 0.975), 3),
    downstr_costs_mean = round(mean(delta.means$mean_total_costs)),
    downstr_costs_p2.5 = round(quantile(delta.means$mean_total_costs, 0.025)),
    downstr_costs_p97.5 = round(quantile(delta.means$mean_total_costs, 0.975)),
    total_qalys_mean = round(mean(delta.means$mean_total_qalys), 3),
    total_qalys_p2.5 = round(quantile(delta.means$mean_total_qalys, 0.025), 3),
    total_qalys_p97.5 = round(quantile(delta.means$mean_total_qalys, 0.975), 3),
    headroom_mean = round(mean(delta.means$headroom)),
    headroom_p2.5 = round(quantile(delta.means$headroom, 0.025)),
    headroom_p97.5 = round(quantile(delta.means$headroom, 0.975))
  )
  rownames(relevant.output) <- NULL
  return(relevant.output)
}


# AUTOMATED PROCESSING OF PSA OUTPUT --------------------------------------
load('output/PSA_All_SMI_vs_UC_basecase.RData', verbose = TRUE)
output.psa <- calc.psaout('Any SMI')

for (i in c(1:13)) {
  print(i)
  load('output/PSA_Usual_care.RData')
  load(paste0('output/PSA_Rank', i, '.RData'))
  output.rank <- calc.psaout(paste('Rank', i))
  output.psa <- rbind(output.psa, output.rank)
}

print(output.psa)

write.csv(output.psa, file = 'output/PSA_results_all_ranks.csv', row.names = FALSE)





# MANUAL PSA OUTPUT ANALYSIS ----------------------------------------------

load('output/PSA_All_SMI_vs_UC_basecase.RData', verbose = TRUE)

load('output/PSA_Usual_care.RData', verbose = TRUE)

load('output/PSA_Rank1.RData', verbose = TRUE)

#PSA Analysis
w.int.means <- frac.fem * psa.results.female + (1 - frac.fem) * psa.results.male
w.uc.means <- frac.fem * psa.results.female.comp + (1 - frac.fem) * psa.results.male.comp

delta.means <- w.int.means - w.uc.means
headrooms <- wtp * delta.means$mean_total_qalys - delta.means$mean_total_costs

psa.mean <- data.frame(
  mean_life_expect = mean(delta.means$mean_life_expectancy),
  mean_total_cost = mean(delta.means$mean_total_costs),
  mean_qalys = mean(delta.means$mean_total_qalys),
  mean_headroom = mean(headrooms)
)

psa.95ci <- data.frame(
  CI_life_expect = quantile(delta.means$mean_life_expectancy, c(0.05, 0.95)),
  CI_total_cost = quantile(delta.means$mean_total_costs, c(0.05, 0.95)),
  CI_qalys = quantile(delta.means$mean_total_qalys, c(0.05, 0.95)),
  CI_headroom = quantile(headrooms, c(0.05, 0.95))
)

