### ANALYSES FOR LOLA MANUSCRIPT ###
library(beepr)

# IMPORT POST TRIAL PATIENT DATA ------------------------------------------

pop.intervention <- read.csv('input/baseline_characteristics_DESMOND_intervention.csv')

pop.control <- read.csv('input/baseline_characteristics_DESMOND_control.csv')

int.cost <- 76*1.25*1.1679

# ADAPT INPUT DATA --------------------------------------------------------
# Import CPI data (source: CBS, through FunCpi form ErasmuscorneR)
inflators <- read.csv('input/CPI_NL_2010_2019.csv')

# Inflate Mt Hood costs inputs from 18 to 19
cost_inputs <- cost_inputs * inflators[9, 2]

# Inflate PAID estimates from 17 to 19
future_medical_cost_inputs[, -1] <- future_medical_cost_inputs[, -1] * inflators[8, 2]
future_nonmedical_cost_inputs[, -1] <- future_nonmedical_cost_inputs[, -1] * inflators[8, 2]

# RUN SIMULATIONS ---------------------------------------------------------

sim.start <- Sys.time()

# Run simulation for the population who received intervention
baseline_characteristics <- pop.intervention

sim.int.male <- SMDMII_model_simulation(1000, 0, 0, 0, 0, 0, 0.04, 0.015, 0, 77)
sim.int.female <- SMDMII_model_simulation(1000, 1, 0, 0, 0, 0, 0.04, 0.015, 0, 77)

# Run simulation for the control population
baseline_characteristics <- pop.control

sim.control.male <- SMDMII_model_simulation(1000, 0, 0, 0, 0, 0, 0.04, 0.015, 0, 77)
sim.control.female <- SMDMII_model_simulation(1000, 1, 0, 0, 0, 0, 0.04, 0.015, 0, 77)

sim.finish <- Sys.time()

print(sim.finish - sim.start)
beep('mario')


# PRESENT RESULTS ---------------------------------------------------------
frac.male <- 0.54 # Percentage males in DESMOND

# Results tables
res.tab.name <- c("Life expectancy", "Total QALYs", "Direct healthcare costs", "Informal care costs", "Productivity costs", "Future medical costs", "Future non-medical consumption", "Total costs")

sim.int.male.tab <- data.frame(sim.int.male$mean_life_expectancy,
                               sim.int.male$mean_total_qalys,
                               sim.int.male$mean_complication_costs + sim.int.male$mean_nocomp_costs + int.cost,
                               sim.int.male$mean_inf_care_costs,
                               sim.int.male$mean_prod_loss_costs, 
                               sim.int.male$mean_future_medical_costs,
                               sim.int.male$mean_future_nonmedical_costs,
                               sim.int.male$mean_total_costs)

names(sim.int.male.tab) <- res.tab.name

sim.int.female.tab <- data.frame(sim.int.female$mean_life_expectancy,
                               sim.int.female$mean_total_qalys,
                               sim.int.female$mean_complication_costs + sim.int.female$mean_nocomp_costs + int.cost,
                               sim.int.female$mean_inf_care_costs,
                               sim.int.female$mean_prod_loss_costs, 
                               sim.int.female$mean_future_medical_costs,
                               sim.int.female$mean_future_nonmedical_costs,
                               sim.int.female$mean_total_costs)

names(sim.int.female.tab) <- res.tab.name

sim.control.male.tab <- data.frame(sim.control.male$mean_life_expectancy,
                               sim.control.male$mean_total_qalys,
                               sim.control.male$mean_complication_costs + sim.control.male$mean_nocomp_costs,
                               sim.control.male$mean_inf_care_costs,
                               sim.control.male$mean_prod_loss_costs, 
                               sim.control.male$mean_future_medical_costs,
                               sim.control.male$mean_future_nonmedical_costs,
                               sim.control.male$mean_total_costs)

names(sim.control.male.tab) <- res.tab.name

sim.control.female.tab <- data.frame(sim.control.female$mean_life_expectancy,
                                 sim.control.female$mean_total_qalys,
                                 sim.control.female$mean_complication_costs + sim.control.female$mean_nocomp_costs,
                                 sim.control.female$mean_inf_care_costs,
                                 sim.control.female$mean_prod_loss_costs, 
                                 sim.control.female$mean_future_medical_costs,
                                 sim.control.female$mean_future_nonmedical_costs,
                                 sim.control.female$mean_total_costs)

names(sim.control.female.tab) <- res.tab.name

# Calculate weighted average m/f ratio
sim.int.weighed.tab <- (frac.male*sim.int.male.tab + (1 - frac.male) * sim.int.female.tab)

sim.control.weighed.tab <- (frac.male*sim.control.male.tab + (1 - frac.male) * sim.control.female.tab)

sim.results <- rbind(sim.control.weighed.tab, sim.int.weighed.tab, sim.int.weighed.tab - sim.control.weighed.tab)
row.names(sim.results) <- c('Control', 'Intervention', 'Differnce')

cost.perspectives <- data.frame('Direct healthcare' = sim.results[, 3],
                                'Direct societal' = rowSums(sim.results[, 3:5]),
                                'Total healthcare' = rowSums(sim.results[, c(3, 6)]),
                                'Total societal' = rowSums(sim.results[, 3:7])
                                )
