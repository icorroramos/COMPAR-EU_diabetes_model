### ANALYSES FOR LOLA MANUSCRIPT ###


# IMPORT POST TRIAL PATIENT DATA ------------------------------------------

pop.intervention <- read.csv('input/baseline_characteristics_DESMOND_intervention.csv')

pop.control <- read.csv('input/baseline_characteristics_DESMOND_control.csv')


# ADAPT INPUT DATA --------------------------------------------------------
# Import CPI data (source: CBS, through FunCpi form ErasmuscorneR)
inflators <- read.csv('input/CPI_NL_2010_2019.csv')

# Inflate Mt Hood costs inputs from 18 to 19
cost_inputs <- cost_inputs * inflators[9, 2]

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



# PRESENT RESULTS ---------------------------------------------------------

res.tab.name <- c("Complication costs", "No complication costs", "Informal care costs", "Productivity costs", "Total costs", "Total QALYs")

sim.int.male.tab <- data.frame(sim.int.male$mean_complication_costs)