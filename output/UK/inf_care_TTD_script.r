# setwd to github folder

# load file
load("Output EuHEA scenario 20k pats.RData")

female_int <- sim.results.female$simulation_patients_history
female_int <- female_int[which(female_int$SDURATION>0),]

female_com <- sim.results.female.comp$simulation_patients_history
female_com <- female_com[which(female_com$SDURATION>0),]

male_int <- sim.results.male$simulation_patients_history
male_int <- male_int[which(male_int$SDURATION>0),]

male_com <- sim.results.male.comp$simulation_patients_history
male_com <- male_com[which(male_com$SDURATION>0),]

# calculate TTD
female_int$TTD <- unlist(aggregate(female_int$SDURATION, list(id = female_int$SIMID), rev)$x)
female_com$TTD <- unlist(aggregate(female_com$SDURATION, list(id = female_com$SIMID), rev)$x)

male_int$TTD <- unlist(aggregate(male_int$SDURATION, list(id = male_int$SIMID), rev)$x)
male_com$TTD <- unlist(aggregate(male_com$SDURATION, list(id = male_com$SIMID), rev)$x)
#fem_1_pat_1_2$TTD <- unlist(aggregate(fem_1_pat_1_2$SDURATION, list(id = fem_1_pat_1_2$SIMID), rev)$x)

# Probability weekly informal care: females
p_week_female_int <- exp(-1.451 + 0.368 + 0.054*(female_int$CURR.AGE-70)+0.0003*(female_int$CURR.AGE-70)^2-0.061*female_int$TTD)/(1+exp(-1.451 + 0.368 + 0.054*(female_int$CURR.AGE-70) + 0.0003*(female_int$CURR.AGE-70)^2 -0.061*female_int$TTD))
p_week_female_com <- exp(-1.451 + 0.368 + 0.054*(female_com$CURR.AGE-70)+0.0003*(female_com$CURR.AGE-70)^2-0.061*female_com$TTD)/(1+exp(-1.451 + 0.368 + 0.054*(female_com$CURR.AGE-70) + 0.0003*(female_com$CURR.AGE-70)^2 -0.061*female_com$TTD))

# Probability weekly informal care: males
p_week_male_int <- exp(-1.451 + 0 + 0.054*(male_int$CURR.AGE-70)+0.0003*(male_int$CURR.AGE-70)^2-0.061*male_int$TTD)/(1+exp(-1.451 + 0 + 0.054*(male_int$CURR.AGE-70) + 0.0003*(male_int$CURR.AGE-70)^2 -0.061*male_int$TTD))
p_week_male_com <- exp(-1.451 + 0 + 0.054*(male_com$CURR.AGE-70)+0.0003*(male_com$CURR.AGE-70)^2-0.061*male_com$TTD)/(1+exp(-1.451 + 0 + 0.054*(male_com$CURR.AGE-70) + 0.0003*(male_com$CURR.AGE-70)^2 -0.061*male_com$TTD))



# Hours per day -- > *7 then hours per week
h_week_female_int <- 7*exp(0.497 + 0.112 + (female_int$CURR.AGE-70)*0.019-0.034*female_int$TTD)
h_week_female_com <- 7*exp(0.497 + 0.112 + (female_com$CURR.AGE-70)*0.019-0.034*female_com$TTD)

h_week_male_int <- 7*exp(0.497 + 0 + (male_int$CURR.AGE-70)*0.019-0.034*male_int$TTD)
h_week_male_com <- 7*exp(0.497 + 0 + (male_com$CURR.AGE-70)*0.019-0.034*male_com$TTD)


# (Average hours per week * cost hour)*number of weeks year: use UK specific data
female_int$INF.CARE.COST.TTD <- p_week_female_int*h_week_female_int*22.79*52
female_com$INF.CARE.COST.TTD <- p_week_female_com*h_week_female_com*22.79*52

male_int$INF.CARE.COST.TTD <- p_week_male_int*h_week_male_int*22.79*52
male_com$INF.CARE.COST.TTD <- p_week_male_com*h_week_male_com*22.79*52

# discount

discount_factor_fem_int <- sim.results.female$simulation_patients_history[which(sim.results.female$simulation_patients_history$SDURATION>0),]$cost_discount_factor
discount_factor_mal_int <- sim.results.male$simulation_patients_history[which(sim.results.male$simulation_patients_history$SDURATION>0),]$cost_discount_factor

discount_factor_fem_com <- sim.results.female.comp$simulation_patients_history[which(sim.results.female.comp$simulation_patients_history$SDURATION>0),]$cost_discount_factor
discount_factor_mal_com <- sim.results.male.comp$simulation_patients_history[which(sim.results.male.comp$simulation_patients_history$SDURATION>0),]$cost_discount_factor


# aggregate per patient and averaged
#mean(aggregate(female_int$INF.CARE.COST/discount_factor_fem, list(id = female_int$SIMID), sum)$x) #previous approach: sim.results.female$mean_inf_care_costs
female_int_infcare <- mean(aggregate(female_int$INF.CARE.COST.TTD/discount_factor_fem_int, list(id = female_int$SIMID), sum)$x)

#mean(aggregate(female_com$INF.CARE.COST, list(id = female_com$SIMID), sum)$x) #previous approach
female_com_infcare <- mean(aggregate(female_com$INF.CARE.COST.TTD/discount_factor_fem_com, list(id = female_com$SIMID), sum)$x)

#mean(aggregate(male_int$INF.CARE.COST, list(id = male_int$SIMID), sum)$x) #previous approach: sim.results.male$mean_inf_care_costs
male_int_infcare <- mean(aggregate(male_int$INF.CARE.COST.TTD/discount_factor_mal_int, list(id = male_int$SIMID), sum)$x)

#mean(aggregate(male_com$INF.CARE.COST, list(id = male_com$SIMID), sum)$x) #previous approach
male_com_infcare <- mean(aggregate(male_com$INF.CARE.COST.TTD/discount_factor_mal_com, list(id = male_com$SIMID), sum)$x)


# Weighted mean based on proportion of males and females (UK case study values)
frac.fem <- 26706/58171 # Fraction of females in UK DM2 population


weighted_int_infcare <- frac.fem * female_int_infcare + (1 - frac.fem) * male_int_infcare
weighted_com_infcare <- frac.fem * female_com_infcare + (1 - frac.fem) * male_com_infcare
delta_infcare <- weighted_int_infcare - weighted_com_infcare
