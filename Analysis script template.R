

# Variable defined to keep track of simulation time (delete afterwards)
init <- Sys.time()

# Control variables of simulation
npats <- 500 # Number of patients in simulation
tx.cost <- 0 # Total treatment cost

treateff.start <- 1 # Cycle in which treatment effect starts
treateff.end <- 4 # Cycle in which treatment effect ends
treateff.decline <- 2 # Cycle in which treatment effect strats to decline linearly

treateff.hba1c <- -1.82 # Treatment effect on HbA1c (in absolute %-points HbA1c)
treateff.hdl <- 0 # Treatment effect on HDL-choleresterol (absolute effect, which unit??)
treateff.ldl <- 0 # Treatment effect on LDL-choleresterol (absolute effect, which unit??)

sim.vars <- list(npats, tx.cost, mget(apropos('treateff.')))



# A patient-level model make use of patient characteristics. 
# QUESTION: still some issues, ask for units.
baseline_characteristics <- read.csv("input/baseline_characteristics_UK.csv", sep=",")
# baseline_characteristics <- read.csv("input/baseline_characteristics_UK_rank_1_study_pop.csv", sep=",")

# Direct costs of diabetes-related complications for UK are age-gender dependent.
male_cost_inputs   <- read.csv("input/Event cost male 2020.csv", sep=",")
female_cost_inputs <- read.csv("input/Event cost female 2020.csv", sep=",")

# When a societal perspective is adopted, we also have future costs. This are obtained from the PAID online tool.
future_medical_cost_inputs <- read.csv("input/UKPAID__AllORUnrelated_Costs_2020.csv", sep=",") # Ingelin: 15/12/2020
future_nonmedical_cost_inputs <- read.csv("input/UK_nonmedical_futurecosts_data.csv", sep=",") # UK costs updated 11/12/2020 - corrected version from Hamraz

# Utilities are UK-based and age/gender dependent
qol_inputs <- read.csv("input/qol_inputs_UK.csv", sep=",")


# TODO: also extract output processing scripts to here.



### Run the model below for females and males separately.


##### Females #####

# Intervention

sim_results_female <- SMDMII_model_simulation(npats,  #patient_size_input: 
                                              1,  #female_input, 1 = female
                                              tx.cost, #tx_cost_input --> Gimon
                                              c(treateff.hba1c,treateff.start,treateff.end,treateff.decline), #treatment_effect_HbA1c_input --> from COMPAR + Assumption
                                              c(treateff.hdl,treateff.start,treateff.end,treateff.decline),  #treatment_effect_HDL_input 
                                              c(treateff.ldl, treateff.start,treateff.end,treateff.decline),  #treatment_effect_LDL_input --> from COMPAR + Assumption
                                              0, #treatment_effect_BMI_input from MH2020
                                              0.035, #cost_disc_rate_input
                                              0.035, #qol_disc_rate_input
                                              0, #run_PSA_input, 0 == no PSA
                                              77 #seed_input
)

# Results tables
sim_CE_results_female_table <- matrix(c(sim_results_female$mean_complication_costs,sim_results_female$mean_nocomp_costs, sim_results_female$mean_tx_costs,
                                        sim_results_female$mean_inf_care_costs, sim_results_female$mean_prod_loss_costs,
                                        sim_results_female$mean_future_medical_costs, sim_results_female$mean_future_nonmedical_costs,
                                        sim_results_female$mean_total_costs, sim_results_female$mean_total_qalys), nrow = 1)

colnames(sim_CE_results_female_table) <- c("Complication costs", "No complication costs", "Tx costs","Informal care costs", "Productivity costs",
                                           "Future medical costs", "Future non-medical costs", "Total costs", "Total QALYs")
rownames(sim_CE_results_female_table) <- "Intervention"
#View(sim_CE_results_female_table)

sim_clinical_results_female_table <- matrix(c(sim_results_female$mean_life_expectancy,
                                              sim_results_female$mean_CHF_rate,
                                              sim_results_female$mean_MI_rate, 
                                              sim_results_female$mean_BLIND_rate,
                                              sim_results_female$mean_ULCER_rate, 
                                              sim_results_female$mean_AMP1_rate,
                                              sim_results_female$mean_AMP2_rate, 
                                              sim_results_female$mean_RENAL_rate,
                                              sim_results_female$mean_STROKE_rate), nrow = 1)

colnames(sim_clinical_results_female_table) <- c("Life expectancy", "CHF rate", "MI rate", "Blindness rate", "Ulcer rate", 
                                                 "1st amputation rate", "2nd amputation rate", "Renal failure rate", "Stroke rate")
rownames(sim_clinical_results_female_table) <- "Intervention"

#View(sim_clinical_results_female_table)


#Comparator


sim_results_female_comp <- SMDMII_model_simulation(npats,  #patient_size_input: run 500 for LOLA
                                                   1,  #female_input, 1 = female
                                                   tx.cost, #tx_cost_input --> Gimon
                                                   rep(0,4), #treatment_effect_HbA1c_input --> from COMPAR + Assumption
                                                   rep(0,4),  #treatment_effect_HDL_input 
                                                   rep(0,4),  #treatment_effect_LDL_input --> from COMPAR + Assumption                                                   0, #treatment_effect_BMI_input from MH2020
                                                   0, #treatment_effect_BMI_input from MH2020
                                                   0.035, #cost_disc_rate_input
                                                   0.035, #qol_disc_rate_input
                                                   0, #run_PSA_input, 0 == no PSA
                                                   77 #seed_input
)

# Results tables
sim_CE_results_female_table_comp <- matrix(c(sim_results_female_comp$mean_complication_costs,
                                             sim_results_female_comp$mean_nocomp_costs, 
                                             sim_results_female_comp$mean_tx_costs,
                                             sim_results_female_comp$mean_inf_care_costs, 
                                             sim_results_female_comp$mean_prod_loss_costs,
                                             sim_results_female_comp$mean_future_medical_costs, 
                                             sim_results_female_comp$mean_future_nonmedical_costs,
                                             sim_results_female_comp$mean_total_costs, 
                                             sim_results_female_comp$mean_total_qalys), nrow = 1)

colnames(sim_CE_results_female_table_comp) <- c("Complication costs", "No complication costs", "Tx costs","Informal care costs", "Productivity costs",
                                                "Future medical costs", "Future non-medical costs", "Total costs", "Total QALYs")
rownames(sim_CE_results_female_table_comp) <- "Comparator"
#View(sim_CE_results_female_table_comp)

sim_clinical_results_female_table_comp <- matrix(c(sim_results_female_comp$mean_life_expectancy,
                                                   sim_results_female_comp$mean_CHF_rate,
                                                   sim_results_female_comp$mean_MI_rate, 
                                                   sim_results_female_comp$mean_BLIND_rate,
                                                   sim_results_female_comp$mean_ULCER_rate, 
                                                   sim_results_female_comp$mean_AMP1_rate,
                                                   sim_results_female_comp$mean_AMP2_rate, 
                                                   sim_results_female_comp$mean_RENAL_rate,
                                                   sim_results_female_comp$mean_STROKE_rate), nrow = 1)

colnames(sim_clinical_results_female_table_comp) <- c("Life expectancy", "CHF rate", "MI rate", "Blindness rate", "Ulcer rate", 
                                                      "1st amputation rate", "2nd amputation rate", "Renal failure rate", "Stroke rate")
rownames(sim_clinical_results_female_table_comp) <- "Comparator"

#View(sim_clinical_results_female_table_comp)



##### Males #####

sim_results_male <- SMDMII_model_simulation(npats, #patient_size_input: run 500 for LOLA
                                            0, #female_input, 1 = female
                                            tx.cost, #tx_cost_input --> Gimon
                                            c(treateff.hba1c,treateff.start,treateff.end,treateff.decline), #treatment_effect_HbA1c_input --> from COMPAR + Assumption
                                            c(treateff.hdl,treateff.start,treateff.end,treateff.decline),  #treatment_effect_HDL_input 
                                            c(treateff.ldl, treateff.start,treateff.end,treateff.decline),  #treatment_effect_LDL_input --> from COMPAR + Assumption                                                   0, #treatment_effect_BMI_input from MH2020
                                            0, #treatment_effect_BMI_input from MH2020
                                            0.035, #cost_disc_rate_input
                                            0.035, #qol_disc_rate_input
                                            0, #run_PSA_input, 0 == no PSA
                                            77 #seed_input
)

# Results tables
sim_CE_results_male_table <- matrix(c(sim_results_male$mean_complication_costs,sim_results_male$mean_nocomp_costs, sim_results_male$mean_tx_costs,
                                      sim_results_male$mean_inf_care_costs, sim_results_male$mean_prod_loss_costs,
                                      sim_results_male$mean_future_medical_costs, sim_results_male$mean_future_nonmedical_costs,
                                      sim_results_male$mean_total_costs, sim_results_male$mean_total_qalys), nrow = 1)

colnames(sim_CE_results_male_table) <- c("Complication costs", "No complication costs", "Tx costs","Informal care costs", "Productivity costs",
                                         "Future medical costs", "Future non-medical costs", "Total costs", "Total QALYs")
rownames(sim_CE_results_male_table) <- "Intervention"
#View(sim_CE_results_male_table)

sim_clinical_results_male_table <- matrix(c(sim_results_male$mean_life_expectancy,
                                            sim_results_male$mean_CHF_rate,
                                            sim_results_male$mean_MI_rate, 
                                            sim_results_male$mean_BLIND_rate,
                                            sim_results_male$mean_ULCER_rate, 
                                            sim_results_male$mean_AMP1_rate,
                                            sim_results_male$mean_AMP2_rate, 
                                            sim_results_male$mean_RENAL_rate,
                                            sim_results_male$mean_STROKE_rate), nrow = 1)

colnames(sim_clinical_results_male_table) <- c("Life expectancy", "CHF rate", "MI rate", "Blindness rate", "Ulcer rate", 
                                               "1st amputation rate", "2nd amputation rate", "Renal failure rate", "Stroke rate")
rownames(sim_clinical_results_male_table) <- "Intervention"

# View(sim_clinical_results_male_table)

# View(sim_results_male$simulation_patients_history[,c("HbA1c","LDL","HDL","SDURATION")])



sim_results_male_comp <- SMDMII_model_simulation(npats, #patient_size_input: run 500 for LOLA
                                                 0, #female_input, 1 = female
                                                 tx.cost, #tx_cost_input --> Gimon
                                                 rep(0,4), #treatment_effect_HbA1c_input --> from COMPAR + Assumption
                                                 rep(0,4),  #treatment_effect_HDL_input 
                                                 rep(0,4),  #treatment_effect_LDL_input --> from COMPAR + Assumption                                                   0, #treatment_effect_BMI_input from MH2020
                                                 0, #treatment_effect_BMI_input from MH2020
                                                 0.035, #cost_disc_rate_input
                                                 0.035, #qol_disc_rate_input
                                                 0, #run_PSA_input, 0 == no PSA
                                                 77 #seed_input
)


# Results tables
sim_CE_results_male_table_comp <- matrix(c(sim_results_male_comp$mean_complication_costs,
                                           sim_results_male_comp$mean_nocomp_costs, 
                                           sim_results_male_comp$mean_tx_costs,
                                           sim_results_male_comp$mean_inf_care_costs, 
                                           sim_results_male_comp$mean_prod_loss_costs,
                                           sim_results_male_comp$mean_future_medical_costs, 
                                           sim_results_male_comp$mean_future_nonmedical_costs,
                                           sim_results_male_comp$mean_total_costs, 
                                           sim_results_male_comp$mean_total_qalys), nrow = 1)

colnames(sim_CE_results_male_table_comp) <- c("Complication costs", "No complication costs", "Tx costs","Informal care costs", "Productivity costs",
                                              "Future medical costs", "Future non-medical costs", "Total costs", "Total QALYs")
rownames(sim_CE_results_male_table_comp) <- "Comparator"
# View(sim_CE_results_male_table_comp)

sim_clinical_results_male_table_comp <- matrix(c(sim_results_male_comp$mean_life_expectancy,
                                                 sim_results_male_comp$mean_CHF_rate,
                                                 sim_results_male_comp$mean_MI_rate, 
                                                 sim_results_male_comp$mean_BLIND_rate,
                                                 sim_results_male_comp$mean_ULCER_rate, 
                                                 sim_results_male_comp$mean_AMP1_rate,
                                                 sim_results_male_comp$mean_AMP2_rate, 
                                                 sim_results_male_comp$mean_RENAL_rate,
                                                 sim_results_male_comp$mean_STROKE_rate), nrow = 1)

colnames(sim_clinical_results_male_table_comp) <- c("Life expectancy", "CHF rate", "MI rate", "Blindness rate", "Ulcer rate", 
                                                    "1st amputation rate", "2nd amputation rate", "Renal failure rate", "Stroke rate")
rownames(sim_clinical_results_male_table_comp) <- "Comaprator"

# View(sim_clinical_results_male_table_comp)


# Variable defined to keep track of simulation time (delete afterwards)
end <- Sys.time()
end - init


# SINK RESULTS DATA -------------------------------------------------------
save(sim.vars,
     sim_CE_results_female_table,
     sim_clinical_results_female_table,
     sim_CE_results_female_table_comp,
     sim_clinical_results_female_table_comp,
     sim_CE_results_male_table,
     sim_clinical_results_male_table,
     sim_CE_results_male_table_comp,
     sim_clinical_results_male_table_comp,
     file = 'Sim_results.Rdata'
)

# Optional - save to csv's
# write.csv(sim_CE_results_female_table, 'sim_CE_female_int.csv', quote = FALSE)
# write.csv(sim_clinical_results_female_table, 'sim_clin_female_int.csv', quote = FALSE)
# write.csv(sim_CE_results_female_table_comp, 'sim_CE_female_comp.csv', quote = FALSE)
# write.csv(sim_clinical_results_female_table_comp, 'sim_clin_female_comp.csv', quote = FALSE)
# 
# write.csv(sim_CE_results_male_table,'sim_CE_male_int.csv', quote = FALSE)
# write.csv(sim_clinical_results_male_table, 'sim_clin_male_int.csv', quote = FALSE)
# write.csv(sim_CE_results_male_table_comp, 'sim_CE_male_comp.csv', quote = FALSE)
# write.csv(sim_clinical_results_male_table_comp, 'sim_clin_male_comp.csv', quote = FALSE)