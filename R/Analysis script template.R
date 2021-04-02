# Variable defined to keep track of simulation time (delete afterwards)
init <- Sys.time()

# Model function: 
source("R/SMI in type II diabetes - HE model v3.R")

# Control variables of simulation
npats_input <- 50 # Number of patients in simulation
tx_cost_input <- 0 # Total treatment cost

treateff_start_input   <- 1 # Cycle in which treatment effect starts
treateff_end_input     <- 4 # Cycle in which treatment effect ends
treateff_decline_input <- 2 # Cycle in which treatment effect starts to decline linearly

treateff_hba1c_input <- -1.82 # Treatment effect on HbA1c (in absolute %-points HbA1c)
treateff_hdl_input   <- 0 # Treatment effect on HDL-cholesterol (absolute effect, which unit??)
treateff_ldl_input   <- 0 # Treatment effect on LDL-cholesterol (absolute effect, which unit??)

# The sim.vars object collects all parameters that define the simulation into one object
# Then, it is saved with the output of the simulation. This way, if we  have multiple output files, we always have the information 
# on the relevant input parameters that were used to produce the output.
sim.vars <- list(npats_input, tx_cost_input, mget(apropos('treateff.')))

# A patient-level model make use of patient characteristics. 
baseline_characteristics <- read.csv("input/UK/baseline_characteristics_UK.csv", sep=",")
# baseline_characteristics <- read.csv("/input/UK/baseline_characteristics_UK_rank_1_study_pop.csv", sep=",")

# Direct costs of diabetes-related complications for UK are age-gender dependent.
male_cost_inputs   <- read.csv("input/UK/Event cost male 2020.csv", sep=",")
female_cost_inputs <- read.csv("input/UK/Event cost female 2020.csv", sep=",")

# When a societal perspective is adopted, we also have future costs. This are obtained from the PAID online tool.
future_medical_cost_inputs    <- read.csv("input/UK/UKPAID__AllORUnrelated_Costs_2020.csv", sep=",") # Ingelin: 15/12/2020
future_nonmedical_cost_inputs <- read.csv("input/UK/UK_nonmedical_futurecosts_data.csv", sep=",") # UK costs updated 11/12/2020 - corrected version from Hamraz

# Utilities are UK-based and age/gender dependent
qol_inputs <- read.csv("input/UK/qol_inputs_UK.csv", sep=",")

# Utility decrements associated to diabetes-related events
qol_events_inputs <- read.csv("input/UK/qol_events_inputs_UK.csv", sep=",")

# TODO: also extract output processing scripts to here.

### Run the model below for females and males separately.

##### Females #####

# Intervention

sim_results_female <- SMDMII_model_simulation(npats_input,  #patient_size_input: 
                                              1,  #female_input, 1 = female
                                              tx_cost_input, #tx_cost_input
                                              c(treateff_hba1c_input,treateff_start_input,treateff_end_input,treateff_decline_input), #treatment_effect_HbA1c_input --> from COMPAR + Assumption
                                              c(treateff_hdl_input,treateff_start_input,treateff_end_input,treateff_decline_input),  #treatment_effect_HDL_input 
                                              c(treateff_ldl_input, treateff_start_input,treateff_end_input,treateff_decline_input),  #treatment_effect_LDL_input --> from COMPAR + Assumption
                                              0, #treatment_effect_BMI_input from MH2020
                                              0.035, #cost_disc_rate_input
                                              0.035, #qol_disc_rate_input
                                              65, #retirement_age_input
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
View(sim_CE_results_female_table)



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

View(sim_clinical_results_female_table)


#Comparator


sim_results_female_comp <- SMDMII_model_simulation(npats_input,  #patient_size_input: run 500 for LOLA
                                                   1,  #female_input, 1 = female
                                                   tx_cost_input, #tx_cost_input --> Gimon
                                                   rep(0,4), #treatment_effect_HbA1c_input --> from COMPAR + Assumption
                                                   rep(0,4),  #treatment_effect_HDL_input 
                                                   rep(0,4),  #treatment_effect_LDL_input --> from COMPAR + Assumption                                                   
                                                   0, #treatment_effect_BMI_input from MH2020
                                                   0.035, #cost_disc_rate_input
                                                   0.035, #qol_disc_rate_input
                                                   65, #retirement_age_input
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
View(sim_CE_results_female_table_comp)

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

View(sim_clinical_results_female_table_comp)



##### Males #####

sim_results_male <- SMDMII_model_simulation(npats_input, #patient_size_input: run 500 for LOLA
                                            0, #female_input, 1 = female
                                            tx_cost_input, #tx_cost_input --> Gimon
                                            c(treateff_hba1c_input,treateff_start_input,treateff_end_input,treateff_decline_input), #treatment_effect_HbA1c_input --> from COMPAR + Assumption
                                            c(treateff_hdl_input,treateff_start_input,treateff_end_input,treateff_decline_input),  #treatment_effect_HDL_input 
                                            c(treateff_ldl_input, treateff_start_input,treateff_end_input,treateff_decline_input),  #treatment_effect_LDL_input --> from COMPAR + Assumption                                                   
                                            0, #treatment_effect_BMI_input from MH2020
                                            0.035, #cost_disc_rate_input
                                            0.035, #qol_disc_rate_input
                                            65, #retirement_age_input
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



sim_results_male_comp <- SMDMII_model_simulation(npats_input, #patient_size_input: run 500 for LOLA
                                                 0, #female_input, 1 = female
                                                 tx_cost_input, #tx_cost_input --> Gimon
                                                 rep(0,4), #treatment_effect_HbA1c_input --> from COMPAR + Assumption
                                                 rep(0,4),  #treatment_effect_HDL_input 
                                                 rep(0,4),  #treatment_effect_LDL_input --> from COMPAR + Assumption                                                   
                                                 0, #treatment_effect_BMI_input from MH2020
                                                 0.035, #cost_disc_rate_input
                                                 0.035, #qol_disc_rate_input
                                                 65, # retirement_age_input
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
