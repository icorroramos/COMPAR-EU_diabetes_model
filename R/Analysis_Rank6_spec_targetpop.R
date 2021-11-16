##################
# Set-up options #
##################

# Clear global environment, flush memory, load packages
source('R/Analysis_start_chunk.R')

# Variable defined to keep track of simulation time (delete afterwards)
init <- Sys.time()

# Adjust target population
baseline_characteristics <- read.csv("input/UK/baseline_characteristics_UK_rank_6_study_pop.csv", sep=",")

# Treatment effect inputs
treateff_start   <- 1 # Cycle in which treatment effect starts
treateff_end     <- 3 # Cycle in which treatment effect ends
treateff_decline <- 2 # Cycle in which treatment effect starts to decline (linearly)

treateff_hba1c <- -0.7037 # Treatment effect on HbA1c (in absolute %-points HbA1c)
treateff_hdl   <- 3.867 * 0.02586 * 10 # TRANSFORMATION *10 FOR MODEL INPUT # Treatment effect on HDL-cholesterol (absolute effect in mmol/l)
treateff_ldl   <- -6.1872 * 0.02586 * 10 # TRANSFORMATION *10 FOR MODEL INPUT# Treatment effect on LDL-cholesterol (absolute effect in mmol/l)
treateff_bmi <- -0.12 # Treatment effect on BMI (in absolute points)
treateff_sbp <- -1.7 / 10 # TRANSFORMATION /10 FOR MODEL INPUT# Treatment effect on SBP (in absolute mmHg)

# Tx effects are vectors: the current assumption is that the same start, end and decline is assumed for all effect modifiers
treatment_effect_HbA1c_input <- c(treateff_hba1c, treateff_start, treateff_end, treateff_decline)
treatment_effect_HDL_input <- c(treateff_hdl, treateff_start, treateff_end, treateff_decline)
treatment_effect_LDL_input <- c(treateff_ldl, treateff_start, treateff_end, treateff_decline)
treatment_effect_BMI_input <- c(treateff_bmi, treateff_start, treateff_end, treateff_decline)
treatment_effect_SBP_input <- c(treateff_sbp, treateff_start, treateff_end, treateff_decline)


# The sim.vars object collects all parameters that define the simulation into one object
# Then, it is saved with the output of the simulation. This way, if we  have multiple output files, we always have the information 
# on the relevant input parameters that were used to produce the output.
sim.vars <- list(seed_input, npats_input, tx_cost_input, mget(apropos('treateff.')))



# RESULT OUTPUT FOR OUTPUT SCRIPT -----------------------------------------
sim.results.female <- SMDMII_model_simulation(patient_size_input = npats_input,
                                              female_input = 1,
                                              tx_cost_input = tx_cost_input,
                                              treatment_effect_HbA1c_input = treatment_effect_HbA1c_input, 
                                              treatment_effect_HDL_input = reatment_effect_HDL_input,
                                              treatment_effect_LDL_input = treatment_effect_LDL_input, 
                                              treatment_effect_BMI_input = treatment_effect_BMI_input,
                                              treatment_effect_SBP_input = treatment_effect_SBP_input,
                                              cost_disc_rate_input = discount_cost_input, 
                                              qol_disc_rate_input = discount_util_input,
                                              retirement_age_input = retirement_age_input, 
                                              inf_care_hours_input = inf_care_hours_input, 
                                              worked_hours_input = worked_hours_input,
                                              working_days_lost_input = working_days_lost_input,
                                              cost_hour_sick_input = cost_hour_sick_input,
                                              friction_period_input = friction_period_input, 
                                              informal_care_coef_input = informal_care_coef_input,
                                              employment_coef_input = employment_coef_input,
                                              prod_costs_coef_input = prod_costs_coef_input, 
                                              inf_care_age_scale_input = inf_care_age_scale_input, 
                                              prod_loss_age_scale_input = prod_loss_age_scale_input,
                                              psa_input = psa_input,
                                              seed_input = seed_input)


sim.results.female.comp <- SMDMII_model_simulation(patient_size_input = npats_input,
                                                   female_input = 1,
                                                   tx_cost_input = tx_cost_input,
                                                   treatment_effect_HbA1c_input = rep(0,4),
                                                   treatment_effect_HDL_input = rep(0,4),
                                                   treatment_effect_LDL_input = rep(0,4),
                                                   treatment_effect_BMI_input = rep(0,4),
                                                   treatment_effect_SBP_input = rep(0,4),
                                                   cost_disc_rate_input = discount_cost_input, 
                                                   qol_disc_rate_input = discount_util_input,
                                                   retirement_age_input = retirement_age_input, 
                                                   inf_care_hours_input = inf_care_hours_input, 
                                                   worked_hours_input = worked_hours_input,
                                                   working_days_lost_input = working_days_lost_input,
                                                   cost_hour_sick_input = cost_hour_sick_input,
                                                   friction_period_input = friction_period_input, 
                                                   informal_care_coef_input = informal_care_coef_input,
                                                   employment_coef_input = employment_coef_input,
                                                   prod_costs_coef_input = prod_costs_coef_input, 
                                                   inf_care_age_scale_input = inf_care_age_scale_input, 
                                                   prod_loss_age_scale_input = prod_loss_age_scale_input,
                                                   psa_input = psa_input,
                                                   seed_input = seed_input)

sim.results.male <- SMDMII_model_simulation(patient_size_input = npats_input,
                                            female_input = 0,
                                            tx_cost_input = tx_cost_input,
                                            treatment_effect_HbA1c_input = treatment_effect_HbA1c_input, 
                                            treatment_effect_HDL_input = reatment_effect_HDL_input,
                                            treatment_effect_LDL_input = treatment_effect_LDL_input, 
                                            treatment_effect_BMI_input = treatment_effect_BMI_input,
                                            treatment_effect_SBP_input = treatment_effect_SBP_input,
                                            cost_disc_rate_input = discount_cost_input, 
                                            qol_disc_rate_input = discount_util_input,
                                            retirement_age_input = retirement_age_input, 
                                            inf_care_hours_input = inf_care_hours_input, 
                                            worked_hours_input = worked_hours_input,
                                            working_days_lost_input = working_days_lost_input,
                                            cost_hour_sick_input = cost_hour_sick_input,
                                            friction_period_input = friction_period_input, 
                                            informal_care_coef_input = informal_care_coef_input,
                                            employment_coef_input = employment_coef_input,
                                            prod_costs_coef_input = prod_costs_coef_input, 
                                            inf_care_age_scale_input = inf_care_age_scale_input, 
                                            prod_loss_age_scale_input = prod_loss_age_scale_input,
                                            psa_input = psa_input,
                                            seed_input = seed_input)

sim.results.male.comp <- SMDMII_model_simulation(patient_size_input = npats_input,
                                                 female_input = 0,
                                                 tx_cost_input = tx_cost_input,
                                                 treatment_effect_HbA1c_input = rep(0,4),
                                                 treatment_effect_HDL_input = rep(0,4),
                                                 treatment_effect_LDL_input = rep(0,4),
                                                 treatment_effect_BMI_input = rep(0,4),
                                                 treatment_effect_SBP_input = rep(0,4),
                                                 cost_disc_rate_input = discount_cost_input, 
                                                 qol_disc_rate_input = discount_util_input,
                                                 retirement_age_input = retirement_age_input, 
                                                 inf_care_hours_input = inf_care_hours_input, 
                                                 worked_hours_input = worked_hours_input,
                                                 working_days_lost_input = working_days_lost_input,
                                                 cost_hour_sick_input = cost_hour_sick_input,
                                                 friction_period_input = friction_period_input, 
                                                 informal_care_coef_input = informal_care_coef_input,
                                                 employment_coef_input = employment_coef_input,
                                                 prod_costs_coef_input = prod_costs_coef_input, 
                                                 inf_care_age_scale_input = inf_care_age_scale_input, 
                                                 prod_loss_age_scale_input = prod_loss_age_scale_input,
                                                 psa_input = psa_input,
                                                 seed_input = seed_input)


# Print simulation duration
end <- Sys.time()
print(end - init)

# Save simulation results
save(sim.vars,
     sim.results.female,
     sim.results.female.comp,
     sim.results.male,
     sim.results.male.comp,
     file = 'output/Rank6_spec_targetpop.RData')


