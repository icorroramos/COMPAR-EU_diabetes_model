##################
# Set-up options #
##################

# Clear global environment, flush memory, load packages
source('R/Analysis_start_chunk.R')

# Variable defined to keep track of simulation time (delete afterwards)
init <- Sys.time()

# The sim.vars object collects all parameters that define the simulation into one object
# Then, it is saved with the output of the simulation. This way, if we  have multiple output files, we always have the information 
# on the relevant input parameters that were used to produce the output.
sim.vars.uc <- list(seed_input, npats_input, tx_cost_input, mget(apropos('treateff.')))



# RESULT OUTPUT FOR OUTPUT SCRIPT -----------------------------------------
sim.results.female.comp <- SMDMII_model_simulation(patient_size_input = npats_input,
                                                   female_input = 1,
                                                   tx_cost_input = tx_cost_input,
                                                   treatment_effect_HbA1c_input = rep(0,4),
                                                   treatment_effect_HDL_input = rep(0,4),
                                                   treatment_effect_LDL_input = rep(0,4),
                                                   treatment_effect_BMI_input = rep(0,4),
                                                   treatment_effect_SBP_input = rep(0,4),
                                                   cost_disc_rate_input = discount_factors[1], 
                                                   qol_disc_rate_input = discount_factors[2],
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
                                                   run_PSA_input = psa_input,
                                                   seed_input = seed_input)


sim.results.male.comp <- SMDMII_model_simulation(patient_size_input = npats_input,
                                                 female_input = 0,
                                                 tx_cost_input = tx_cost_input,
                                                 treatment_effect_HbA1c_input = rep(0,4),
                                                 treatment_effect_HDL_input = rep(0,4),
                                                 treatment_effect_LDL_input = rep(0,4),
                                                 treatment_effect_BMI_input = rep(0,4),
                                                 treatment_effect_SBP_input = rep(0,4),
                                                 cost_disc_rate_input = discount_factors[1], 
                                                 qol_disc_rate_input = discount_factors[2],
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
                                                 run_PSA_input = psa_input,
                                                 seed_input = seed_input)


# Print simulation duration
end <- Sys.time()
print(paste(end, 'End of Usual care analysis'))
print('Duration of Usual care analysis')
print(end-init)

# Save simulation results
save(sim.vars.uc,
     sim.results.female.comp,
     sim.results.male.comp,
     file = paste0('output/', country.id, '/Usual_care_outcomes_', npats_input, '-pats_seed-', seed_input, '.RData'))

