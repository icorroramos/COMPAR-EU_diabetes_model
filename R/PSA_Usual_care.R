##################
# Set-up options #
##################

# Clear global environment, flush memory, load packages
source('R/PSA_start_chunk.R')

# Variable defined to keep track of simulation time (delete afterwards)
init <- Sys.time()

#
# Create empty objects to store PSA data in
sim.var.data <- matrix(nrow = n_psa_input, ncol = 2)
female.comp.data <- list()
male.comp.data <- list()

for (i in 1:n_psa_input){

        # RESULT OUTPUT FOR OUTPUT SCRIPT -----------------------------------------
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
                                                           seed_input = psa.seed[i])

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
                                                         seed_input = psa.seed[i])
        
        # Store iteration results in data objects

       sim.var.data[i, ] <- c(psa.seed[i], npats_input) 
       female.comp.data[[i]] <- unlist(sim.results.female.comp[-1])
       male.comp.data[[i]] <- unlist(sim.results.male.comp[-1])
       
}

colnames(sim.var.data) <- c('Seed', 'Npats')
psa.results.female.comp <- as.data.frame(do.call(rbind, female.comp.data))
psa.results.male.comp <- as.data.frame(do.call(rbind, male.comp.data))


# Print simulation duration
end <- Sys.time()
print(end - init)


# Save simulation results
save(sim.var.data,
     psa.results.female.comp,
     psa.results.male.comp,
     file = 'output/PSA_Usual_care.RData')
