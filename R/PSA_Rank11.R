##################
# Set-up options #
##################

# Clear global environment, flush memory, load packages
source('R/PSA_start_chunk.R')

# Variable defined to keep track of simulation time (delete afterwards)
init <- Sys.time()


# Get treatment effect draws
draws.hba1c <- rpert(n_psa_input, min = -1.0843, mode = -0.5986, max = -0.113)
draws.hdl <- rpert(n_psa_input, min = 0.7186, mode = 3.579, max = 6.4407)
draws.ldl <- rpert(n_psa_input, min = -21.7249, mode = -16.632, max = -11.539)
draws.bmi <- rpert(n_psa_input, min = -1.5356, mode = -0.4006, max = 0.7344)
draws.sbp <- rpert(n_psa_input, min = -9.6507, mode = -1.1086, max = 7.4336)


# Create empty objects to store PSA data in
sim.var.data <- matrix(nrow = n_psa_input, ncol = 11)
female.int.data <- list()
female.comp.data <- list()
male.int.data <- list()
male.comp.data <- list()

for (i in 1:n_psa_input){
        # Treatment effect inputs
        treateff_start   <- 1 # Cycle in which treatment effect starts
        treateff_end     <- 3 # Cycle in which treatment effect ends
        treateff_decline <- 2 # Cycle in which treatment effect starts to decline (linearly)
        
        treateff_hba1c <- draws.hba1c[i] # Treatment effect on HbA1c (in absolute %-points HbA1c)
        treateff_hdl   <- draws.hdl[i] * 0.02586 * 10 # TRANSFORMATION *10 FOR MODEL INPUT # Treatment effect on HDL-cholesterol (absolute effect in mmol/l)
        treateff_ldl   <- draws.ldl[i] * 0.02586 * 10 # TRANSFORMATION *10 FOR MODEL INPUT# Treatment effect on LDL-cholesterol (absolute effect in mmol/l)
        treateff_bmi <- draws.bmi[i] # Treatment effect on BMI (in absolute points)
        treateff_sbp <- draws.sbp[i] / 10 # TRANSFORMATION /10 FOR MODEL INPUT# Treatment effect on SBP (in absolute mmHg)
        

        # Tx effects are vectors: the current assumption is that the same start, end and decline is assumed for all effect modifiers
        treatment_effect_HbA1c_input <- c(treateff_hba1c, treateff_start, treateff_end, treateff_decline)
        treatment_effect_HDL_input <- c(treateff_hdl, treateff_start, treateff_end, treateff_decline)
        treatment_effect_LDL_input <- c(treateff_ldl, treateff_start, treateff_end, treateff_decline)
        treatment_effect_BMI_input <- c(treateff_bmi, treateff_start, treateff_end, treateff_decline)
        treatment_effect_SBP_input <- c(treateff_sbp, treateff_start, treateff_end, treateff_decline)

        
        # RESULT OUTPUT FOR OUTPUT SCRIPT -----------------------------------------
        sim.results.female <- SMDMII_model_simulation(patient_size_input = npats_input,
                                                      female_input = 1,
                                                      tx_cost_input = tx_cost_input,
                                                      treatment_effect_HbA1c_input = treatment_effect_HbA1c_input, 
                                                      treatment_effect_HDL_input = treatment_effect_HDL_input,
                                                      treatment_effect_LDL_input = treatment_effect_LDL_input, 
                                                      treatment_effect_BMI_input = treatment_effect_BMI_input,
                                                      treatment_effect_SBP_input = treatment_effect_SBP_input,
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
                                                      seed_input = psa.seed[i])
        
        sim.results.male <- SMDMII_model_simulation(patient_size_input = npats_input,
                                                    female_input = 0,
                                                    tx_cost_input = tx_cost_input,
                                                    treatment_effect_HbA1c_input = treatment_effect_HbA1c_input, 
                                                    treatment_effect_HDL_input = treatment_effect_HDL_input,
                                                    treatment_effect_LDL_input = treatment_effect_LDL_input, 
                                                    treatment_effect_BMI_input = treatment_effect_BMI_input,
                                                    treatment_effect_SBP_input = treatment_effect_SBP_input,
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
                                                    seed_input = psa.seed[i])
        
        # Store iteration results in data objects
        sim.var.data[i, ] <- c(psa.seed[i], npats_input, tx_cost_input, unlist(mget(apropos('treateff.'))))
        female.int.data[[i]] <- unlist(sim.results.female[-1])
        male.int.data[[i]] <- unlist(sim.results.male[-1])
}

colnames(sim.var.data) <- c('Seed', 'Npats', 'Tx_Cost', 'BMI', 'Decline', 'End', 'HbA1c', 'HDL', 'LDL', 'SBP', 'Start')
psa.results.female <- as.data.frame(do.call(rbind, female.int.data))
psa.results.male <- as.data.frame(do.call(rbind, male.int.data))

# Print simulation duration
end <- Sys.time()
print(end - init)


# Save simulation results
save(sim.var.data,
     psa.results.female,
     psa.results.male,
     file = 'output/PSA_Rank11.RData')
