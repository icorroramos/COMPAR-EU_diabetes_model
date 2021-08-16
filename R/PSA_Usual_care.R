##################
# Set-up options #
##################

# Clear global environment, flush memory, load packages
source('R/PSA_start_chunk.R')

# Variable defined to keep track of simulation time (delete afterwards)
init <- Sys.time()

#
# Create empty objects to store PSA data in
sim.var.data <- matrix(nrow = n_psa_input, ncol = 11)
female.comp.data <- list()
male.comp.data <- list()

for (i in 1:n_psa_input){

        # RESULT OUTPUT FOR OUTPUT SCRIPT -----------------------------------------
        sim.results.female.comp <- SMDMII_model_simulation(npats_input,
                                                           female_input = 1,
                                                           tx_cost_input,
                                                           treatment_effect_HbA1c_input = rep(0,4), 
                                                           treatment_effect_HDL_input = rep(0,4),
                                                           treatment_effect_LDL_input = rep(0,4), 
                                                           treatment_effect_BMI_input = rep(0,4),
                                                           treatment_effect_SBP_input = rep(0,4),
                                                           discount_cost_input, 
                                                           discount_util_input, 
                                                           retirement_age_input, 
                                                           psa_input,
                                                           seed_input)

        sim.results.male.comp <- SMDMII_model_simulation(npats_input,
                                                         female_input = 0,
                                                         tx_cost_input,
                                                         treatment_effect_HbA1c_input = rep(0,4),
                                                         treatment_effect_HDL_input = rep(0,4),
                                                         treatment_effect_LDL_input = rep(0,4),
                                                         treatment_effect_BMI_input = rep(0,4),
                                                         treatment_effect_SBP_input = rep(0,4),
                                                         discount_cost_input,
                                                         discount_util_input,
                                                         retirement_age_input,
                                                         psa_input,
                                                         seed_input)
        
        # Store iteration results in data objects

        
       female.comp.data[[i]] <- unlist(sim.results.female.comp[-1])
       male.comp.data[[i]] <- unlist(sim.results.male.comp[-1])
       
}

psa.results.female.comp <- as.data.frame(do.call(rbind, female.comp.data))
psa.results.male.comp <- as.data.frame(do.call(rbind, male.comp.data))


# Print simulation duration
end <- Sys.time()
print(end - init)


# Save simulation results
save(psa.results.female.comp,
     psa.results.male.comp,
     file = 'output/PSA_Usual_care.RData')
