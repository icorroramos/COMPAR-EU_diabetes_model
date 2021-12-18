### Function to combine the results of two model runs.
### Maintains names and structure of output objects so that analysis scripts can be used

stitcher <- function(comp.name, location.1, location.2){
  load(paste0(location.2, comp.name, '.RData'))
  npats.loc2 <- ifelse(exists('sim.vars'), sim.vars[[2]], sim.vars.uc[[2]])
  
  sim.names <- ls()[startsWith(ls(), 'sim.results')]
  
  for (i in 1:length(sim.names)) {
    assign(paste0(sim.names[i], '.extend'), get(sim.names[i]))
  }
  
  load(paste0(location.1, comp.name, '.RData'))
  
  npats.loc1 <- ifelse(exists('sim.vars'), sim.vars[[2]], sim.vars.uc[[2]])
  
  for (i in 1:length(sim.names)) {
    loc2.data <- get(paste0(sim.names[i], '.extend'))
    loc2.data[[1]]$SIMID <- loc2.data[[1]]$SIMID + npats.loc1
    loc1.data <- get(sim.names[i])
    loc1.data[[1]] <- rbind(loc1.data[[1]], loc2.data[[1]])
    loc1.data[-1] <- (npats.loc1 * as.data.frame(loc1.data[-1]) + npats.loc2 * as.data.frame(loc2.data[-1])) / (npats.loc1 + npats.loc2)
    assign(sim.names[i],  loc1.data, envir = .GlobalEnv)
  }
}


### Function to combine n model runs stored in same directory
### Maintains names and structure of output objects so that analysis scripts can be used

sewing_machine <- function(n.run, filename) {
  n.pats.run <- rep(NA, n.run)
  run.seed <-rep(NA, n.run)
  
  load(paste0(filename, '1.RData')) # Get first output to check which simulations were run
  sim.names <- ls()[startsWith(ls(), 'sim.results')]
  list.names <- names(get(sim.names[1]))
  n.sims <- length(sim.names)
  
  all.means <- lapply(1:n.sims, matrix, data = NA, nrow = 27, ncol = n.run)
  all.patients <- list() 
  
  for (i in 1:n.run) {
    load(paste0(filename, i, '.RData'))
    sim.names.run <- ls()[startsWith(ls(), 'sim.results')]
    if (!identical(sim.names,sim.names.run[1:n.sims])) stop(paste('Simulations in run', i, 'not the same as in first run'))
    n.pats.run[i] <- ifelse(exists('sim.vars'), sim.vars[[2]], sim.vars.uc[[2]])
    run.seed[i] <- ifelse(exists('sim.vars'), sim.vars[[1]], sim.vars.uc[[1]])
    for (j in 1:n.sims) {
      curr.sim <- get(sim.names[j])
      all.means[[j]][, i] <- unlist(curr.sim[2:28])
      all.patients[[j]] <- rbind.fill(curr.sim[1])
    }
  }
  
  run.weight <- n.pats.run / sum(n.pats.run)
  weighed.means <- lapply(all.means, '*', run.weight)
  unnamed.means <- lapply(weighed.means, rowSums)
  results.means <- lapply(unnamed.means, function(x) {names(x) <- list.names[2:28]; return(x)})
  
  for (i in 1:n.sims){
    results.patient <- list(all.patients[[i]])
    names(results.patient) <- list.names[1]
    assign(sim.names[i], c(results.patient, as.list(results.means[[i]])), envir = .GlobalEnv)
    print(paste('Created object', sim.names[i], 'in global environment'))
  }
  assign('sewing.seeds', run.seed, envir = .GlobalEnv)
  print('Created object sewing.seeds containing all seed values')
}
