
stitcher <- function(comp.name, location.1, location.2){
  load(paste0(location.2, comp.name, '.RData'))
  npats.loc2 <- sim.vars[[2]]
  
  sim.names <- ls()[startsWith(ls(), 'sim.results')]
  
  for (i in 1:length(sim.names)) {
    assign(paste0(sim.names[i], '.extend'), get(sim.names[i]))
  }
  
  load(paste0(location.1, comp.name, '.RData'))
  
  npats.loc1 <- sim.vars[[2]]
  
  for (i in 1:length(sim.names)) {
    loc2.data <- get(paste0(sim.names[i], '.extend'))
    loc2.data[[1]]$SIMID <- loc2.data[[1]]$SIMID + npats.loc1
    loc1.data <- get(sim.names[i])
    loc1.data[[1]] <- rbind(loc1.data[[1]], loc2.data[[1]])
    loc1.data[-1] <- (npats.loc1 * as.data.frame(loc1.data[-1]) + npats.loc2 * as.data.frame(loc2.data[-1])) / (npats.loc1 + npats.loc2)
    assign(sim.names[i],  loc1.data, envir = .GlobalEnv)
  }
}
