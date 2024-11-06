

VSLite.iterative <- function(parameters, obs.chron, chron.variant = "std",
                           syear, eyear, phi, Tinput, Pinput,
                           ramp = "orig",
                           integration = "orig") {
  
  # Subsetting observed chronology
  obs.chron <- obs.chron[as.numeric(row.names(obs.chron)) %in% c(syear : eyear), ]
  
  # Selection of chronology variant
  obs.chron <- obs.chron[,chron.variant]
  
  # Z-transformation of observed chronology
  obs.chron.z <- (obs.chron - mean(obs.chron))/sd(obs.chron)
  
  for (iter in c(1:nrow(parameters))){
  
  sim <- VSLite(syear, eyear, phi, Tinput, Pinput,
                T1 = parameters[iter , "T1"], T2 = parameters[iter , "T2"], T3 = parameters[iter , "T3"], T4 = parameters[iter , "T4"],
                M1 = parameters[iter , "M1"], M2 = parameters[iter , "M2"], M3 = parameters[iter , "M3"], M4 = parameters[iter , "M4"],
                Acor = parameters[iter, "Acor"], I_0 =  parameters[iter, "I_0"],
                ramp = ramp, integration = integration)
  
  parameters[iter, "correlation"] <- cor(t(sim$mod.trw), obs.chron.z) # Correlation between observed and simulated chronology
  if(min(colSums(sim$gINT)) == 0){parameters[iter, "correlation"] <- -1} #NEW
  
  rm(sim)
    
  }
  
  return(list(random.table = parameters, 
              best = parameters[parameters$correlation == max(parameters$correlation, na.rm = T) & !is.na(parameters$correlation), ],
              ramp = ramp,
              integration = integration,
              obs.trw = obs.chron))
}
