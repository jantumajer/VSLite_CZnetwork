
randomization <- function(iter = 5000){
  
  random.par <- data.frame(iter = c(1:iter), T1 = NA, T2 = NA, T3 = NA, T4 = NA,
                                             M1 = NA, M2 = NA, M3 = NA, M4 = NA,
                                             I_0 = NA, Acor = NA)
  
  # Temperature response function
  random.par[,"T1"] <- runif(iter, min = 2, max = 9) # min = 0
  random.par[,"T2"] <- runif(iter, min = 9, max = 20)
  random.par[,"T3"] <- runif(iter, min = (random.par[,"T2"] + 0.5), max = 30)
  random.par[,"T4"] <- runif(iter, min = pmax(random.par[,"T3"], 20), max = 35)
  
  # Moisture response function
  random.par[,"M1"] <- runif(iter, min = 0.0, max = 0.3) # max = 0.4
  random.par[,"M2"] <- runif(iter, min = (random.par[,"M1"] + 0.025), max = 0.5) # max = 0.6
  random.par[,"M3"] <- 0.80 # runif(iter, min = pmax(random.par[,"M2"], 0.5), max = 0.8)
  random.par[,"M4"] <- 0.85 # Alternatively, play with runif(...)
  
  # Autocorrelation in tree-ring width series
  random.par[,"I_0"] <- sample(-12:-2, iter, replace = T)
  random.par[,"Acor"] <- runif(iter, min = -1.5, max = 1.5)
  
  return(random.par)
  
}
