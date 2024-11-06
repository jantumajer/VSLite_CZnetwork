# Implementing the approach from:  
# Ols, C., Klesse, S., Girardin, M.P., Evans, M.E.K., DeRose, R.J., Trouet, V., 2023. Detrending climate data prior to climate–growth analyses in dendroecology: A common best practice? Dendrochronologia 79, 126094. https://doi.org/10.1016/j.dendro.2023.126094

library(dplR)

climate.detrend <- function(climate, var, spline){
  
  if (var == "temp"){
    climate.d <- (dplR::detrend(climate+100, method = "Spline", nyrs = spline, return.info = T, difference = T)$curves) - 100 # Trend fits
    climate.d2 <- climate - climate.d # Subtraction (!)
    
    for (l in c(2:13)){
      climate.d2[,l] <- climate.d2[,l] + colMeans(climate)[l] # Scaling
      }
  }
  
  if (var == "prec"){
    climate.d <- (dplR::detrend(climate+100, method = "Spline", nyrs = spline, return.info = T, difference = T)$curves) - 100 # Trend fits
    climate.d2 <- climate / climate.d # Division (!)
    
    for (l in c(2:13)){
      climate.d2[,l] <- climate.d2[,l] * colMeans(climate)[l] # Scaling
    }
  }
  
  climate.d2[,1] <- climate[,1]
  
  return(climate.d2)
}