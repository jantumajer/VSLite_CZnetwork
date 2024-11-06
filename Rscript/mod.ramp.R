
mod.ramp <- function(x,x1,x2,x3,x4){
  
  growth.rate <- matrix(nrow = 12, ncol = ncol(x))
  
  for (i in c(1:ncol(x))){
    for (j in c(1:12)){
      
      if(x[j,i] < x1){growth.rate[j,i] <- 0} # No growth due to LOW limitation
        else {if(x[j,i] > x1 & x[j,i] < x2){growth.rate[j,i] <- (x[j,i] - x1)/(x2-x1)} # Suboptimal growth due to LOW limitation
          else {if(x[j,i] > x2 & x[j,i] < x3){growth.rate[j,i] <- 1} # Optimal growth
            else {if(x[j,i] > x3 & x[j,i] < x4){growth.rate[j,i] <- (x4 - x[j,i])/(x4-x3)} # Suboptimal growth due to HIGH limitation
              else {if(x[j,i] > x4){growth.rate[j,i] <- 0}}}}} # No growth due to HIGH limitation
    }
  }
  
return(growth.rate)
}


