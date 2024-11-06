# Compute integral growth rate from the partial growth rates of temperature and soil moisture
# Multiplicative version according to Tumajer et al. 2017: Agricultural and Forest Meteorology, 247: 56-64.

integrate.multiplic <- function(gt, gm, ge){
  
  gint <- matrix(nrow = 12, ncol = ncol(gt))
  
  for (i in c(1:ncol(gint))){
    for(j in c(1:12)){
      
    gint[j,i] <- ge[j,1] * gt[j,i] * gm[j,i]
    
    }
  }
  
  return(gint)
}