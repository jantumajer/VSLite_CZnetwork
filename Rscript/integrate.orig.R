# Compute integral growth rate from the partial growth rates of temperature and soil moisture
# Original version according to Tolwinski-Ward et al. 2011

integrate.orig <- function(gt, gm, ge){
  
  gint <- matrix(nrow = 12, ncol = ncol(gt))
  
  for (i in c(1:ncol(gint))){
    gint[,i] <- ge * pmin(gt[,i], gm[,i])
  }
  
  return(gint)
}