#' VS-Lite model of tree ring width growth.
#' 
#' \code{VSLite} simulates tree ring width growth.
#' 
#' R port of VS-Lite Model of Tree Ring Width by Suz TOlwinski-Ward, 2015. For more references,
#' see xxxxyyyyyzzzz.
#' 
#' @param syear Start year of simulation.
#' @param eyear End year of simulation.
#' @param phi Latitude of site (in degrees N).
#' @param T (12 x Nyrs) Matrix of ordered mean monthly temperatures (in degEes C).
#' @param P (12 x Nyrs) Matrix of ordered accumulated monthly precipitation (in mm).
#' @param T1 Lower temperature threshold for growth to begin (scalar, deg. C).
#' @param T2 Upper temperature threshold for growth sensitivity to temp (scalar, deg. C).
#' @param M1 Lower moisture threshold for growth to begin (scalar, v.v).
#' @param M2 Upper moisture threshold for growth sensitivity to moisture (scalar, v/v).
#' @param Mmax Scalar maximum soil moisture held by the soil (in v/v).
#' @param Mmin Scalar minimum soil moisture (for error-catching) (in v/v).
#' @param alph Scalar runoff parameter 1 (in inverse months).
#' @param m.th Scalar runoff parameter 3 (unitless).
#' @param mu.th Scalar runoff parameter 2 (unitless).
#' @param rootd Scalar root/"bucket" depth (in mm).
#' @param M0 Initial value for previous month's soil moisture at t = 1 (in v/v).
#' @param substep Use leaky bucket code with sub-monthly time-stepping? (TRUE/FALSE) 
#' @param I_0 lower bound of integration window (months before January in NH)
#' @param I_f upper bound of integration window (months after January in NH)
#' @param hydroclim Switch; value is either "P" (default) or "M" depending on whether the 
#' second input climate variable is precipitation, in which case soil moisture is estimated
#' using the Leaky Bucket model of the CPC, or soil moisture, in which case the inputs are 
#' used directly to compute the growth response.
#' 
#' @return trw
#' @return gT
#' @return gM
#' @return gE
#' @return M
#' @return potEv
#' @return sample.mean.width 
#' @return sample.std.width 
#' 
#' @seealso \code{\link{compute.gE}},\code{\link{std.ramp}},\code{\link{leakybucket.monthly}},\code{\link{leakybucket.submonthly}}
#'
#' @export
####################################################################################################


VSLite <- function(syear, eyear, phi, Tinput, Pinput,
                        T1 = 8, T2 = 23, T3 = 28, T4= 35,
                        M1 = 0.01, M2 = 0.05, M3 = 0.7, M4 = 0.85,
                        Mmax = 0.76, Mmin = 0.01, alph = 0.093,
                        m.th = 4.886, mu.th = 5.8, rootd = 1000, M0 = 0.2,
                        substep = 0, hydroclim = "P",
                        I_0 = 1, I_f = 12, Acor = 1,
                        ramp = "orig",
                        integration = "orig",
                        corr = NULL # Pouze pro ulozeni hodnoty korelacniho koeficientu nejlepsiho modelu
                        ){
  #############################################################################
  nyrs <- length(syear:eyear)
  Gr <- gT <- gM <- M <- matrix(NA,12,nyrs);
  
  # Subseting climatic data
  T <- t(Tinput[Tinput[,1] %in% c(syear : eyear), c(2:13)])
  P <- t(Pinput[Pinput[,1] %in% c(syear : eyear), c(2:13)])
  #############################################################################
  
  ## Load in soil moisture, or estimate it with the Leaky Bucket model:
  if(hydroclim == "M"){
    ## Read in soil moisture:
    M = P;
  }else{# Compute soil moisture:
    if(substep == 1){
      M <- leakybucket.submonthly(syear,eyear,phi,T,P,
                                  Mmax,Mmin,alph,m.th,mu.th,rootd,M0);
      potEV <- NULL
      BAL <- NULL
    }else{
      M <- leakybucket.monthly(syear,eyear,phi,T,P,
                               Mmax,Mmin,alph,m.th,mu.th,rootd,M0)$M;
      potEV <- leakybucket.monthly(syear,eyear,phi,T,P,
                               Mmax,Mmin,alph,m.th,mu.th,rootd,M0)$potEV;
      BAL <- P - potEV
    }
    if(substep !=1 && substep != 0){
      cat("'substep' param must either be set to 1 or 0.");
      return
    }
  }
  
  # Compute gE, the scaled monthly proxy for insolation:
  gE <- compute.gE(phi);
  
  #############################################################################
  
  ### Calculate Growth Response functions gT and gM using one of available ramp functions
  if(ramp == "orig") {  
    gT <- std.ramp(T, T1, T2)
    gM <- std.ramp(M, M1, M2)}
  
  if(ramp == "modif") {  
    gT <- mod.ramp(T, T1, T2, T3, T4)
    gM <- mod.ramp(M, M1, M2, M3, M4)}
  
  ### Compute overall growth rate using one of available integration functions
  if(integration == "orig") { 
    Gr <- integrate.orig(gT, gM, gE)}
  
  if(integration == "modif") { 
    Gr <- integrate.multiplic(gT, gM, gE)}
 
  ############## Compute proxy quantity from growth responses #################
  width <- matrix(NA,nyrs,1);
  if (phi>0){ # Site in Northern Hemisphere:
    if (I_0<0){ # if we include part of the previous year in each year's modeled growth:
      startmo <- 13+I_0;
      endmo <- I_f;
      # use average of growth data across modeled years to estimate first year's growth due
      # to previous year:
      width[1] <- sum(Gr[1:endmo,1]) + Acor * sum(rowMeans(Gr[startmo:12,]));
      for(cyear in 2:nyrs){
        width[cyear] <- sum(Gr[1:endmo,cyear]) + Acor * sum(Gr[startmo:12,cyear-1]);
      }
    }else{ # no inclusion of last year's growth conditions in estimates of this year's growth:
      startmo <- I_0+1;
      endmo <- I_f;
      width <- colSums(Gr[startmo:endmo,])
    }
  }
  if(phi<0){ # if site is in the Southern Hemisphere:
    # (Note: in the Southern Hemisphere, ring widths are dated to the year in which growth began!)
    startmo <- 7+I_0; # (eg. I_0 = -4 in SH corresponds to starting integration in March of cyear)
    endmo <- I_f-6; # (eg. I_f = 12 in SH corresponds to ending integraion in June of next year)
    for (cyear in 1:(nyrs-1)){
      width(cyear) <- sum(Gr[startmo:12,cyear]) + sum(Gr[1:endmo,cyear+1]);
    }
    # use average of growth data across modeled years to estimate last year's growth due
    # to the next year:
    width[nyrs] <- sum(Gr[startmo:12,nyrs])+sum(rowMeans(Gr[1:endmo,]));
  }
  
  # Simulated proxy series standardized width:
  trw <- t((width-mean(width))/sd(width)); 

  #############################################################################
  # Return output:
  out <- list(mod.trw = trw, gT = gT, gM = gM, gINT = Gr, gE = gE,
              Moist = M, Temp = T, potEV = potEV, BAL = BAL,
              syear = syear, eyear = eyear,
              par = list(T1 = T1, T2 = T2, T3 = T3, T4 = T4,
                         M1 = M1, M2 = M2, M3 = M3, M4 = M4,
                         Acor = Acor, I_0 = I_0,
                         corr = corr))
  return(out)

}


