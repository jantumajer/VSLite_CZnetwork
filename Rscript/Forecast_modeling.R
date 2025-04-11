library(RColorBrewer); library(plot.matrix); library(ggplot2); library(reshape2); library(readxl)
library(e1071); library(multimode); library(pracma)

library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis

#################

# Model definitions
setwd("E:/VSLite_R/R/")

## 0a] Model sub-algorithms ##
# Calculation of partial growth rates to photoperiod
source("compute.gE.R") 
source("daylength.factor.from.lat.R") 

# Soil moisture model
source("leakybucket.monthly_annualTm.R") 
source("leakybucket.submonthly.R") 

# Ramp functions
source("std.ramp.R") # Original non-declining ramp functions
source("mod.ramp.R") # Modified increasing-stable-decreasing ramp functions

# Integration functions
source("integrate.orig.R") # Integration based on MINIMUM of growth rates (following Liebig's law)
source("integrate.multiplic.R") # Integration based on PRODUCT of growth rates (following initial TRACH model)

# Model definitions
source("VSLite.R")

## 0b] Functions to calibrate the model against local site chronology ##
source("randomization_MonikaDuby.R")
source("VSLite.iterative.R")

## 0c] Graphical functions
source("charts.R")

## 0d] detrending climatic data
source("climate.detrend.R")

#################

SITE <- read.csv("E:/TACR/zchron_TestData2/META/data_april2_withCZU.csv", row.names = 1)
ADD <- read_excel("E:/TACR/zchron_TestData2/ADD/additions_202303.xlsx", sheet = "site")
ADD2 <- read_excel("E:/TACR/zchron_TestData2/ADD2/additions_202401.xlsx", sheet = "site")
ADD_Senfeldr <- read.csv("E:/TACR/zchron_TestData2/ADD_Senfeldr/add_senfeldr.csv")
ADD_Mendelu <- read_excel("E:/TACR/zchron_TestData2/ADD_Mendelu/Mendelu_site_data_nove2023.xlsx")
ADD_IFER <- read_excel("E:/TACR/VS_for_GACR/IFER/metadata_IFER.xlsx")
SITE <- rbind(SITE, ADD[,c(2:50)], ADD2[,c(2:50)], ADD_Senfeldr, ADD_Mendelu, ADD_IFER)
rm(ADD, ADD2, ADD_Senfeldr, ADD_Mendelu, ADD_IFER)
SITE <- SITE[!(SITE$raw_data_file_name == ""),] # Series from Cada and Rydval still missing

SITEb <- SITE[SITE$last_year > 1994 | is.na(SITE$last_year),]
SITE <- SITE[SITE$site_code %in% SITEb$site_code,]

SITE2 <- read_xlsx("E:/TACR/VS_for_GACR/GACR_Model_Revision/metadata_GACR_SK_REMOTEF.xlsx")
SITE <- rbind(SITE, SITE2)

nyr <- 60

#################

for (i in c(1:nrow(SITE))){
  site <- SITE[i, "site_code"] # Which site is processed
  genus <- substr(SITE[i, "species"], 1, 2)
  phi <- SITE[i, "site_lat_decimal"]
  lon <- SITE[i, "site_long_decimal"]
  
  
  if(paste(site, ".Rda", sep ="") %in% list.files(paste("e:/TACR/VS_for_GACR/Results/", genus, sep = ""))){
  
  # Model parameters
  model <- readRDS(paste("e:/TACR/VS_for_GACR/Results/", genus,"/", site, ".Rda", sep = ""))
  par <- model$par

  ### Loading and preprocessing climatic data
  clim <- read.csv(paste("E:/TACR/Climate_Grids/Climate_tables/", site, "_clim.csv", sep = ""))
  clim_sub <- clim[clim$Year %in% c(1995:2014),]
  
  # Detrending climatic data
  temperature <- dcast(clim[clim$Year < 2024, ], formula = Year ~ Month, value.var = "Temp")
  precipitation <- dcast(clim[clim$Year < 2024, ], formula = Year ~ Month, value.var = "Prec")
  temperature.det <- climate.detrend(temperature, var = "temp", spline = nyr)
  precipitation.det <- climate.detrend(precipitation, var = "prec", spline = nyr)
  temperature.det_sub <- temperature.det[temperature.det$Year %in% c(1995:2014),]
  precipitation.det_sub <- precipitation.det[precipitation.det$Year %in% c(1995:2014),]
  
  clim_agg <- cbind(aggregate(clim_sub[,c("Temp", "Prec")], by = list(Month = clim_sub$Month), FUN = mean),
                    Temp_15perc =   sapply(temperature.det_sub[,c(2:13)], quantile, probs = c(0.15)),
                    Temp_85perc =   sapply(temperature.det_sub[,c(2:13)], quantile, probs = c(0.85)),
                    Prec_25perc =   sapply(precipitation.det_sub[,c(2:13)], quantile, probs = c(0.25)),
                    Prec_75perc =   sapply(precipitation.det_sub[,c(2:13)], quantile, probs = c(0.75)))
  
  # clim_agg <- data.frame(Month = c(1:12), Temp = colMeans(temperature.det[temperature.det$Year %in% c(1995:2014),c(2:13)]),
  #                                         Prec = colMeans(precipitation.det[precipitation.det$Year %in% c(1995:2014),c(2:13)]))
  
  # Load SSP anomalies
  anomaliesT_85 <- read.csv(paste("E:/TACR/WBG_projekce/Ensemble_nc/Climate_tables/", site, "_anomaliesT_85.csv", sep =""))
  anomaliesT_70 <- read.csv(paste("E:/TACR/WBG_projekce/Ensemble_nc/Climate_tables/", site, "_anomaliesT_70.csv", sep =""))
  anomaliesT_45 <- read.csv(paste("E:/TACR/WBG_projekce/Ensemble_nc/Climate_tables/", site, "_anomaliesT_45.csv", sep =""))
  anomaliesT_26 <- read.csv(paste("E:/TACR/WBG_projekce/Ensemble_nc/Climate_tables/", site, "_anomaliesT_26.csv", sep =""))
  anomaliesP_85 <- read.csv(paste("E:/TACR/WBG_projekce/Ensemble_nc/Climate_tables/", site, "_anomaliesP_85.csv", sep =""))
  anomaliesP_70 <- read.csv(paste("E:/TACR/WBG_projekce/Ensemble_nc/Climate_tables/", site, "_anomaliesP_70.csv", sep =""))
  anomaliesP_45 <- read.csv(paste("E:/TACR/WBG_projekce/Ensemble_nc/Climate_tables/", site, "_anomaliesP_45.csv", sep =""))
  anomaliesP_26 <- read.csv(paste("E:/TACR/WBG_projekce/Ensemble_nc/Climate_tables/", site, "_anomaliesP_26.csv", sep =""))
  
  
  # Time series of climatic data
  anomaliesT_85_mean <- cbind(year = c(2021:2100), anomaliesT_85[,c(2:13)] + rep(clim_agg$Temp, each = 80))
  anomaliesT_85_dry <- cbind(year = c(2021:2100), anomaliesT_85[,c(2:13)] + rep(clim_agg$Temp_85perc, each = 80))
  anomaliesT_85_cool <- cbind(year = c(2021:2100), anomaliesT_85[,c(2:13)] + rep(clim_agg$Temp_15perc, each = 80))
  anomaliesP_85_mean <- cbind(year = c(2021:2100), anomaliesP_85[,c(2:13)] + rep(clim_agg$Prec, each = 80))
  anomaliesP_85_dry <- cbind(year = c(2021:2100), anomaliesP_85[,c(2:13)] + rep(clim_agg$Prec_25perc, each = 80))
  anomaliesP_85_cool <- cbind(year = c(2021:2100), anomaliesP_85[,c(2:13)] + rep(clim_agg$Prec_75perc, each = 80))
  anomaliesP_85_mean[anomaliesP_85_mean < 0] <- 0
  anomaliesP_85_dry[anomaliesP_85_dry < 0] <- 0
  anomaliesP_85_cool[anomaliesP_85_cool < 0] <- 0
 
  anomaliesT_70_mean <- cbind(year = c(2021:2100), anomaliesT_70[,c(2:13)] + rep(clim_agg$Temp, each = 80))
  anomaliesT_70_dry <- cbind(year = c(2021:2100), anomaliesT_70[,c(2:13)] + rep(clim_agg$Temp_85perc, each = 80))
  anomaliesT_70_cool <- cbind(year = c(2021:2100), anomaliesT_70[,c(2:13)] + rep(clim_agg$Temp_15perc, each = 80))
  anomaliesP_70_mean <- cbind(year = c(2021:2100), anomaliesP_70[,c(2:13)] + rep(clim_agg$Prec, each = 80))
  anomaliesP_70_dry <- cbind(year = c(2021:2100), anomaliesP_70[,c(2:13)] + rep(clim_agg$Prec_25perc, each = 80))
  anomaliesP_70_cool <- cbind(year = c(2021:2100), anomaliesP_70[,c(2:13)] + rep(clim_agg$Prec_75perc, each = 80))
  anomaliesP_70_mean[anomaliesP_70_mean < 0] <- 0
  anomaliesP_70_dry[anomaliesP_70_dry < 0] <- 0
  anomaliesP_70_cool[anomaliesP_70_cool < 0] <- 0
  
  anomaliesT_45_mean <- cbind(year = c(2021:2100), anomaliesT_45[,c(2:13)] + rep(clim_agg$Temp, each = 80))
  anomaliesT_45_dry <- cbind(year = c(2021:2100), anomaliesT_45[,c(2:13)] + rep(clim_agg$Temp_85perc, each = 80))
  anomaliesT_45_cool <- cbind(year = c(2021:2100), anomaliesT_45[,c(2:13)] + rep(clim_agg$Temp_15perc, each = 80))
  anomaliesP_45_mean <- cbind(year = c(2021:2100), anomaliesP_45[,c(2:13)] + rep(clim_agg$Prec, each = 80))
  anomaliesP_45_dry <- cbind(year = c(2021:2100), anomaliesP_45[,c(2:13)] + rep(clim_agg$Prec_25perc, each = 80))
  anomaliesP_45_cool <- cbind(year = c(2021:2100), anomaliesP_45[,c(2:13)] + rep(clim_agg$Prec_75perc, each = 80))
  anomaliesP_45_mean[anomaliesP_45_mean < 0] <- 0
  anomaliesP_45_dry[anomaliesP_45_dry < 0] <- 0
  anomaliesP_45_cool[anomaliesP_45_cool < 0] <- 0
  
  anomaliesT_26_mean <- cbind(year = c(2021:2100), anomaliesT_26[,c(2:13)] + rep(clim_agg$Temp, each = 80))
  anomaliesT_26_dry <- cbind(year = c(2021:2100), anomaliesT_26[,c(2:13)] + rep(clim_agg$Temp_85perc, each = 80))
  anomaliesT_26_cool <- cbind(year = c(2021:2100), anomaliesT_26[,c(2:13)] + rep(clim_agg$Temp_15perc, each = 80))
  anomaliesP_26_mean <- cbind(year = c(2021:2100), anomaliesP_26[,c(2:13)] + rep(clim_agg$Prec, each = 80))
  anomaliesP_26_dry <- cbind(year = c(2021:2100), anomaliesP_26[,c(2:13)] + rep(clim_agg$Prec_25perc, each = 80))
  anomaliesP_26_cool <- cbind(year = c(2021:2100), anomaliesP_26[,c(2:13)] + rep(clim_agg$Prec_75perc, each = 80))
  anomaliesP_26_mean[anomaliesP_26_mean < 0] <- 0
  anomaliesP_26_dry[anomaliesP_26_dry < 0] <- 0
  anomaliesP_26_cool[anomaliesP_26_cool < 0] <- 0
  
  ### VS-Lite model in forcasting mode
  simulation.forecast_85 <- VSLite(phi = phi, Pinput = anomaliesP_85_mean, Tinput = anomaliesT_85_mean,
                                syear = 2021, eyear = 2100,
                                T1 = par$T1, T2 = par$T2, T3 = par$T3, T4 = par$T4,
                                M1 = par$M1, M2 = par$M2, M3 = par$M3, M4 = par$M4,
                                Acor = par$Acor, I_0 = par$I_0,
                                ramp = "modif", integration = "orig",
                                corr = NA)
  
  simulation.forecast_85_dry <- VSLite(phi = phi, Pinput = anomaliesP_85_dry, Tinput = anomaliesT_85_dry,
                                    syear = 2021, eyear = 2100,
                                    T1 = par$T1, T2 = par$T2, T3 = par$T3, T4 = par$T4,
                                    M1 = par$M1, M2 = par$M2, M3 = par$M3, M4 = par$M4,
                                    Acor = par$Acor, I_0 = par$I_0,
                                    ramp = "modif", integration = "orig",
                                    corr = NA)
  
  simulation.forecast_85_cool <- VSLite(phi = phi, Pinput = anomaliesP_85_cool, Tinput = anomaliesT_85_cool,
                                     syear = 2021, eyear = 2100,
                                     T1 = par$T1, T2 = par$T2, T3 = par$T3, T4 = par$T4,
                                     M1 = par$M1, M2 = par$M2, M3 = par$M3, M4 = par$M4,
                                     Acor = par$Acor, I_0 = par$I_0,
                                     ramp = "modif", integration = "orig",
                                     corr = NA)
  
  simulation.forecast_70 <- VSLite(phi = phi, Pinput = anomaliesP_70_mean, Tinput = anomaliesT_70_mean,
                                   syear = 2021, eyear = 2100,
                                   T1 = par$T1, T2 = par$T2, T3 = par$T3, T4 = par$T4,
                                   M1 = par$M1, M2 = par$M2, M3 = par$M3, M4 = par$M4,
                                   Acor = par$Acor, I_0 = par$I_0,
                                   ramp = "modif", integration = "orig",
                                   corr = NA)
  
  simulation.forecast_70_dry <- VSLite(phi = phi, Pinput = anomaliesP_70_dry, Tinput = anomaliesT_70_dry,
                                       syear = 2021, eyear = 2100,
                                       T1 = par$T1, T2 = par$T2, T3 = par$T3, T4 = par$T4,
                                       M1 = par$M1, M2 = par$M2, M3 = par$M3, M4 = par$M4,
                                       Acor = par$Acor, I_0 = par$I_0,
                                       ramp = "modif", integration = "orig",
                                       corr = NA)
  
  simulation.forecast_70_cool <- VSLite(phi = phi, Pinput = anomaliesP_70_cool, Tinput = anomaliesT_70_cool,
                                        syear = 2021, eyear = 2100,
                                        T1 = par$T1, T2 = par$T2, T3 = par$T3, T4 = par$T4,
                                        M1 = par$M1, M2 = par$M2, M3 = par$M3, M4 = par$M4,
                                        Acor = par$Acor, I_0 = par$I_0,
                                        ramp = "modif", integration = "orig",
                                        corr = NA)
  
  simulation.forecast_45 <- VSLite(phi = phi, Pinput = anomaliesP_45_mean, Tinput = anomaliesT_45_mean,
                                   syear = 2021, eyear = 2100,
                                   T1 = par$T1, T2 = par$T2, T3 = par$T3, T4 = par$T4,
                                   M1 = par$M1, M2 = par$M2, M3 = par$M3, M4 = par$M4,
                                   Acor = par$Acor, I_0 = par$I_0,
                                   ramp = "modif", integration = "orig",
                                   corr = NA)
  
  simulation.forecast_45_dry <- VSLite(phi = phi, Pinput = anomaliesP_45_dry, Tinput = anomaliesT_45_dry,
                                       syear = 2021, eyear = 2100,
                                       T1 = par$T1, T2 = par$T2, T3 = par$T3, T4 = par$T4,
                                       M1 = par$M1, M2 = par$M2, M3 = par$M3, M4 = par$M4,
                                       Acor = par$Acor, I_0 = par$I_0,
                                       ramp = "modif", integration = "orig",
                                       corr = NA)
  
  simulation.forecast_45_cool <- VSLite(phi = phi, Pinput = anomaliesP_45_cool, Tinput = anomaliesT_45_cool,
                                        syear = 2021, eyear = 2100,
                                        T1 = par$T1, T2 = par$T2, T3 = par$T3, T4 = par$T4,
                                        M1 = par$M1, M2 = par$M2, M3 = par$M3, M4 = par$M4,
                                        Acor = par$Acor, I_0 = par$I_0,
                                        ramp = "modif", integration = "orig",
                                        corr = NA)

  
  simulation.forecast_26 <- VSLite(phi = phi, Pinput = anomaliesP_26_mean, Tinput = anomaliesT_26_mean,
                                   syear = 2021, eyear = 2100,
                                   T1 = par$T1, T2 = par$T2, T3 = par$T3, T4 = par$T4,
                                   M1 = par$M1, M2 = par$M2, M3 = par$M3, M4 = par$M4,
                                   Acor = par$Acor, I_0 = par$I_0,
                                   ramp = "modif", integration = "orig",
                                   corr = NA)
  
  simulation.forecast_26_dry <- VSLite(phi = phi, Pinput = anomaliesP_26_dry, Tinput = anomaliesT_26_dry,
                                       syear = 2021, eyear = 2100,
                                       T1 = par$T1, T2 = par$T2, T3 = par$T3, T4 = par$T4,
                                       M1 = par$M1, M2 = par$M2, M3 = par$M3, M4 = par$M4,
                                       Acor = par$Acor, I_0 = par$I_0,
                                       ramp = "modif", integration = "orig",
                                       corr = NA)
  
  simulation.forecast_26_cool <- VSLite(phi = phi, Pinput = anomaliesP_26_cool, Tinput = anomaliesT_26_cool,
                                        syear = 2021, eyear = 2100,
                                        T1 = par$T1, T2 = par$T2, T3 = par$T3, T4 = par$T4,
                                        M1 = par$M1, M2 = par$M2, M3 = par$M3, M4 = par$M4,
                                        Acor = par$Acor, I_0 = par$I_0,
                                        ramp = "modif", integration = "orig",
                                        corr = NA)
  ### MODEL - plotting results
  # growth.rates(simulation.forecast_85)
  # growth.rates.cumul(simulation.forecast_85)
  # plot(t(simulation.forecast_85_dry$mod.trw)[,1], type = "l")
  # growth.matrix(simulation.forecast_85_dry)
  # growth.rates(simulation.forecast_85_dry)
  
  saveRDS(simulation.forecast_85, paste("e:/TACR/VS_for_GACR/Results_Forecast8_5/Mean/", genus, "/", site, ".Rda", sep = ""))
  saveRDS(simulation.forecast_85_dry, paste("e:/TACR/VS_for_GACR/Results_Forecast8_5/Dry/", genus, "/", site, ".Rda", sep = ""))
  saveRDS(simulation.forecast_85_cool, paste("e:/TACR/VS_for_GACR/Results_Forecast8_5/Cool/", genus, "/", site, ".Rda", sep = ""))
  saveRDS(simulation.forecast_70, paste("e:/TACR/VS_for_GACR/Results_Forecast7_0/Mean/", genus, "/", site, ".Rda", sep = ""))
  saveRDS(simulation.forecast_70_dry, paste("e:/TACR/VS_for_GACR/Results_Forecast7_0/Dry/", genus, "/", site, ".Rda", sep = ""))
  saveRDS(simulation.forecast_70_cool, paste("e:/TACR/VS_for_GACR/Results_Forecast7_0/Cool/", genus, "/", site, ".Rda", sep = ""))
  saveRDS(simulation.forecast_45, paste("e:/TACR/VS_for_GACR/Results_Forecast4_5/Mean/", genus, "/", site, ".Rda", sep = ""))
  saveRDS(simulation.forecast_45_dry, paste("e:/TACR/VS_for_GACR/Results_Forecast4_5/Dry/", genus, "/", site, ".Rda", sep = ""))
  saveRDS(simulation.forecast_45_cool, paste("e:/TACR/VS_for_GACR/Results_Forecast4_5/Cool/", genus, "/", site, ".Rda", sep = ""))
  saveRDS(simulation.forecast_26, paste("e:/TACR/VS_for_GACR/Results_Forecast2_6/Mean/", genus, "/", site, ".Rda", sep = ""))
  saveRDS(simulation.forecast_26_dry, paste("e:/TACR/VS_for_GACR/Results_Forecast2_6/Dry/", genus, "/", site, ".Rda", sep = ""))
  saveRDS(simulation.forecast_26_cool, paste("e:/TACR/VS_for_GACR/Results_Forecast2_6/Cool/", genus, "/", site, ".Rda", sep = ""))
  
  }
  print(paste(i))
}

