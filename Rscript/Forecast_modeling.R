library(RColorBrewer); library(plot.matrix); library(ggplot2); library(reshape2); library(readxl)
library(e1071); library(multimode); library(pracma)

#################

# Model definitions
setwd("E:/VSLite_R/R/")

## 0a] Model sub-algorithms ##
# Calculation of partial growth rates to photoperiod
source("compute.gE.R") 
source("daylength.factor.from.lat.R") 

# Soil moisture model
source("leakybucket.monthly.R") 
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
source("randomization.R")
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

# SSP5-8.5
TTrend <- read.delim("E:/TACR/WBG_projekce/T8_5.txt")
PTrend <- read.delim("E:/TACR/WBG_projekce/P8_5.txt")
ssp <- "8_5"

# SSP3-7.0
TTrend <- read.delim("E:/TACR/WBG_projekce/T7_0.txt")
PTrend <- read.delim("E:/TACR/WBG_projekce/P7_0.txt")
ssp <- "7_0"

# SSP2-4.5
TTrend <- read.delim("E:/TACR/WBG_projekce/T4_5.txt")
PTrend <- read.delim("E:/TACR/WBG_projekce/P4_5.txt")
ssp <- "4_5"

# SSP1-2.6
TTrend <- read.delim("E:/TACR/WBG_projekce/T2_6.txt")
PTrend <- read.delim("E:/TACR/WBG_projekce/P2_6.txt")
ssp <- "2_6"

nyr <- 60

#################

for (i in c(1:nrow(SITE))){
  site <- SITE[i, "site_code"] # Which site is processed
  genus <- substr(SITE[i, "species"], 1, 2)
  phi <- SITE[i, "site_lat_decimal"]
  
  if(paste(site, ".Rda", sep ="") %in% list.files(paste("e:/TACR/VS_for_GACR/Results/", genus, sep = ""))){
  
  # Model parameters
  model <- readRDS(paste("e:/TACR/VS_for_GACR/Results/", genus,"/", site, ".Rda", sep = ""))
  par <- model$par

  ### Loading and preprocessing climatic data
  clim <- read.csv(paste("E:/TACR/Climate_Grids/Climate_tables/", site, "_clim.csv", sep = ""))
  clim_sub <- clim[clim$Year %in% c(1995:2014),]
  
  # Detrending climatic data
  # temperature <- dcast(clim, formula = Year ~ Month, value.var = "Temp")
  # precipitation <- dcast(clim, formula = Year ~ Month, value.var = "Prec")
  # temperature.det <- climate.detrend(temperature, var = "temp", spline = nyr)
  # precipitation.det <- climate.detrend(precipitation, var = "prec", spline = nyr)
  
  clim_agg <- aggregate(clim_sub[,c("Temp", "Prec")], by = list(Month = clim_sub$Month), FUN = mean)
  # clim_agg <- data.frame(Month = c(1:12), Temp = colMeans(temperature.det[temperature.det$Year %in% c(1995:2014),c(2:13)]),
  #                                         Prec = colMeans(precipitation.det[precipitation.det$Year %in% c(1995:2014),c(2:13)]))
  
  Tbaseline <- as.data.frame(rbind(clim_agg$Temp, clim_agg$Temp, clim_agg$Temp, clim_agg$Temp, clim_agg$Temp))
  Pbaseline <- as.data.frame(rbind(clim_agg$Prec, clim_agg$Prec, clim_agg$Prec, clim_agg$Prec, clim_agg$Prec))
  Tbaseline <- cbind(YEAR = c(1:5), Tbaseline)
  Pbaseline <- cbind(YEAR = c(1:5), Pbaseline)
  
  Tforecast <- Pforecast <- TTrend
  

  for (m in c(1:nrow(Tforecast))){
    for (n in c(3:14)){
    Tforecast[m, n] <- Tbaseline[1,(n-1)] + TTrend[m, n]
    Pforecast[m, n] <- Pbaseline[1,(n-1)] + PTrend[m, n]
    }
  }
  
  Pforecast[Pforecast < 0] <- 0
  
  ### VS-Lite model in forcasting mode
  simulation.forecast <- VSLite(phi = phi, Pinput = Pforecast[,c(2:14)], Tinput = Tforecast[,c(2:14)],
                                syear = 2021, eyear = 2100,
                                T1 = par$T1, T2 = par$T2, T3 = par$T3, T4 = par$T4,
                                M1 = par$M1, M2 = par$M2, M3 = par$M3, M4 = par$M4,
                                Acor = par$Acor, I_0 = par$I_0,
                                ramp = "modif", integration = "orig",
                                corr = NA)

  
  ### MODEL - plotting results
  # growth.rates(simulation.forecast)
  # growth.rates.cumul(simulation.forecast)

  
  saveRDS(simulation.forecast, paste("e:/TACR/VS_for_GACR/Results_Forecast", ssp, "/", genus, "/", site, ".Rda", sep = ""))

  
  }
  print(paste(i))
}

