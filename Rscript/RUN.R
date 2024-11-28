library(dplR)
library(ggplot2)
library(reshape2)
library(readxl)


###########################################
#### 0] Loading all required functions ####
###########################################
setwd("E:/VSLite_R/R/")

## 0a] Model sub-algorithms ##
# Calculation of partial growth rates to photoperiod
source("compute.gE.R") # Source Suzan Tolwinski-Ward
source("daylength.factor.from.lat.R") # Source Suzan Tolwinski-Ward

# Soil moisture model
source("leakybucket.monthly.R") # Source Suzan Tolwinski-Ward 
source("leakybucket.submonthly.R") # Source Suzan Tolwinski-Ward

# Ramp functions
source("std.ramp.R") # Source Suzan Tolwinski-Ward - Original non-declining ramp functions
source("mod.ramp.R") # Modified increasing-stable-decreasing ramp functions

# Integration functions
source("integrate.orig.R") # Source Suzan Tolwinski-Ward - Integration based on MINIMUM of growth rates (following Liebig's law)
source("integrate.multiplic.R") # Integration based on PRODUCT of growth rates (following initial TRACH model)

# Model definitions
source("VSLite.R") # Source Suzan Tolwinski-Ward

## 0b] Functions to calibrate the model against local site chronology ##
source("randomization.R")
source("VSLite.iterative.R")

## 0c] Graphical functions
source("charts.R")

## 0d] detrending climatic data
source("climate.detrend.R")

#########################################
### 0] Loading the list of input data ###
#########################################
SITE <- read.csv("E:/TACR/zchron_TestData2/META/data_april2_withCZU.csv", row.names = 1)
ADD <- read_excel("E:/TACR/zchron_TestData2/ADD/additions_202303.xlsx", sheet = "site")
ADD2 <- read_excel("E:/TACR/zchron_TestData2/ADD2/additions_202401.xlsx", sheet = "site")
ADD_Senfeldr <- read.csv("E:/TACR/zchron_TestData2/ADD_Senfeldr/add_senfeldr.csv")
ADD_Mendelu <- read_excel("E:/TACR/zchron_TestData2/ADD_Mendelu/Mendelu_site_data_nove2023.xlsx")
ADD_IFER <- read_excel("E:/TACR/VS_for_GACR/IFER/metadata_IFER.xlsx")
SITE <- rbind(SITE, ADD[,c(2:50)], ADD2[,c(2:50)], ADD_Senfeldr, ADD_Mendelu, ADD_IFER)
rm(ADD, ADD2, ADD_Senfeldr, ADD_Mendelu, ADD_IFER)
SITE <- SITE[!(SITE$raw_data_file_name == ""),] # Series from Cada and Rydval still missing


### Spline and chronology variant
nyr = 60
variant = "std"

for (i in c(1:nrow(SITE))){
  site <- SITE[i, "site_code"]
  rwl.code <- SITE[i, "raw_data_file_name"]
  genus <- substr(SITE[i, "species"], 1, 2)
  
  ### Loading geographical data
  phi <- SITE[i, "site_lat_decimal"]
  
  ### Loading and preprocessing climatic data
  clim <- read.csv(paste("E:/TACR/Climate_Grids/Climate_tables/", site, "_clim.csv", sep = ""))
  temperature <- dcast(clim, formula = Year ~ Month, value.var = "Temp")
  precipitation <- dcast(clim, formula = Year ~ Month, value.var = "Prec")
  # Detrending
  # temperature.det <- climate.detrend(temperature, var = "temp", spline = nyr)
  # precipitation.det <- climate.detrend(precipitation, var = "prec", spline = nyr)
  
  ### Loading and preprocessing growth data
  series <- read.rwl(paste("E:/TACR/zchron_TestData2/RWL/", rwl.code, sep = ""))
  det.series <- detrend(series, method = "Spline", nyrs = nyr)
  chronology <- chron(det.series, biweight = T, prewhiten = T)
  # Period
  chronology.sub <- chronology[chronology$samp.depth > 4, ]
  if (nrow(chronology.sub) > 40 & "1961" %in% row.names(chronology.sub)){
  syear <- 1961
  eyear <- min(max(as.numeric(row.names(chronology.sub))), 2020)

  ### MODEL - randomization of parameters
  random.par <- randomization(iter = 10000)
  tuning <- VSLite.iterative(parameters = random.par, obs.chron = chronology.sub, chron.variant = variant,
                             syear = syear, eyear = eyear,
                             phi = phi, Pinput = precipitation, Tinput = temperature,
                             ramp = "modif", integration = "orig")
  
  if(max(tuning$best$correlation) > -1){
  
    ### MODEL - simulation based on the optimal parameters
    simulation <- VSLite(phi = phi, Pinput = precipitation, Tinput = temperature,
                        syear = syear, eyear = eyear,
                        T1 = tuning$best$T1, T2 = tuning$best$T2, T3 = tuning$best$T3, T4 = tuning$best$T4,
                        M1 = tuning$best$M1, M2 = tuning$best$M2, M3 = tuning$best$M3, M4 = tuning$best$M4,
                        Acor = tuning$best$Acor, I_0 = tuning$best$I_0,
                        ramp  = tuning$ramp, integration = tuning$integration,
                        corr = tuning$best$correlation)
  
    ### MODEL - plotting results
    # growth.matrix(simulation, tuning)
    # obsmod.chronologies(simulation, tuning)
    # growth.rates(simulation)
    # growth.rates.trends(simulation)
    # growth.rates.cumul(simulation)
  
    saveRDS(simulation, paste("e:/TACR/VS_for_GACR/Results/", genus,"/", site, ".Rda", sep = ""))
    }
  }
}

