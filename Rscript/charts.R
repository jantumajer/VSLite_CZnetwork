library(reshape2); library(ggplot2)

######################################################################################
### Matrix of monthly integral growth rates with dominant climatic limiting factor ###
######################################################################################

growth.matrix <- function(simul, tun){
  
  ## 1] Extracting key variables from the simulation
  # Growth rates
  Gr <- as.matrix(simul$gINT)
  GrM <- simul$gM
  GrT <- simul$gT
  
  # Climatic variables
  sm <- simul$Moist
  temp <- simul$Temp
  
  # Model parameters
  par <- tun$best
  
  ## 2] Reformating data into matrixes
  matrix <- matrix(0, ncol=ncol(Gr), nrow = 12)
  colnames(matrix) <- c(simul$syear : simul$eyear)
  rownames(matrix) <- c(1:12)

  ## 3] Conditions to identify a dominant limiting factor for each month
  # Different versions depending on shape of ramp function and integration equation
  
  # 3a] Original ramp function and original integration
  if(tun$integration == "orig" & tun$ramp == "orig"){
    matrix[((GrT > GrM) & GrM > 0)] <- "Drought" 
    matrix[((GrT < GrM) & GrT > 0)] <- "Cold" 
    matrix[GrT == 1 & GrM == 1] <- "Optimal" 
    matrix[Gr == 0] <- "Dormancy"
    }
  
  # 3b] Modified ramp function and original integration
  if(tun$integration == "orig" & tun$ramp == "modif"){
    matrix[((GrT > GrM) & GrM > 0 & sm > par[1,"M3"])] <- "Moist" 
    matrix[((GrT > GrM) & GrM > 0 & sm < par[1,"M2"])] <- "Drought" 
    matrix[((GrT < GrM) & GrT > 0 & temp > par[1,"T3"])] <- "Warm" 
    matrix[((GrT < GrM) & GrT > 0 & temp < par[1,"T2"])] <- "Cold" 
    matrix[GrT == 1 & GrM == 1] <- "Optimal" 
    matrix[Gr == 0] <- "Dormancy"
    }
  
  # 3c] Original ramp function and modified integration
  if(tun$integration == "modif" & tun$ramp == "orig"){
    matrix[(GrM > 0 & GrT == 1)] <- "Drought" 
    matrix[(GrT > 0 & GrM == 1)] <- "Cold" 
    matrix[GrT == 1 & GrM == 1] <- "Optimal"
    matrix[Gr > 0 & GrT < 1 & GrM < 1] <- "Mixed" 
    matrix[Gr == 0] <- "Dormancy"
    }
  
  # 3d] Modified ramp function and modified integration
  if(tun$integration == "modif" & tun$ramp == "modif"){
    matrix[(GrT == 1 & GrM > 0 & sm > par[1,"M3"])] <- "Moist" 
    matrix[(GrT == 1 & GrM > 0 & sm < par[1,"M2"])] <- "Drought" 
    matrix[(GrM == 1 & GrT > 0 & temp > par[1,"T3"])] <- "Warm" 
    matrix[(GrM == 1 & GrT > 0 & temp < par[1,"T2"])] <- "Cold" 
    matrix[GrT == 1 & GrM == 1] <- "Optimal" 
    matrix[Gr > 0 & GrT < 1 & GrM < 1] <- "Mixed" 
    matrix[Gr == 0] <- "Dormancy"
    }
  
  ## 4] Reformating (melting)
  matrix.Gr.melt <- melt(Gr)
  matrix.limit.melt <- melt(matrix)
  
  INPUT.TO.PLOT <- cbind(matrix.limit.melt, matrix.Gr.melt$value); colnames(INPUT.TO.PLOT) <- c("MONTH", "YEAR", "Limit", "GrINT")

  ## 5] Plotting
  # Define colors for each climatic growth-limiting factors
  colvec<- c("Drought" = "blue",
             "Cold" = "red",
             "Warm" = "purple",
             "Moist" = "grey20",
             "Optimal" = "green4",
             "Mixed" = "orange2",
             "Dormancy" = "white")

  matrix.plot <- ggplot(data = INPUT.TO.PLOT) + 
    geom_raster(aes(fill = Limit, alpha = GrINT, x = YEAR, y = MONTH)) + 
    scale_fill_manual(values = colvec)+
    scale_x_continuous(name = "Year", limits = c(simul$syear-1, simul$eyear+1))+
    scale_y_continuous(name = "Month", limits = c(0, 13), breaks = c(1:12), labels = c("Jan", "", "Mar", "", "May", "", "Jul", "", "Sep", "", "Nov", ""))+
    theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 8, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          axis.ticks.length=unit(.1, "cm"),
          strip.text = element_text(size = 10, colour = "black"),
          # legend.position = "none",
          # plot.title = element_text(color = "black", hjust = 0.5, vjust = 0, face = "bold", size = 18)
          )

return(matrix.plot)

}

###########################################
### Observed and simulated chronologies ###
###########################################

obsmod.chronologies <- function(simul, tun){
  
  # 1] Merging simulated and observed chronologies from their source tables
  INPUT.TO.PLOT <- as.data.frame(cbind(c(simul$syear : simul$eyear),
                                       t(simul$mod.trw),
                                       (tun$obs.trw - mean(tun$obs.trw))/sd(tun$obs.trw)))
  colnames(INPUT.TO.PLOT) <- c("YEAR", "MODEL", "OBSERVED")
  
  # 2] Correlation coefficient between observed and simulated chronologies
  cor.value <- cor(INPUT.TO.PLOT$MODEL, INPUT.TO.PLOT$OBSERVED)
  
  # 3] Plotting
  chronologies <- ggplot(data = INPUT.TO.PLOT) +
    geom_hline(yintercept = 0, col = "grey20", linetype = "dashed", linewidth = 0.1, alpha = 0.4)+
    geom_line(aes(x = YEAR, y = MODEL), col = "red", linetype = "solid", linewidth = 0.2, alpha = 1)+
    geom_line(aes(x = YEAR, y = OBSERVED), col = "blue", linetype = "solid", linewidth = 0.2, alpha = 1)+
    geom_text(aes(label = round(as.numeric(cor.value), 2)), x = (simul$eyear + simul$syear)/2, y = max(INPUT.TO.PLOT[,c(2,3)]) - 0.25, size = 5)+
    xlim(simul$syear, simul$eyear) +
    xlab("Year")+
    ylab(" ")+
    theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 8, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          axis.ticks.length=unit(.3, "cm"),
          strip.text.y = element_text(size = 10, colour = "black"),
          strip.text.x = element_text(size = 10, colour = "black"),  
          # legend.position = "none",
          plot.title = element_text(color = "black", hjust = 0.5, vjust = 0, face = "bold", size = 18))

return(chronologies)
}

##########################################
### Intra-annual chart of growth rates ###
##########################################

growth.rates <- function(simul){
  
  ## 1] Extracting key variables from the simulation
  # Growth rates
  Gr <- simul$gINT
  GrM <- simul$gM
  GrT <- simul$gT
  GrE <- simul$gE
  
  # 2] Calculating long-term means of partial (GrT, GrM) and integral (GrINT) growth rates
  res <- data.frame(MONTH = c(1:12), Temp = NA, Moist = NA, Sol = NA, Int=NA, Intplus = NA, Intminus = NA)
  for (i in c(1:12)){
    res[i,"Temp"] <- mean(GrT[i,])
    res[i,"Moist"] <- mean(GrM[i,])
    res[i,"Int"] <- mean(Gr[i,])
    res[i,"Sol"] <- GrE[i,1]
    res[i,"Intminus"] <- max(0, mean(Gr[i,]) - sd(Gr[i,]))
    res[i,"Intplus"] <- min(1, mean(Gr[i,]) + sd(Gr[i,]))
  }

  clim <- ggplot(data=res, aes(x = MONTH)) +
    # geom_ribbon(aes(ymin=Intminus, ymax=Intplus), alpha = 0.25, fill = "black")+
    geom_line(aes(y=Sol), col = "orange", linetype = "solid", linewidth = 0.25)+
    geom_line(aes(y=Temp), col = "red", linetype = "solid", linewidth = 0.25)+
    geom_line(aes(y=Moist), col = "blue", linetype = "solid", linewidth = 0.25)+
    geom_line(aes(y=Int), col = "black", linetype = "solid", linewidth = 0.5)+
    ylim(0,1) +
    scale_x_continuous(name = "Month", limits = c(1, 12), breaks = c(1:12), labels = c("Jan", "", "Mar", "", "May", "", "Jul", "", "Sep", "", "Nov", ""))+
    ylab("Partial and integral growth rates")+
    scale_color_grey()+ 
    theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, colour = "black"),
          axis.ticks.length=unit(.3, "cm"),
          strip.text = element_text(size = 14, colour = "black"),
          # legend.position = "none",
          plot.title = element_text(color = "black", hjust = 0.5, vjust = 0, face = "bold", size = 18))

  return(clim)

}

############################################
### Cumulative growth rate over the year ###
############################################

growth.rates.cumul <- function(simul){
  
  ## 1] Extracting key variables from the simulation
  # Growth rates
  Gr <- simul$gINT

  # 2] Calculating cumulative integral growth rates (GrINT)
  Gr.cumul <- Gr
  for (i in c(2:12)){
    Gr.cumul[i,] <- Gr.cumul[(i-1), ] + Gr[i, ]
  }
  
  Gr.cumul.melt <- melt(Gr.cumul) # Individual years
  Gr.cumul.mean <- data.frame(month = c(1:12), value = rowMeans(Gr.cumul))
  
  clim <- ggplot() +
    geom_line(aes(x = Var1, y = value, group = Var2), col = "grey20", alpha = 0.25, linetype = "solid", linewidth = 0.5, data = Gr.cumul.melt) +
    geom_line(aes(x = month, y = value), col = "red", linetype = "solid", linewidth = 1.5, data = Gr.cumul.mean) +
    scale_x_continuous(name = "Month", limits = c(1, 12), breaks = c(1:12), labels = c("Jan", "", "Mar", "", "May", "", "Jul", "", "Sep", "", "Nov", ""))+
    ylab("Cumulative integral growth rates")+
    scale_color_grey()+ 
    theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 8, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          axis.ticks.length=unit(.3, "cm"),
          strip.text = element_text(size = 10, colour = "black"),
          # legend.position = "none",
          plot.title = element_text(color = "black", hjust = 0.5, vjust = 0, face = "bold", size = 18))
  
  return(clim)
  
}

#################
### Phenology ###
#################

phenology <- function(simul){
  
  ## 1] Extracting key variables from the simulation
  # Growth rates
  Gr <- simul$gINT
  rownames(Gr) <- c(1:12)
  colnames(Gr) <- c(simul$syear : simul$eyear)
  
  INPUT.TO.PLOT <- data.frame(YEAR = c(simul$syear : simul$eyear),
                              SOS = NA, # start of growing season
                              EOS = NA) # end of growing season
  
  # 2] Identifying the first (SOS) and the last (EOS) month with GrINT > 0 
  for (j in c(1:ncol(Gr))){
    if(sum(Gr[,j] == 0) == 12) {INPUT.TO.PLOT[j, c("SOS", "EOS")] <- NA} # No growth
    if(sum(Gr[,j] == 0) == 11) {a <- as.data.frame(Gr[,j] == 0);
                                k <- 2; while(a[k,1] == a[k+1,1]){k <- k + 1}
                                INPUT.TO.PLOT[j, c("SOS", "EOS")] <- k + 1}
    if(sum(Gr[,j] == 0) < 11) {INPUT.TO.PLOT[j, "SOS"] <- min(as.numeric(rownames(Gr[!(Gr[,j] == 0),])))
                               INPUT.TO.PLOT[j, "EOS"] <- max(as.numeric(rownames(Gr[!(Gr[,j] == 0),])))} 
  }
  
  # 3] Plotting
  phenology <- ggplot(data = INPUT.TO.PLOT) +
    geom_line(aes(x = YEAR, y = SOS), col = "green", linetype = "solid", linewidth = 0.2, alpha = 1)+
    geom_line(aes(x = YEAR, y = EOS), col = "brown", linetype = "solid", linewidth = 0.2, alpha = 1)+
    xlim(simul$syear, simul$eyear) +
    scale_y_continuous(name = "Month", limits = c(1, 12), breaks = c(1:12), labels = c("Jan", "", "Mar", "", "May", "", "Jul", "", "Sep", "", "Nov", ""))+
    xlab("Year")+
    theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 8, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          axis.ticks.length=unit(.3, "cm"),
          strip.text.y = element_text(size = 10, colour = "black"),
          strip.text.x = element_text(size = 10, colour = "black"),  
          # legend.position = "none",
          plot.title = element_text(color = "black", hjust = 0.5, vjust = 0, face = "bold", size = 18))
  
  return(phenology)
}


##############################
### Trends in growth rates ###
##############################
growth.rates.trends <- function(simul){
  
  ## 1] Extracting key variables from the simulation
  # Growth rates
  Gr <- simul$gINT
  GrM <- simul$gM
  GrT <- simul$gT
  GrE <- simul$gE

  # 1b] Output file
  Output <- data.frame(MONTH = c(1:12), Temp.Slope = NA, Temp.P = NA, Moist.Slope = NA, Moist.P = NA, Int.Slope = NA, Int.P = NA)

  # 2a] Calculation of trends
  for (i in c(1:12)){
    DATA <- data.frame(YEAR=c(simul$syear : simul$eyear), T=(GrT[i,]), M=(GrM[i,]), GR=(Gr[i,])); colnames(DATA) <- c("YEAR", "T", "M", "GR")
  
    modelT <- summary(lm(DATA$T ~ DATA$YEAR))
    Output[i, "Temp.Slope"] <- modelT$coefficients[2,1] # p-value of slope
    Output[i, "Temp.P"] <- modelT$coefficients[2,4] # Estimate of slope

    modelM <- summary(lm(DATA$M ~ DATA$YEAR))
    Output[i, "Moist.Slope"] <- modelM$coefficients[2,1] # p-value of slope
    Output[i, "Moist.P"] <- modelM$coefficients[2,4] # Estimate of slope

    modelGR <- summary(lm(DATA$GR ~ DATA$YEAR))
    Output[i, "Int.Slope"] <- modelGR$coefficients[2,1] # p-value of slope
    Output[i, "Int.P"] <- modelGR$coefficients[2,4] # Estimate of slope

  
    rm(modelM, modelT, modelGR, DATA)
  
  }
  
  Temp.dat <- Output[,c(1:3)]
  Moist.dat <- Output[,c(1, 4, 5)]
  Int.dat <- Output[,c(1, 6, 7)]
  colnames(Temp.dat) <- colnames(Moist.dat) <- colnames(Int.dat) <- c("Month", "Slope", "P")
  Output.melt <- rbind(cbind(Temp.dat, VAR = "Temperature (GrT)"),
                       cbind(Moist.dat, VAR = "Moisture (GrM)"),
                       cbind(Int.dat, VAR = "Integral (GrINT)"))
  
  Output.melt[is.nan(Output.melt$P), "P"] <- 1
  Output.melt[Output.melt$P < 0.05, "Sig"] <- 1
  Output.melt[Output.melt$P >= 0.05, "Sig"] <- 0.25
  
  trends <- ggplot(data=Output.melt) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_col(aes(x = Month, y = Slope, alpha = Sig, fill = VAR))+
    scale_x_continuous(name = "Month", limits = c(0, 13), breaks = c(1:12), labels = c("Jan", "", "Mar", "", "May", "", "Jul", "", "Sep", "", "Nov", ""))+
    ylab("Trend in partial and integral growth rates [1/year]")+
    scale_fill_manual(values = c("Temperature (GrT)" = "red", "Moisture (GrM)" = "blue", "Integral (GrINT)" = "black"))+ 
    facet_grid(VAR ~ .) +
    theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 8, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          axis.ticks.length=unit(.3, "cm"),
          strip.text = element_text(size = 10, colour = "black"),
          legend.position = "none",
          plot.title = element_text(color = "black", hjust = 0.5, vjust = 0, face = "bold", size = 18))
  
  return(trends)
}