## Water PCB concentrations data analysis per site
## Passaic River

# Install packages
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("robustbase")
install.packages("dplyr")
install.packages("tibble")
install.packages("Matrix")
install.packages("lme4")
install.packages("MuMIn")
install.packages("lmerTest")
install.packages("zoo")
install.packages("dataRetrieval")
install.packages("reshape")
install.packages("tidyr")
install.packages('patchwork')
install.packages("scales")

# Load libraries
{
  library(ggplot2)
  library(scales) # function trans_breaks
  library(stringr) # str_detect
  library(robustbase) # function colMedians
  library(dplyr) # performs %>%
  library(tibble) # adds a column
  library(lme4) # performs lme
  library(MuMIn) # gets Rs from lme
  library(lmerTest) # gets the p-value from lme
  library(zoo) # yields seasons
  library(dataRetrieval) # read data from USGS
  library(reshape)
  library(tidyr) # function gather
  library(patchwork) # combine plots
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# Select Passaic River data ---------------------------------------------------
pas <- wdc[str_detect(wdc$LocationName, 'Passaic River'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  pas$SampleDate <- as.Date(pas$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(pas$SampleDate) - min(as.Date(pas$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- pas$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(pas$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  pas.tpcb <- cbind(factor(pas$SiteID), pas$SampleDate,
                    pas$Latitude, pas$Longitude, as.matrix(pas$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(pas.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Include USGS flow and temperature data --------------------------------------------------
{
  # Include flow data from USGS station Passaic River
  sitepasN1 <- "01381900" # No temp
  sitepasN2 <- "01379500" # No temp
  sitepasN3 <- "01389005" # No flow
  sitepasN4 <- "01389010" # No temp
  sitepasN5 <- "01389500" # No temp
  sitepasN6 <- "01389890" # No temp

  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow.1 <- readNWISdv(sitepasN1, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  flow.2 <- readNWISdv(sitepasN2, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  flow.3 <- readNWISdv(sitepasN4, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  flow.4 <- readNWISdv(sitepasN5, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  flow.5 <- readNWISdv(sitepasN6, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  temp <- readNWISdv(sitepasN3, paramtemp,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  
  # Add USGS data to pass.tpcb.2, matching dates, conversion to m3/s
  pas.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.1$Date)]
  pas.tpcb$flow.2 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.2$Date)]
  pas.tpcb$flow.3 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.3$Date)]
  pas.tpcb$flow.4 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.4$Date)]
  pas.tpcb$flow.5 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.5$Date)]
  pas.tpcb$temp <- 273.15 + temp$X_.from.middle.intake_00010_00003[match(pas.tpcb$date,
                                                       temp$Date)]
}

# Remove site -------------------------------------------------------------
# Remove site located in the ocean.Possible typo in original coordinates.
pas.tpcb.1 <- subset(pas.tpcb, SiteID != c("WCPCB-PASS022"))

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- pas.tpcb$tPCB
time <- pas.tpcb$time
site <- pas.tpcb$site.code
season <- pas.tpcb$season
flow <- pas.tpcb$flow.5
tem <- pas.tpcb$temp
# tPCB vs. time + season + flow + temp + site
lme.pas.tpcb <- lmer(log10(tpcb) ~ 1 + time + flow + tem + season + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"))

# See results
summary(lme.pas.tpcb)
# Look at residuals
{
  res.pas.tpcb <- resid(lme.pas.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/PassaicRiverQ-QtPCB.pdf")
  qqnorm(res.pas.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.pas.tpcb)
  dev.off()
}
# Shapiro test
shapiro.test(resid(lme.pas.tpcb)) # Lme doesn't work.

# LME for individual PCBs -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  pas.pcb <- subset(pas, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  pas.pcb <- subset(pas.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  pas.pcb <- log10(pas.pcb)
  # Replace -inf to NA
  pas.pcb <- do.call(data.frame,
                     lapply(pas.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  pas.pcb.1 <- pas.pcb[,
                       -which(colSums(is.na(pas.pcb))/nrow(pas.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(pas$SiteID)
  # Change date format
  SampleDate <- as.Date(pas$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Create individual code for each site sampled
  site.numb <- pas$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(pas$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to pass.pcb
  pas.pcb.1 <- cbind(pas.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     site.numb, season.s)
  # Include flow data from USGS station Passaic River
  sitepasN1 <- "01381900" # No temp
  sitepasN2 <- "01379500" # No temp
  sitepasN3 <- "01389005" # No flow
  sitepasN4 <- "01389010" # No temp
  sitepasN5 <- "01389500" # No temp
  sitepasN6 <- "01389890" # No temp
  
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow.1 <- readNWISdv(sitepasN1, paramflow,
                       min(pas.pcb.1$SampleDate), max(pas.pcb.1$SampleDate))
  flow.2 <- readNWISdv(sitepasN2, paramflow,
                       min(pas.pcb.1$SampleDate), max(pas.pcb.1$SampleDate))
  flow.3 <- readNWISdv(sitepasN4, paramflow,
                       min(pas.pcb.1$SampleDate), max(pas.pcb.1$SampleDate))
  flow.4 <- readNWISdv(sitepasN5, paramflow,
                       min(pas.pcb.1$SampleDate), max(pas.pcb.1$SampleDate))
  flow.5 <- readNWISdv(sitepasN6, paramflow,
                       min(pas.pcb.1$SampleDate), max(pas.pcb.1$SampleDate))
  temp <- readNWISdv(sitepasN3, paramtemp,
                     min(pas.pcb.1$SampleDate), max(pas.pcb.1$SampleDate))
  
  # Add USGS data to pass.tpcb.1, matching dates, conversion to m3/s
  pas.pcb.1$flow.1 <- 0.03*flow.1$X_00060_00003[match(pas.pcb.1$SampleDate,
                                                      flow.1$Date)]
  pas.pcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(pas.pcb.1$SampleDate,
                                                      flow.2$Date)]
  pas.pcb.1$flow.3 <- 0.03*flow.3$X_00060_00003[match(pas.pcb.1$SampleDate,
                                                      flow.3$Date)]
  pas.pcb.1$flow.4 <- 0.03*flow.4$X_00060_00003[match(pas.pcb.1$SampleDate,
                                                      flow.4$Date)]
  pas.pcb.1$flow.5 <- 0.03*flow.5$X_00060_00003[match(pas.pcb.1$SampleDate,
                                                      flow.5$Date)]
  pas.pcb.1$temp <- 273.15 + temp$X_.from.middle.intake_00010_00003[match(pas.pcb.1$SampleDate,
                                                                          temp$Date)]
  # Remove metadata
  pas.pcb.2 <- subset(pas.pcb.1, select = -c(SiteID:temp))
}

# Get covariates
time <- pas.pcb.1$time
season <- pas.pcb.1$season
site <- pas.pcb.1$site.numb
flow <- pas.pcb.1$flow.3 # flow.3 yiels 6 PCB congeners
tem <- pas.pcb.1$temp

# Create matrix to store results
lme.pcb <- matrix(nrow = length(pas.pcb.2[1,]), ncol = 25)

# Perform LME
# Need to remove season 3, not enough data.
for (i in 1:length(pas.pcb.2[1,])) {
  fit <- lmer(pas.pcb.2[,i] ~ 1 + time + flow + tem + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"))
  lme.pcb[i,1] <- fixef(fit)[1] # intercept
  lme.pcb[i,2] <- summary(fit)$coef[1,"Std. Error"] # intercept error
  lme.pcb[i,3] <- summary(fit)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.pcb[i,4] <- fixef(fit)[2] # time
  lme.pcb[i,5] <- summary(fit)$coef[2,"Std. Error"] # time error
  lme.pcb[i,6] <- summary(fit)$coef[2,"Pr(>|t|)"] # time p-value
  lme.pcb[i,7] <- fixef(fit)[3] # flow
  lme.pcb[i,8] <- summary(fit)$coef[3,"Std. Error"] # flow error
  lme.pcb[i,9] <- summary(fit)$coef[3,"Pr(>|t|)"] # flow p-value
  lme.pcb[i,10] <- fixef(fit)[4] # temp
  lme.pcb[i,11] <- summary(fit)$coef[4,"Std. Error"] # temp error
  lme.pcb[i,12] <- summary(fit)$coef[4,"Pr(>|t|)"] # temp p-value
  lme.pcb[i,13] <- fixef(fit)[4] # # season 1
  lme.pcb[i,14] <- summary(fit)$coef[5,"Std. Error"] # season 1 error
  lme.pcb[i,15] <- summary(fit)$coef[5,"Pr(>|t|)"] # # season 1 p-value
  lme.pcb[i,16] <- fixef(fit)[6] # season 2
  lme.pcb[i,17] <- summary(fit)$coef[6,"Std. Error"] # season 2 error
  lme.pcb[i,18] <- summary(fit)$coef[6,"Pr(>|t|)"] # season 2 p-value
  lme.pcb[i,19] <- -log(2)/lme.pcb[i,4]/365 # t0.5
  lme.pcb[i,20] <- abs(-log(2)/lme.pcb[i,4]/365)*lme.pcb[i,5]/abs(lme.pcb[i,4]) # t0.5 error
  lme.pcb[i,21] <- as.data.frame(VarCorr(fit))[1,'sdcor']
  lme.pcb[i,22] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
  lme.pcb[i,23] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
  lme.pcb[i,24] <- shapiro.test(resid(fit))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(fit)
  # Calculate residuals and RMSE
  residuals <- pas.pcb.2[, i] - predictions
  non_na_indices <- !is.na(residuals)
  lme.pcb[i, 25] <- sqrt(mean(residuals[non_na_indices]^2))
  
}

# Transform result to data.frame so factor 2 can be included
lme.pcb <- as.data.frame(lme.pcb)

# Add factor of 2
# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(pas.pcb.2[,1]),
                      ncol = length(pas.pcb.2[1,]))

# Create a vector to store factor 2 for each congener
factor2_vector <- numeric(length = length(pas.pcb.2[1,]))

for (i in 1:length(pas.pcb.2[1,])) {
  fit <- lmer(pas.pcb.2[,i] ~ 1 + time + flow + tem + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  
  lme.fit.pcb[,i] <- fitted(fit)
  
  # Calculate factor2 for each congener
  factor2 <- 10^(lme.fit.pcb[, i])/10^(pas.pcb.2[, i])
  factor2_vector[i] <- sum(factor2 > 0.5 & factor2 < 2,
                           na.rm = TRUE) / (sum(!is.na(factor2))) * 100
}

# Add factor 2 to lme.pcb data.frame
lme.pcb$factor2 <- factor2_vector

# Change number format of factor 2 to 3 significant figures
lme.pcb$factor2 <- formatC(signif(lme.pcb$factor2, digits = 3))

# Add congener names
congeners <- colnames(pas.pcb.2)
lme.pcb <- as.data.frame(cbind(congeners, lme.pcb))

# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "flow", "flow.error", "flow.pv", "temp", "temp.error",
                       "temp.pv", "season1", "season1.error", "season1.pv",
                       "season2", "season2.error", "season2, pv", "t05",
                       "t05.error", "RandonEffectSiteStdDev", "R2nR", "R2R",
                       "Normality", "RMSE", "Factor2")

# Remove congeners with no normal distribution
# Shapiro test p-value < 0.05
lme.pcb$Normality <- as.numeric(lme.pcb$Normality)
# Get the congeners that are not showing normality
lme.pcb.out <- lme.pcb[lme.pcb$Normality < 0.05, ]
lme.pcb <- lme.pcb[lme.pcb$Normality > 0.05, ]

# Export results
write.csv(lme.pcb, file = "Output/Data/Sites/csv/PassaicRiver/PassaicLmePCB.csv",
          row.names = FALSE)

# Obtain observations vs predictions
# Select congeners that are not showing normality to be remove from pass.pcb.2
df <- data.frame(names_to_remove = lme.pcb.out$Congeners)
# Get column indices to remove
cols_to_remove <- which(names(pas.pcb.2) %in% df$names_to_remove)
# Remove columns from che.pcb.2 with congeners that don't show normality
pas.pcb.3 <- pas.pcb.2[, -cols_to_remove]

# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(pas.pcb.3[,1]),
                      ncol = length(pas.pcb.3[1,]))

for (i in 1:length(pas.pcb.3[1,])) {
  fit <- lmer(pas.pcb.3[,i] ~ 1 + time + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  lme.fit.pcb[,i] <- fitted(fit)
}

# Individual PCB congener plots -------------------------------------------
# (1) Plot 1:1 for all congeners
# Transform lme.fit.pcb to data.frame
lme.fit.pcb <- as.data.frame(lme.fit.pcb)
# Add congener names to lme.fit.pcb columns
colnames(lme.fit.pcb) <- colnames(pas.pcb.3)
# Add code number to first column
df1 <- cbind(code = row.names(pas.pcb.3), pas.pcb.3)
df2 <- cbind(code = row.names(lme.fit.pcb), lme.fit.pcb)

for (i in 2:length(df1)) {
  col_name <- if (i == 1) {
    ""  # leave the name empty for the first plot
  } else {
    names(df1)[i] # use the column name for other plots
  }
  
  # create plot for each pair of columns
  p <- ggplot(data = data.frame(x = df1$code, y1 = 10^(df1[, i]), y2 = 10^(df2[, i])),
              aes(x = y1, y = y2)) +
    geom_point(shape = 21, size = 3, fill = "white") +
    scale_y_log10(limits = c(0.01, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(limits = c(0.01, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
    ylab(expression(bold("Predicted lme concentration PCBi (pg/L)"))) +
    theme_bw() +
    theme(aspect.ratio = 15/15) +
    annotation_logticks(sides = "bl") +
    geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
    geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
    geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
    annotate('text', x = 10^1, y = 10^4, label = gsub("\\.", "+", names(df1)[i]),
             size = 3, fontface = 2)
  # save plot
  ggsave(paste0("Output/Plots/Sites/ObsPred/PassaicRiver/", col_name, ".png"),
         plot = p, width = 6, height = 6, dpi = 500)
}

# (2) All plots in one page
# Create a list to store all the plots
plot_list <- list()

# Loop over the columns of df1 and df2
for (i in 2:length(df1)) {
  col_name <- paste(names(df1)[i], sep = "")  # use the column name for plot title
  # Create plot for each pair of columns and add to plot_list
  p <- ggplot(data = data.frame(x = df1$code, y1 = 10^(df1[, i]), y2 = 10^(df2[, i])),
              aes(x = y1, y = y2)) +
    geom_point(shape = 21, size = 3, fill = "white") +
    scale_y_log10(limits = c(0.01, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(limits = c(0.01, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
    ylab(expression(bold("Predicted lme concentration PCBi (pg/L)"))) +
    theme_bw() +
    theme(aspect.ratio = 15/15, 
          axis.title = element_text(size = 8)) +
    annotation_logticks(sides = "bl") +
    annotate('text', x = 25, y = 10^4, label = gsub("\\.", "+", col_name),
             size = 2.5, fontface = 2) +
    geom_abline(intercept = 0, slope = 1, col = "white", linewidth = 0.7) +
    geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
    geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7)
  
  plot_list[[i-1]] <- p  # add plot to list
}
# Combine all the plots using patchwork
combined_plot <- wrap_plots(plotlist = plot_list, ncol = 4)
# Save the combined plot
ggsave("Output/Plots/Sites/ObsPred/PassaicRiver/LmeCombined_plot.png", combined_plot,
       width = 15, height = 15, dpi = 500)

# (3) Create a list to store all the cleaned data frames
cleaned_df_list <- list()
# Loop over the columns of df1 and df2
for (i in 2:length(df1)) {
  # Create a new data frame by binding the columns of df1 and df2 for each pair of columns
  df_pair <- cbind(df1[,1], df1[,i], df2[,i])
  colnames(df_pair) <- c("code", "observed", "predicted")
  # Remove the rows with missing values
  cleaned_df_pair <- na.omit(df_pair)
  # Add the cleaned data frame to the list
  cleaned_df_list[[i-1]] <- cleaned_df_pair
}

{
  # Modify data to be plotted
  # Combine all the cleaned data frames using rbind
  combined_cleaned_df <- do.call(rbind, cleaned_df_list)
  # Convert the matrix to a data frame
  combined_cleaned_df <- as.data.frame(combined_cleaned_df)
  # Convert the code column to a factor
  combined_cleaned_df$code <- as.factor(combined_cleaned_df$code)
  # Convert the observed and predicted columns to numeric
  combined_cleaned_df[,2:3] <- apply(combined_cleaned_df[,2:3], 2, as.numeric)
}

# Export results for plotting
# Add column LocationName
combined_cleaned_df$LocationName <- "Passaic River"
write.csv(combined_cleaned_df,
          file = "Output/Data/Sites/csv/PassaicRiver/PassaicLmeObsPredPCB.csv",
          row.names = FALSE)

# Plot all the pairs together
p <- ggplot(combined_cleaned_df, aes(x = 10^(observed), y = 10^(predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(0.01, 10^4), 
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.01, 10^4), 
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration PCBi (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 15/15, 
        axis.title = element_text(size = 10)) +
  annotation_logticks(sides = "bl") +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  annotate("text", x = 0.1, y = 10^3.6,
           label = expression(atop("Passaic River",
                                   paste("6 PCB congeners (n = 77 pairs)"))),
           size = 4, fontface = 2)

# See plot
print(p)

# Save plot
ggsave("Output/Plots/Sites/ObsPred/PassaicRiver/PassaicLmeObsPredPCB.png",
       plot = p, width = 8, height = 8, dpi = 500)

         