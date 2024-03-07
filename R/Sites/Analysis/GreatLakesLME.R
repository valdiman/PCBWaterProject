## Water PCB concentrations data analysis
## Lake Michigan Mass Balance & Great Lakes
## Only Linear Mixed-Effects Model (lme)

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

# Select LMMB and Great Lakes data ---------------------------------------------------
grl <- wdc[str_detect(wdc$LocationName, 'Lake Michigan Mass Balance|Great Lakes'), ]

# (1) Just get lake data, remove data from tributaries
grl <- grl[!grepl("^Tributary", grl$SiteName), ]

# Data preparation --------------------------------------------------------
{
  # Change date format
  grl$SampleDate <- as.Date(grl$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(grl$SampleDate) - min(as.Date(grl$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- grl$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(grl$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  grl.tpcb <- cbind(factor(grl$SiteID), grl$SampleDate,
                    grl$Latitude, grl$Longitude, as.matrix(grl$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(grl.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Add water temperature data ----------------------------------------------
# See code: R/ExtractingData/LMMB/WaterTemp.R
{
  # Read water temperature
  wtp <- read.csv("Output/Data/Sites/csv/GreatLakes/WaterTemp/LakeMichiganWT.csv")
  # Convert date columns to Date format
  wtp$Date <- as.Date(wtp$Date)
  # Add water temperature to grl.tpcb
  grl.tpcb$temp <- wtp$WTMP_K[match(grl.tpcb$date, wtp$Date)]
  # Remove samples with temp = NA
  grl.tpcb <- na.omit(grl.tpcb)
}

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- grl.tpcb$tPCB
time <- grl.tpcb$time
wtemp <- grl.tpcb$temp
site <- grl.tpcb$site.code
season <- grl.tpcb$season
lme.grl.tpcb <- lmer(log10(tpcb) ~ 1 + time + wtemp + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.grl.tpcb)

# Shapiro test
shapiro.test(resid(lme.grl.tpcb))  # Lme doesn't work, p-value < 0.05

# (2) Samples from lake Michigan
grl.tpcb.1 <- subset(grl.tpcb, grepl("LMM", SiteID))

# Using grl.tpcb.1
# Get variables
tpcb <- grl.tpcb.1$tPCB
time <- grl.tpcb.1$time
wtemp <- grl.tpcb.1$temp
site <- grl.tpcb.1$site.code
season <- grl.tpcb.1$season
lme.grl.tpcb <- lmer(log10(tpcb) ~ 1 + time + wtemp + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.grl.tpcb)

# Shapiro test
shapiro.test(resid(lme.grl.tpcb))  # Lme doesn't work, p-value < 0.05

# LME for individual PCBs -------------------------------------------------
# Prepare data.frame
{
  # Select only Lake Michigan data
  grl.pcb.0 <- subset(grl, grepl("LMM", SiteID))
  # Remove metadata
  grl.pcb <- subset(grl.pcb.0, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  grl.pcb <- subset(grl.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  grl.pcb <- log10(grl.pcb)
  # Replace -inf to NA
  grl.pcb <- do.call(data.frame,
                     lapply(grl.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  grl.pcb.1 <- grl.pcb[,
                       -which(colSums(is.na(grl.pcb))/nrow(grl.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(grl.pcb.0$SiteID)
  # Change date format
  SampleDate <- as.Date(grl.pcb.0$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Create individual code for each site sampled
  site.numb <- grl.pcb.0$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(grl.pcb.0$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to grl.pcb.1
  grl.pcb.1 <- cbind(grl.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     site.numb, season.s)
  # Add water temperature to grl.tpcb
  grl.pcb.1$temp <- wtp$WTMP_K[match(grl.pcb.1$SampleDate, wtp$Date)]
  # Remove samples with temp = NA
  grl.pcb.1 <- grl.pcb.1[!is.na(grl.pcb.1$temp), ]
  # Remove metadata
  grl.pcb.2 <- subset(grl.pcb.1, select = -c(SiteID:temp))
}

# Get covariates
time <- grl.pcb.1$time
wtemp <- grl.pcb.1$temp
season <- grl.pcb.1$season
site <- grl.pcb.1$site.numb

# Create matrix to store results
lme.pcb <- matrix(nrow = length(grl.pcb.2[1,]), ncol = 22)

# Perform LME
for (i in 1:length(grl.pcb.2[1,])) {
  fit <- lmer(grl.pcb.2[,i] ~ 1 + time + wtemp + season + (1|site),
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
  lme.pcb[i,7] <- fixef(fit)[2] # water temperature
  lme.pcb[i,8] <- summary(fit)$coef[2,"Std. Error"] # water temperature error
  lme.pcb[i,9] <- summary(fit)$coef[2,"Pr(>|t|)"] # water temperature p-value
  lme.pcb[i,10] <- fixef(fit)[3] # # season 2
  lme.pcb[i,11] <- summary(fit)$coef[3,"Std. Error"] # season 2 error
  lme.pcb[i,12] <- summary(fit)$coef[3,"Pr(>|t|)"] # # season 2 p-value
  lme.pcb[i,13] <- fixef(fit)[4] # season 3
  lme.pcb[i,14] <- summary(fit)$coef[4,"Std. Error"] # season 3 error
  lme.pcb[i,15] <- summary(fit)$coef[4,"Pr(>|t|)"] # season 3 p-value
  lme.pcb[i,16] <- -log(2)/lme.pcb[i,4]/365 # t0.5
  lme.pcb[i,17] <- abs(-log(2)/lme.pcb[i,4]/365)*lme.pcb[i,5]/abs(lme.pcb[i,4]) # t0.5 error
  lme.pcb[i,18] <- as.data.frame(VarCorr(fit))[1,'sdcor']
  lme.pcb[i,19] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
  lme.pcb[i,20] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
  lme.pcb[i,21] <- shapiro.test(resid(fit))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(fit)
  # Calculate residuals and RMSE
  residuals <- grl.pcb.2[, i] - predictions
  non_na_indices <- !is.na(residuals)
  lme.pcb[i, 22] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Transform result to data.frame so factor 2 can be included
lme.pcb <- as.data.frame(lme.pcb)

# Add factor of 2
# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(grl.pcb.2[,1]),
                      ncol = length(grl.pcb.2[1,]))

# Create a vector to store factor 2 for each congener
factor2_vector <- numeric(length = length(grl.pcb.2[1,]))

for (i in 1:length(grl.pcb.2[1,])) {
  fit <- lmer(grl.pcb.2[,i] ~ 1 + time + wtemp + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  
  lme.fit.pcb[,i] <- fitted(fit)
  
  # Calculate factor2 for each congener
  factor2 <- 10^(lme.fit.pcb[, i])/10^(grl.pcb.2[, i])
  factor2_vector[i] <- sum(factor2 > 0.5 & factor2 < 2,
                           na.rm = TRUE) / (sum(!is.na(factor2))) * 100
}

# Add factor 2 to lme.pcb data.frame
lme.pcb$factor2 <- factor2_vector

# Change number format of factor 2 to 3 significant figures
lme.pcb$factor2 <- formatC(signif(lme.pcb$factor2, digits = 3))

# Add congener names
congeners <- colnames(grl.pcb.2)
lme.pcb <- as.data.frame(cbind(congeners, lme.pcb))

# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "temperature", "temperature.error", "temperature.pv",
                       "season2", "season2.error", "season2, pv", "season3",
                       "season3.error", "season3.pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality",
                       "RMSE", "Factor2")

# Remove congeners with no normal distribution
# Shapiro test p-value < 0.05
lme.pcb$Normality <- as.numeric(lme.pcb$Normality)
# Get the congeners that are not showing normality
lme.pcb.out <- lme.pcb[lme.pcb$Normality < 0.05, ]
lme.pcb <- lme.pcb[lme.pcb$Normality > 0.05, ]

# Export results
write.csv(lme.pcb, file = "Output/Data/Sites/csv/GreatLakes/GreatLakesLmePCB.csv",
          row.names = FALSE)

# Obtain observations vs predictions
# Generate predictions
# Select congeners that are not showing normality to be remove from grl.pcb.2
df <- data.frame(names_to_remove = lme.pcb.out$Congeners)
# Get column indices to remove
cols_to_remove <- which(names(grl.pcb.2) %in% df$names_to_remove)
# Remove columns from fox.pcb.2 with congeners that don't show normality
grl.pcb.3 <- grl.pcb.2[, -cols_to_remove]

# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(grl.pcb.3[,1]),
                      ncol = length(grl.pcb.3[1,]))

for (i in 1:length(grl.pcb.3[1,])) {
  fit <- lmer(grl.pcb.3[,i] ~ 1 + time + wtemp + season + (1|site),
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
colnames(lme.fit.pcb) <- colnames(grl.pcb.3)
# Add code number to first column
df1 <- cbind(code = row.names(grl.pcb.3), grl.pcb.3)
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
    annotate('text', x = 0.5, y = 10^4, label = gsub("\\.", "+", names(df1)[i]),
             size = 3, fontface = 2)
  # save plot
  ggsave(paste0("Output/Plots/Sites/ObsPred/GreatLakes/", col_name, ".png"), plot = p,
         width = 6, height = 6, dpi = 500)
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
    annotate('text', x = 0.5, y = 10^4, label = gsub("\\.", "+", col_name),
             size = 2.5, fontface = 2) +
    geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
    geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
    geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7)
  
  plot_list[[i-1]] <- p  # add plot to list
}
# Combine all the plots using patchwork
combined_plot <- wrap_plots(plotlist = plot_list, ncol = 4)
# Save the combined plot
ggsave("Output/Plots/Sites/ObsPred/GreatLakes/LmeCombined_plot.png", plot = combined_plot,
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
combined_cleaned_df$LocationName <- "Great Lakes"
write.csv(combined_cleaned_df,
          file = "Output/Data/Sites/csv/GreatLakes/GreatLakesLmeObsPredPCB.csv",
          row.names = FALSE)

# Plot all the pairs together
p <- ggplot(combined_cleaned_df, aes(x = 10^(observed), y = 10^(predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(0.01, 10^3), 
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.01, 10^3), 
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
  annotate("text", x = 1, y = 10^3.7,
           label = expression(atop("Great Lakes",
                                   paste("5 PCB congeners (n = 381 pairs)"))),
           size = 4, fontface = 2)

# See plot
print(p)

# Save plot
ggsave("Output/Plots/Sites/ObsPred/GreatLakes/GreatLakesLmeObsPredPCB.png",
       plot = p, width = 8, height = 8, dpi = 500)


