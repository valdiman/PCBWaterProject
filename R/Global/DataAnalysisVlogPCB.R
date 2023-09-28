## Water PCB concentrations analysis.
## Data were obtained from EPA and contractors from PCB Superfund
## sites in USA. Using log10 of individual PCBs and then sum them
## to get total PCB.

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
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor08052022.csv")

# Total Concentration Analysis --------------------------------------------
# Data preparation
{
  # Remove samples (rows) with total PCBs  = 0
  wdc.2 <- wdc[!(rowSums(wdc[, c(14:117)], na.rm = TRUE)==0),]
  # Calculate total PCB
  tpcb <- rowSums(wdc.2[, c(14:117)], na.rm = T)
  # Calculate total log PCB
  # Remove metadata
  wdc.3 <- subset(wdc.2, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  pcbi <- subset(wdc.3, select = -c(A1016:A1260))
  # Log10 individual PCBs 
  pcbi.log <- log10(pcbi)
  # Replace -inf to NA
  pcbi.log <- do.call(data.frame,
                      lapply(pcbi.log,
                             function(x) replace(x, is.infinite(x), NA)))
  # Sum individual log 10 PCBs
  tpcb.log <- rowSums(pcbi.log, na.rm = T)
  # Change date format
  wdc.2$SampleDate <- as.Date(wdc.2$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(wdc.2$SampleDate) - min(as.Date(wdc.2$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- wdc.2$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(wdc.2$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  tpcb.log <- cbind(factor(wdc.2$SiteID), wdc.2$SampleDate,
                wdc.2$Latitude, wdc.2$Longitude,
                as.matrix(tpcb.log), data.frame(time.day),
                site.numb, season.s)
  # Add column names
  colnames(tpcb.log) <- c("SiteID", "date", "Latitude", "Longitude",
                      "logtPCB", "time", "site.code",
                      "season")
}

# Global plots ------------------------------------------------------------
# Histogram
hist(tpcb.log$logtPCB)
hist(log10(tpcb.log$logtPCB))

## Total PCBs in 1 box plot
## include 64 pg/L from EPA
ggplot(tpcb.log, aes(x = "", y = logtPCB)) + 
  theme_classic() +
  theme(aspect.ratio = 14/2) +
  xlab(expression(bold(Sigma*"PCB (n = 5250)")))+
  ylab(expression(bold("Water Concentration 1990 - 2020 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 8),
        axis.title.x = element_text(face = "bold", size = 8, vjust = 5)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(lwd = 1.2, width = 0.7, outlier.shape = NA, alpha = 0)
  
# Regression analysis and plots---------------------------------------------
# Plots
# (1) Time trend plots
ggplot(tpcb.log, aes(y = logtPCB,
                 x = format(date,'%Y'))) +
  geom_point(shape = 21, cex = 1.2, fill = "#66ccff") +
  xlab("") +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB 1990 - 2020 (pg/L)"))))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black"))

# (2) Seasonality
ggplot(tpcb.log, aes(x = season, y = logtPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB 1990 - 2020 (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# Regressions -------------------------------------------------------------
# Get variables
log.tPCB <- tpcb.log$logtPCB
time <- tpcb.log$time
site <- tpcb.log$site.code
season <- tpcb.log$season
# (1) Perform linear regression (lr)
# log.tPCB vs. time
lr.log.tpcb.t <- lm(log.tPCB ~ time)
# See results
summary(lr.log.tpcb.t)
# Look at residuals
{
  res <- resid(lr.log.tpcb.t) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) log.tPCB vs. season
lr.log.tpcb.s <- lm(log.tPCB ~ season)
# See results
summary(lr.log.tpcb.s)
# Look at residuals
{
  res <- resid(lr.log.tpcb.s) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) MLR
# log.tPCB vs. time + season
mlr.log.tpcb <- lm(log.tPCB ~ time + season)
# See results
summary(mlr.log.tpcb)
# Look at residuals
{
  res <- resid(mlr.log.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (3) Perform Linear Mixed-Effects Model (lme)
# log.tPCB vs. time + season + site (wdc.log.tpcb)
lmem.log.tpcb <- lmer(log.tPCB ~ 1 + time + season + season + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.log.tpcb)
# Look at residuals
res.log.tpcb <- resid(lmem.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
{
  qqnorm(res.log.tpcb, main = "log10(C)")
  # Add a straight diagonal line to the plot
  qqline(res.log.tpcb)
}
# One-sample Kolmogorov-Smirnov test
ks.test(res.log.tpcb, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.log.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.log.tpcb))[1, 'R2c']

# (2) Get predicted values log.tpcb
fit.values.log.tpcb <- as.data.frame(fitted(lmem.log.tpcb))
# Add column name
colnames(fit.values.log.tpcb) <- c("predictedlog")
# Add predicted values to data.frame
tpcb.log$lmepredictedlog <- (fit.values.log.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(tpcb.log, aes(x = logtPCB, y = lmepredictedlog)) +
  geom_point(shape = 21, size = 2, fill = "#66ccff") +
  xlim(-200, 300) +
  ylim(-200, 300) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", size = 1) +
  theme_bw() +
  theme(aspect.ratio = 15/15)

