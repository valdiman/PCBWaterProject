## Water PCB concentrations data analysis per site
# Housatonic River

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

# Load libraries
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

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("WaterDataCongenerAroclor08052022.csv")

# Select Housatonic River data ---------------------------------------------------
hou.0 <- wdc[str_detect(wdc$SiteName, 'HousatonicRiver'),]

# Data preparation --------------------------------------------------------
# Remove samples (rows) with total PCBs  = 0
hou.1 <- hou.0[!(rowSums(hou.0[, c(12:115)], na.rm = TRUE)==0),]
# Calculate total PCB
tpcb.hou <- rowSums(hou.1[, c(12:115)], na.rm = T)
# Change date format
hou.1$SampleDate <- as.Date(hou.1$SampleDate, format = "%m/%d/%y")
# Calculate sampling time
time.day <- data.frame(as.Date(hou.1$SampleDate) - min(as.Date(hou.1$SampleDate)))
# Create individual code for each site sampled
site.numb <- hou.1$SiteSampled %>% as.factor() %>% as.numeric
# Include season
yq.s <- as.yearqtr(as.yearmon(hou.1$SampleDate, "%m/%d/%Y") + 1/12)
season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                   labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
# Create data frame
hou.tpcb <- cbind(factor(hou.1$SiteSampled), hou.1$SampleDate,
                  hou.1$Latitude, hou.1$Longitude, as.matrix(tpcb.hou),
                  data.frame(time.day), site.numb, season.s)
# Add column names
colnames(hou.tpcb) <- c("site", "date", "Latitude", "Longitude",
                        "tPCB", "time", "site.code", "season")

# Get coordinates per site to plot in Google Earth
hou.location <- hou.tpcb[c('site', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
hou.location <- aggregate(tPCB ~ site + Latitude + Longitude,
                            data = hou.location, mean)

# (2) Calculate total log PCB
# Remove metadata
hou.log <- subset(hou.1, select = -c(ID:AroclorCongener))
# Remove Aroclor data
hou.log <- subset(hou.log, select = -c(A1016:A1260))
# Log 10 individual PCBs 
hou.log <- log10(hou.log)
# Replace -inf to NA
hou.log <- do.call(data.frame,
                     lapply(hou.log,
                            function(x) replace(x, is.infinite(x), NA)))
# Sum individual log 10 PCBs
hou.log.tpcb <- rowSums(hou.log, na.rm = T)
# Generate data.frame for analysis and plots
hou.log.tpcb <- cbind(factor(hou.1$SiteSampled), hou.1$SampleDate,
                      as.matrix(hou.log.tpcb), data.frame(time.day),
                      site.numb, season.s)
colnames(hou.log.tpcb) <- c("site", "date", "logtPCB", "time",
                            "site.code", "season")

# General plots -------------------------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(hou.tpcb$tPCB)
hist(log10(hou.tpcb$tPCB))
# (1.2) log.tPCB
hist(hou.log.tpcb$logtPCB)
hist(log10(hou.log.tpcb$logtPCB))

# (2) Time trend plots
# (2.1) tPCB
ggplot(hou.tpcb, aes(y = tPCB,
                     x = format(date,'%Y'))) +
  geom_point() +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black"))

# (2.2) log.tPCB
ggplot(hou.log.tpcb, aes(y = logtPCB,
                         x = format(date,'%Y'))) +
  geom_point() +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black"))

# (3) Seasonality
# (3.1) tPCB
ggplot(hou.tpcb, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 2012 - 2018 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# (3.2) log.tPCB
ggplot(hou.log.tpcb, aes(x = season, y = logtPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 2012 - 2018 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# (4) Sites
# (4.1) tPCB
ggplot(hou.tpcb, aes(x = factor(site), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 2012 - 2018 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# (4.2) log.tPCB
ggplot(hou.log.tpcb, aes(x = factor(site), y = logtPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 2012 - 2018 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Housatonic River
siteHouRN1 <- "01197000" # EAST BRANCH HOUSATONIC RIVER AT COLTSVILLE, MA
siteHouRN2 <- "01197500" # HOUSATONIC RIVER NEAR GREAT BARRINGTON, MA
siteHouRN3 <- "01198125" # HOUSATONIC RIVER NEAR ASHLEY FALLS, MA
siteHouRN4 <- "01199000" # HOUSATONIC RIVER AT FALLS VILLAGE, CT
# Codes to retrieve data
paramflow <- "00060" # discharge, ft3/s
paramtemp <- "00010" # water temperature, C
# Retrieve USGS data
flow.1 <- readNWISdv(siteHouN1, paramflow,
                     min(hou.tpcb$date), max(hou.tpcb$date))
flow.2 <- readNWISdv(siteHouN2, paramflow,
                     min(hou.tpcb$date), max(hou.tpcb$date))
flow.3 <- readNWISdv(siteHouN3, paramflow,
                     min(hou.tpcb$date), max(hou.tpcb$date))
flow.4 <- readNWISdv(siteHouN4, paramflow,
                     min(hou.tpcb$date), max(hou.tpcb$date))

# Add USGS data to hou.tpcb, matching dates
hou.tpcb$flow.1 <- flow.1$X_00060_00003[match(hou.tpcb$date,
                                              flow.1$Date)]
hou.tpcb$flow.2 <- flow.2$X_00060_00003[match(hou.tpcb$date,
                                              flow.2$Date)]
hou.tpcb$flow.3 <- flow.3$X_00060_00003[match(hou.tpcb$date,
                                              flow.3$Date)]
hou.tpcb$flow.4 <- flow.4$X_00060_00003[match(hou.tpcb$date,
                                              flow.4$Date)]
hou.tpcb$flow.5 <- flow.5$X_00060_00003[match(hou.tpcb$date,
                                              flow$Date)]

# Add USGS data to hou.log.tpcb, matching dates
hou.log.tpcb$flow.1 <- flow.1$X_00060_00003[match(hou.log.tpcb$date,
                                                  flow.1$Date)]
hou.log.tpcb$flow.2 <- flow.2$X_00060_00003[match(hou.log.tpcb$date,
                                                  flow.2$Date)]
hou.log.tpcb$flow.3 <- flow.3$X_00060_00003[match(hou.log.tpcb$date,
                                                  flow.3$Date)]
hou.log.tpcb$flow.4 <- flow.4$X_00060_00003[match(hou.log.tpcb$date,
                                                  flow.4$Date)]

# Remove site -------------------------------------------------------------
# Silver Lake
hou.tpcb.2 <- subset(hou.tpcb, site != c("SilverLake"))


# Regressions -------------------------------------------------------------
# (1) Perform linear regression (lr)
# (1.1) tPCB vs. time
lr.hou.tpcb.t <- lm(log10(tPCB) ~ time, data = hou.tpcb)
# See results
summary(lr.hou.tpcb.t)
# Look at residuals
res <- resid(lr.hou.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. time
lr.hou.log.tpcb.t <- lm(logtPCB ~ time, data = hou.log.tpcb)
# See results
summary(lr.hou.log.tpcb.t)
# Look at residuals
res <- resid(lr.hou.log.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB vs. season
lr.hou.tpcb.s <- lm(log10(tPCB) ~ season, data = hou.tpcb)
# See results
summary(lr.hou.tpcb.s)
# Look at residuals
res <- resid(lr.hou.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. season
lr.hou.log.tpcb.s <- lm(logtPCB ~ season, data = hou.log.tpcb)
# See results
summary(lr.hou.log.tpcb.s)
# Look at residuals
res <- resid(lr.hou.log.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.5) tPCB vs. flow
lr.hou.tpcb.f <- lm(log10(tPCB) ~ flow.4, data = hou.tpcb)
# See results
summary(lr.hou.tpcb.f)
# Look at residuals
res <- resid(lr.hou.tpcb.f) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.6) log.tPCB vs. flow
lr.hou.log.tpcb.f <- lm(logtPCB ~ flow.4, data = hou.log.tpcb)
# See results
summary(lr.hou.log.tpcb.f)
# Look at residuals
res <- resid(lr.hou.log.tpcb.f) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) MLR
# (2.1) tPCB vs. time + season + flow
mlr.hou.tpcb <- lm(log10(tPCB) ~ time + season,
                   data = hou.tpcb)
# See results
summary(mlr.hou.tpcb)
# Look at residuals
res <- resid(mlr.hou.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2.2) log.tPCB vs. time + season + flow
mlr.hou.log.tpcb <- lm(logtPCB ~ time + season + flow.1,
                       data = hou.log.tpcb)
# See results
summary(mlr.hou.log.tpcb)
# Look at residuals
res <- resid(mlr.hou.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (3) Perform Linear Mixed-Effects Model (LMEM)
# (3.1) tPCB vs. time + season + flow + site
tpcb <- hou.tpcb$tPCB
time <- hou.tpcb$time
site <- hou.tpcb$site.code
season <- hou.tpcb$season
flow <- hou.tpcb$flow.1

lmem.hou.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + season + flow + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.hou.tpcb)
# Look at residuals
res.hou.tpcb <- resid(lmem.hou.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.hou.tpcb, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.hou.tpcb)
# Shapiro test
shapiro.test(res.hou.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.hou.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.hou.tpcb))[1, 'R2c']

# Extract coefficient values
time.coeff <- summary(lmem.hou.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lmem.hou.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -log(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# (3.2) log.tPCB vs. time + season + flow + temp + site (hou.log.tpcb.2)
log.tpcb <- hou.log.tpcb.2$logtPCB
time <- hou.log.tpcb.2$time
site <- hou.log.tpcb.2$site.code
season <- hou.log.tpcb.2$season
flow <- hou.log.tpcb.2$flow
tem <- hou.log.tpcb.2$temp

lmem.hou.log.tpcb <- lmer(log.tpcb ~ 1 + time + season + season + flow + tem + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.hou.log.tpcb)
# Look at residuals
res.hou.log.tpcb <- resid(lmem.hou.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.hou.log.tpcb, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.hou.log.tpcb)
# Shapiro test
shapiro.test(res.hou.log.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res.hou.log.tpcb, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.hou.log.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.hou.log.tpcb))[1, 'R2c']

# Modeling plots
# (1) Get predicted values tpcb
fit.values.hou.tpcb <- as.data.frame(fitted(lmem.hou.tpcb))
# Add column name
colnames(fit.values.hou.tpcb) <- c("predicted")
# Add predicted values to data.frame
hou.tpcb$predicted <- 10^(fit.values.hou.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(hou.tpcb, aes(x = tPCB, y = predicted)) +
  geom_point() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 9),
        axis.title.x = element_text(face = "bold", size = 9)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "bl") +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 1.3)

ggplot(hou.tpcb, aes(x = tPCB, y = predicted)) +
  geom_point() +
  scale_x_log10(limits = c(10, 1e5)) +
  scale_y_log10(limits = c(10, 1e5)) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 1.3) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 40, y = 10000,
           label = 'houkane River', colour = 'black', size = 4,
           fontface = 2)

# Plot residuals vs. predictions
plot(log10(hou.tpcb.2$predicted), res.hou.tpcb)
abline(0, 0)

# (2) Get predicted values log.tpcb
fit.values.hou.log.tpcb <- as.data.frame(fitted(lmem.hou.log.tpcb))
# Add column name
colnames(fit.values.hou.log.tpcb) <- c("predicted")
# Add predicted values to data.frame
hou.log.tpcb.2$predicted <- fit.values.hou.log.tpcb$predicted

# Plot prediction vs. observations, 1:1 line
ggplot(hou.log.tpcb.2, aes(x = logtPCB, y = predicted)) +
  geom_point() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 9),
        axis.title.x = element_text(face = "bold", size = 9)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "bl") +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 1.3)

ggplot(hou.log.tpcb.2, aes(x = logtPCB, y = predicted)) +
  geom_point() +
  scale_x_log10(limits = c(5, 1e2)) +
  scale_y_log10(limits = c(5, 1e2)) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 1.3) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Plot residuals vs. predictions
plot(hou.log.tpcb.2$predicted, res.hou.log.tpcb)
abline(0, 0)
