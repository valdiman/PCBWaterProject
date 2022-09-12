## Water PCB concentrations data analysis per site
# Portland Harbor

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

# Select Portland Harbor data ---------------------------------------------------
por.0 <- wdc[str_detect(wdc$SiteName, 'PortlandHarbor'),]

# Data preparation --------------------------------------------------------
# Remove samples (rows) with total PCBs  = 0
por.1 <- por.0[!(rowSums(por.0[, c(12:115)], na.rm = TRUE)==0),]
# Calculate total PCB
por.tpcb <- rowSums(por.1[, c(12:115)], na.rm = T)
# Change date format
por.1$SampleDate <- as.Date(por.1$SampleDate, format = "%m/%d/%y")
# Calculate sampling time
time.day <- data.frame(as.Date(por.1$SampleDate) - min(as.Date(por.1$SampleDate)))
# Create individual code for each site sampled
site.numb <- por.1$SiteSampled %>% as.factor() %>% as.numeric
# Include season
yq.s <- as.yearqtr(as.yearmon(por.1$SampleDate, "%m/%d/%Y") + 1/12)
season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                   labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
# Create data frame
por.tpcb <- cbind(factor(por.1$SiteSampled), por.1$SampleDate,
                por.1$Latitude, por.1$Longitude, as.matrix(por.tpcb),
                data.frame(time.day), site.numb, season.s)
# Add column manes
colnames(por.tpcb) <- c("site", "date", "Latitude", "Longitude",
                        "tPCB", "time", "site.code", "season")

# Get coordinates per site
por.location <- por.tpcb[c('site', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
por.location <- aggregate(tPCB ~ site + Latitude + Longitude,
                            data = por.location, mean)

# (2) Calculate total log PCB
# Remove metadata
por.log <- subset(por.1, select = -c(ID:AroclorCongener))
# Remove Aroclor data
por.log <- subset(por.log, select = -c(A1016:A1260))
# Log 10 individual PCBs 
por.log <- log10(por.log)
# Replace -inf to NA
por.log <- do.call(data.frame,
                   lapply(por.log,
                          function(x) replace(x, is.infinite(x), NA)))
# Sum individual log 10 PCBs
por.log.tpcb <- rowSums(por.log, na.rm = T)
# Generate data.frame for analysis and plots
por.log.tpcb <- cbind(factor(por.1$SiteSampled), por.1$SampleDate,
                      as.matrix(por.log.tpcb), data.frame(time.day),
                      site.numb, season.s)
colnames(por.log.tpcb) <- c("site", "date", "logtPCB", "time",
                            "site.code", "season")

# General plots -------------------------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(por.tpcb$tPCB)
hist(log10(por.tpcb$tPCB))
# (1.2) log.tPCB
hist(por.log.tpcb$logtPCB)
hist(log10(por.log.tpcb$logtPCB))

# (2) Time trend plots
# (2.1) tPCB
ggplot(por.tpcb, aes(y = tPCB,
                     x = format(date,'%Y%m'))) +
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
ggplot(por.log.tpcb, aes(y = logtPCB,
                         x = format(date,'%Y%m'))) +
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
ggplot(por.tpcb, aes(x = season, y = tPCB)) +
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
ggplot(por.log.tpcb, aes(x = season, y = logtPCB)) +
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
ggplot(por.tpcb, aes(x = factor(site), y = tPCB)) + 
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
ggplot(por.log.tpcb, aes(x = factor(site), y = logtPCB)) + 
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

# Remove site -------------------------------------------------------------
# ?

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Portland Harbor
sitePorN1 <- "14211720" # WILLAMETTE RIVER AT PORTLAND, OR
sitePorN2 <- "14211820" # COLUMBIA SLOUGH AT PORTLAND, OR No!

# Codes to retrieve data
paramflow <- "00060" # discharge, ft3/s
paramtemp <- "00010" # water temperature, C Not available

# Flow (ft3/s)
flow.1 <- readNWISdv(sitePorN1, paramflow,
                   min(por.tpcb$date), max(por.tpcb$date))
flow.2 <- readNWISdv(sitePorN2, paramflow,
                     min(por.tpcb$date), max(por.tpcb$date))
temp.1 <- readNWISdv(sitePorN1, paramtemp,
                   min(por.tpcb$date), max(por.tpcb$date))

# Add USGS data to fox.tpcb, matching dates
por.tpcb$flow.1 <- flow.1$X_00060_00003[match(por.tpcb$date, flow.1$Date)]
por.tpcb$flow.2 <- flow.2$X_00060_00003[match(por.tpcb$date, flow.2$Date)]
por.tpcb$temp.1 <- temp.1$X_00010_00003[match(por.tpcb$date, temp.1$Date)]
# Remove samples with temp = NA
por.tpcb.2 <- na.omit(por.tpcb)

# Add USGS data to fox.log.tpcb, matching dates
por.log.tpcb$flow.1 <- flow.1$X_00060_00003[match(por.log.tpcb$date,
                                                  flow.1$Date)]
por.log.tpcb$flow.2 <- flow.2$X_00060_00003[match(por.log.tpcb$date,
                                                  flow.2$Date)]
por.log.tpcb$temp.1 <- temp.1$X_00010_00003[match(por.log.tpcb$date,
                                                  temp.1$Date)]
# Remove samples with temp = NA
por.log.tpcb <- na.omit(por.log.tpcb)

# Regressions -------------------------------------------------------------
# (1) Perform linear regression (lr)
# (1.1) tPCB vs. time
lr.por.tpcb.t <- lm(log10(tPCB) ~ time, data = por.tpcb)
# See results
summary(lr.por.tpcb.t)
# Look at residuals
res <- resid(lr.por.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. time
lr.por.log.tpcb.t <- lm(logtPCB ~ time, data = por.log.tpcb)
# See results
summary(lr.por.log.tpcb.t)
# Look at residuals
res <- resid(lr.por.log.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB vs. season
lr.por.tpcb.s <- lm(log10(tPCB) ~ season, data = por.tpcb)
# See results
summary(lr.por.tpcb.s)
# Look at residuals
res <- resid(lr.por.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. season
lr.por.log.tpcb.s <- lm(logtPCB ~ season, data = por.log.tpcb)
# See results
summary(lr.por.log.tpcb.s)
# Look at residuals
res <- resid(lr.por.log.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.5) tPCB vs. flow (fox.tpcb.2)
lr.por.tpcb.f <- lm(log10(tPCB) ~ flow.1, data = por.tpcb)
# See results
summary(lr.por.tpcb.f)
# Look at residuals
res <- resid(lr.por.tpcb.f) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.6) log.tPCB vs. flow (fox.log.tpcb.2)
lr.por.log.tpcb.f <- lm(logtPCB ~ flow.1, data = por.log.tpcb)
# See results
summary(lr.por.log.tpcb.f)
# Look at residuals
res <- resid(lr.por.log.tpcb.f) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.7) tPCB vs. water temperature (fox.tpcb.2)
lr.por.tpcb.te <- lm(log10(tPCB) ~ temp.1, data = por.tpcb)
# See results
summary(lr.por.tpcb.te)
# Look at residuals
res <- resid(lr.por.tpcb.te) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.8) log.tPCB vs. temperature (por.log.tpcb.2)
lr.por.log.tpcb.te <- lm(logtPCB ~ temp.1, data = por.log.tpcb)
# See results
summary(lr.por.log.tpcb.te)
# Look at residuals
res <- resid(lr.por.log.tpcb.te) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) MLR
# (2.1) tPCB vs. time + season + flow + temp (por.tpcb)
mlr.por.tpcb <- lm(log10(tPCB) ~ time + season + flow.1 + temp.1,
                   data = por.tpcb)
# See results
summary(mlr.por.tpcb)
# Look at residuals
res <- resid(mlr.por.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2.2) log.tPCB vs. time + season + flow + temp (por.log.tpcb.2)
mlr.por.log.tpcb <- lm(logtPCB ~ time + season + flow.1 + temp.1,
                       data = por.log.tpcb)
# See results
summary(mlr.por.log.tpcb)
# Look at residuals
res <- resid(mlr.por.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (3) Perform Linear Mixed-Effects Model (LMEM)
# (3.1) tPCB vs. time + season + flow + temp + site (por.tpcb.2)
tpcb <- por.tpcb.2$tPCB
time <- por.tpcb.2$time
site <- por.tpcb.2$site.code
season <- por.tpcb.2$season
flow <- por.tpcb.2$flow.1
tem <- por.tpcb.2$temp.1

lmem.por.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + season + flow + tem + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.por.tpcb)
# Look at residuals
res.por.tpcb <- resid(lmem.por.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.por.tpcb, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.por.tpcb)
# Shapiro test
shapiro.test(res.por.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.por.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.por.tpcb))[1, 'R2c']

# Extract coefficient values
time.coeff <- summary(lmem.por.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lmem.por.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -log(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# Modeling plots
# (1) Get predicted values tpcb
fit.values.por.tpcb <- as.data.frame(fitted(lmem.por.tpcb))
# Add column name
colnames(fit.values.por.tpcb) <- c("predicted")
# Add predicted values to data.frame
por.tpcb.2$predicted <- 10^(fit.values.por.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(por.tpcb.2, aes(x = tPCB, y = predicted)) +
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

ggplot(por.tpcb.2, aes(x = tPCB, y = predicted)) +
  geom_point() +
  scale_x_log10(limits = c(10, 1e4)) +
  scale_y_log10(limits = c(10, 1e4)) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 1.3) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 50, y = 10000,
           label = 'Portland Harbor', colour = 'black', size = 4,
           fontface = 2)

# Plot residuals vs. predictions
plot(log10(por.tpcb.2$predicted), res.por.tpcb)
abline(0, 0)

