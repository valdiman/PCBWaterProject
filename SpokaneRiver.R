## Water PCB concentrations data analysis per site
# Spokane River

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

# Select spo River data ---------------------------------------------------
spo.0 <- wdc[str_detect(wdc$SiteName, 'SpokaneRiver'),]

# Data preparation --------------------------------------------------------
# Remove samples (rows) with total PCBs  = 0
spo.1 <- spo.0[!(rowSums(spo.0[, c(12:115)], na.rm = TRUE)==0),]
# Calculate total PCB
tpcb.spo <- rowSums(spo.1[, c(12:115)], na.rm = T)
# Change date format
spo.1$SampleDate <- as.Date(spo.1$SampleDate, format = "%m/%d/%y")
# Calculate sampling time
time.day <- data.frame(as.Date(spo.1$SampleDate) - min(as.Date(spo.1$SampleDate)))
# Create individual code for each site sampled
site.numb <- spo.1$SiteSampled %>% as.factor() %>% as.numeric
# Include season
yq.s <- as.yearqtr(as.yearmon(spo.1$SampleDate, "%m/%d/%Y") + 1/12)
season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                   labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
# Create data frame
spo.tpcb <- cbind(factor(spo.1$SiteSampled), spo.1$SampleDate,
                  spo.1$Latitude, spo.1$Longitude, as.matrix(tpcb.spo),
                  data.frame(time.day), site.numb, season.s)
# Add column names
colnames(spo.tpcb) <- c("site", "date", "Latitude", "Longitude",
                        "tPCB", "time", "site.code", "season")

# Get coordinates per site to plot in Google Earth
spo.location <- spo.tpcb[c('site', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
spo.location <- aggregate(tPCB ~ site + Latitude + Longitude,
                            data = spo.location, mean)

# (2) Calculate total log PCB
# Remove metadata
spo.log <- subset(spo.1, select = -c(ID:AroclorCongener))
# Remove Aroclor data
spo.log <- subset(spo.log, select = -c(A1016:A1260))
# Log 10 individual PCBs 
spo.log <- log10(spo.log)
# Replace -inf to NA
spo.log <- do.call(data.frame,
                     lapply(spo.log,
                            function(x) replace(x, is.infinite(x), NA)))
# Sum individual log 10 PCBs
spo.log.tpcb <- rowSums(spo.log, na.rm = T)
# Generate data.frame for analysis and plots
spo.log.tpcb <- cbind(factor(spo.1$SiteSampled), spo.1$SampleDate,
                      as.matrix(spo.log.tpcb), data.frame(time.day),
                      site.numb, season.s)
colnames(spo.log.tpcb) <- c("site", "date", "logtPCB", "time",
                            "site.code", "season")

# General plots -------------------------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(spo.tpcb$tPCB)
hist(log10(spo.tpcb$tPCB))
# (1.2) log.tPCB
hist(spo.log.tpcb$logtPCB)
hist(log10(spo.log.tpcb$logtPCB))

# (2) Time trend plots
# (2.1) tPCB
ggplot(spo.tpcb, aes(y = tPCB,
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
ggplot(spo.log.tpcb, aes(y = logtPCB,
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
ggplot(spo.tpcb, aes(x = season, y = tPCB)) +
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
ggplot(spo.log.tpcb, aes(x = season, y = logtPCB)) +
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
ggplot(spo.tpcb, aes(x = factor(site), y = tPCB)) + 
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
ggplot(spo.log.tpcb, aes(x = factor(site), y = logtPCB)) + 
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
# Remove site Lake Winnebago (background site)
# ?

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Spokane River
siteSpoN1 <- "12417650" # SPOKANE RIVER BLW BLACKWELL NR COEUR D ALENE ID
siteSpoN2 <- "12419000" # Spokane River near Post Falls, ID
siteSpoN3 <- "12422500" # Spokane River at Spokane, WA
siteSpoN4 <- "12424000" # Hangman Creek at Spokane, WA
siteSpoN5 <- "12426000" # SPOKANE RIVER BELOW NINE MILE DAM AT SPOKANE, WA
# Codes to retrieve data
paramflow <- "00060" # discharge, ft3/s
# paramtemp <- "00010" # water temperature, C
# Retrieve USGS data
flow.1 <- readNWISdv(siteSpoN1, paramflow,
                   min(spo.tpcb$date), max(spo.tpcb$date))
flow.2 <- readNWISdv(siteSpoN2, paramflow,
                     min(spo.tpcb$date), max(spo.tpcb$date))
flow.3 <- readNWISdv(siteSpoN3, paramflow,
                     min(spo.tpcb$date), max(spo.tpcb$date))
flow.4 <- readNWISdv(siteSpoN4, paramflow,
                     min(spo.tpcb$date), max(spo.tpcb$date))
flow.5 <- readNWISdv(siteSpoN5, paramflow,
                     min(spo.tpcb$date), max(spo.tpcb$date))

# Add USGS data to spo.tpcb, matching dates
spo.tpcb$flow.1 <- flow.1$X_00060_00003[match(spo.tpcb$date,
                                              flow$Date)]
spo.tpcb$flow.2 <- flow.2$X_00060_00003[match(spo.tpcb$date,
                                              flow$Date)]
spo.tpcb$flow.3 <- flow.3$X_00060_00003[match(spo.tpcb$date,
                                              flow$Date)]
spo.tpcb$flow.4 <- flow.4$X_00060_00003[match(spo.tpcb$date,
                                              flow$Date)]
spo.tpcb$flow.5 <- flow.5$X_00060_00003[match(spo.tpcb$date,
                                              flow$Date)]

# Remove samples with temp = NA
spo.tpcb.2 <- na.omit(spo.tpcb.2)

# Add USGS data to spo.log.tpcb, matching dates
spo.log.tpcb.2$flow <- flow$X_.Primary.Stream.Flow._00060_00003[match(spo.log.tpcb.2$date,
                                                                    flow$Date)]
spo.log.tpcb.2$temp <- temp$X_00010_00003[match(spo.log.tpcb.2$date, temp$Date)]
# Remove samples with temp = NA
spo.log.tpcb.2 <- na.omit(spo.log.tpcb.2)

# Regressions -------------------------------------------------------------
# (1) Perform linear regression (lr)
# (1.1) tPCB vs. time
lr.spo.tpcb.t <- lm(log10(tPCB) ~ time, data = spo.tpcb)
# See results
summary(lr.spo.tpcb.t)
# Look at residuals
res <- resid(lr.spo.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. time
lr.spo.log.tpcb.t <- lm(logtPCB ~ time, data = spo.log.tpcb.2)
# See results
summary(lr.spo.log.tpcb.t)
# Look at residuals
res <- resid(lr.spo.log.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB vs. season
lr.spo.tpcb.s <- lm(log10(tPCB) ~ season, data = spo.tpcb)
# See results
summary(lr.spo.tpcb.s)
# Look at residuals
res <- resid(lr.spo.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. season
lr.spo.log.tpcb.s <- lm(logtPCB ~ season, data = spo.log.tpcb.2)
# See results
summary(lr.spo.log.tpcb.s)
# Look at residuals
res <- resid(lr.spo.log.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.5) tPCB vs. flow (spo.tpcb.2)
lr.spo.tpcb.f <- lm(log10(tPCB) ~ flow.5, data = spo.tpcb)
# See results
summary(lr.spo.tpcb.f)
# Look at residuals
res <- resid(lr.spo.tpcb.f) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.6) log.tPCB vs. flow (spo.log.tpcb.2)
lr.spo.log.tpcb.f <- lm(logtPCB ~ flow, data = spo.log.tpcb.2)
# See results
summary(lr.spo.log.tpcb.f)
# Look at residuals
res <- resid(lr.spo.log.tpcb.f) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.7) tPCB vs. water temperature (spo.tpcb.2)
lr.spo.tpcb.te <- lm(log10(tPCB) ~ temp, data = spo.tpcb.2)
# See results
summary(lr.spo.tpcb.te)
# Look at residuals
res <- resid(lr.spo.tpcb.te) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.8) log.tPCB vs. temperature (spo.log.tpcb.2)
lr.spo.log.tpcb.te <- lm(logtPCB ~ temp, data = spo.log.tpcb.2)
# See results
summary(lr.spo.log.tpcb.te)
# Look at residuals
res <- resid(lr.spo.log.tpcb.te) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) MLR
# (2.1) tPCB vs. time + season + flow (spo.tpcb)
mlr.spo.tpcb <- lm(log10(tPCB) ~ time + season + flow.1,
                   data = spo.tpcb)
# See results
summary(mlr.spo.tpcb)
# Look at residuals
res <- resid(mlr.spo.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2.2) log.tPCB vs. time + season + flow + temp (spo.log.tpcb.2)
mlr.spo.log.tpcb <- lm(logtPCB ~ time + season + flow + temp,
                       data = spo.log.tpcb.2)
# See results
summary(mlr.spo.log.tpcb)
# Look at residuals
res <- resid(mlr.spo.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (3) Perform Linear Mixed-Effects Model (LMEM)
# (3.1) tPCB vs. time + season + flow + site (spo.tpcb.2)
tpcb <- spo.tpcb$tPCB
time <- spo.tpcb$time
site <- spo.tpcb$site.code
season <- spo.tpcb$season
flow <- spo.tpcb$flow.4

lmem.spo.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + season + flow + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.spo.tpcb)
# Look at residuals
res.spo.tpcb <- resid(lmem.spo.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.spo.tpcb, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.spo.tpcb)
# Shapiro test
shapiro.test(res.spo.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.spo.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.spo.tpcb))[1, 'R2c']

# Extract coefficient values
time.coeff <- summary(lmem.spo.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lmem.spo.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -log(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# (3.2) log.tPCB vs. time + season + flow + temp + site (spo.log.tpcb.2)
log.tpcb <- spo.log.tpcb.2$logtPCB
time <- spo.log.tpcb.2$time
site <- spo.log.tpcb.2$site.code
season <- spo.log.tpcb.2$season
flow <- spo.log.tpcb.2$flow
tem <- spo.log.tpcb.2$temp

lmem.spo.log.tpcb <- lmer(log.tpcb ~ 1 + time + season + season + flow + tem + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.spo.log.tpcb)
# Look at residuals
res.spo.log.tpcb <- resid(lmem.spo.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.spo.log.tpcb, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.spo.log.tpcb)
# Shapiro test
shapiro.test(res.spo.log.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res.spo.log.tpcb, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.spo.log.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.spo.log.tpcb))[1, 'R2c']

# Modeling plots
# (1) Get predicted values tpcb
fit.values.spo.tpcb <- as.data.frame(fitted(lmem.spo.tpcb))
# Add column name
colnames(fit.values.spo.tpcb) <- c("predicted")
# Add predicted values to data.frame
spo.tpcb$predicted <- 10^(fit.values.spo.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(spo.tpcb, aes(x = tPCB, y = predicted)) +
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

ggplot(spo.tpcb, aes(x = tPCB, y = predicted)) +
  geom_point() +
  scale_x_log10(limits = c(10, 1e4)) +
  scale_y_log10(limits = c(10, 1e4)) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 1.3) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 40, y = 10000,
           label = 'Spokane River', colour = 'black', size = 4,
           fontface = 2)

# Plot residuals vs. predictions
plot(log10(spo.tpcb.2$predicted), res.spo.tpcb)
abline(0, 0)

# (2) Get predicted values log.tpcb
fit.values.spo.log.tpcb <- as.data.frame(fitted(lmem.spo.log.tpcb))
# Add column name
colnames(fit.values.spo.log.tpcb) <- c("predicted")
# Add predicted values to data.frame
spo.log.tpcb.2$predicted <- fit.values.spo.log.tpcb$predicted

# Plot prediction vs. observations, 1:1 line
ggplot(spo.log.tpcb.2, aes(x = logtPCB, y = predicted)) +
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

ggplot(spo.log.tpcb.2, aes(x = logtPCB, y = predicted)) +
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
plot(spo.log.tpcb.2$predicted, res.spo.log.tpcb)
abline(0, 0)
