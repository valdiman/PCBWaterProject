## Water PCB concentrations data analysis per site
# Fox River

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

# Select Fox River data ---------------------------------------------------
fox.0 <- wdc[str_detect(wdc$SiteName, 'FoxRiver'),]
# Lake Winnebago is a background site.

# Data preparation --------------------------------------------------------
# Remove samples (rows) with total PCBs  = 0
fox.1 <- fox.0[!(rowSums(fox.0[, c(12:115)], na.rm = TRUE)==0),]
# Calculate total PCB
tpcb.fox <- rowSums(fox.1[, c(12:115)], na.rm = T)
# Change date format
fox.1$SampleDate <- as.Date(fox.1$SampleDate, format = "%m/%d/%y")
# Calculate sampling time
time.day <- data.frame(as.Date(fox.1$SampleDate) - min(as.Date(fox.1$SampleDate)))
# Create individual code for each site sampled
site.numb <- fox.1$SiteSampled %>% as.factor() %>% as.numeric
# Include season
yq.s <- as.yearqtr(as.yearmon(fox.1$SampleDate, "%m/%d/%Y") + 1/12)
season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                   labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
# Create data frame
fox.tpcb <- cbind(factor(fox.1$SiteSampled), fox.1$SampleDate,
                  fox.1$Latitude, fox.1$Longitude, as.matrix(tpcb.fox),
                  data.frame(time.day), site.numb, season.s)
# Add column names
colnames(fox.tpcb) <- c("site", "date", "Latitude", "Longitude",
                        "tPCB", "time", "site.code", "season")

# Get coordinates per site to plot in Google Earth
fox.location <- fox.tpcb[c('site', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
fox.location <- aggregate(tPCB ~ site + Latitude + Longitude,
                            data = fox.location, mean)

# (2) Calculate total log PCB
# Remove metadata
fox.log <- subset(fox.1, select = -c(ID:AroclorCongener))
# Remove Aroclor data
fox.log <- subset(fox.log, select = -c(A1016:A1260))
# Log 10 individual PCBs 
fox.log <- log10(fox.log)
# Replace -inf to NA
fox.log <- do.call(data.frame,
                     lapply(fox.log,
                            function(x) replace(x, is.infinite(x), NA)))
# Sum individual log 10 PCBs
fox.log.tpcb <- rowSums(fox.log, na.rm = T)
# Generate data.frame for analysis and plots
fox.log.tpcb <- cbind(factor(fox.1$SiteSampled), fox.1$SampleDate,
                      as.matrix(fox.log.tpcb), data.frame(time.day),
                      site.numb, season.s)
colnames(fox.log.tpcb) <- c("site", "date", "logtPCB", "time",
                            "site.code", "season")

# General plots -------------------------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(fox.tpcb$tPCB)
hist(log10(fox.tpcb$tPCB))
# (1.2) log.tPCB
hist(fox.log.tpcb$logtPCB)
hist(log10(fox.log.tpcb$logtPCB))

# (2) Time trend plots
# (2.1) tPCB
ggplot(fox.tpcb, aes(y = tPCB,
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
ggplot(fox.log.tpcb, aes(y = logtPCB,
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
ggplot(fox.tpcb, aes(x = season, y = tPCB)) +
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
ggplot(fox.log.tpcb, aes(x = season, y = logtPCB)) +
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
ggplot(fox.tpcb, aes(x = factor(site), y = tPCB)) + 
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
ggplot(fox.log.tpcb, aes(x = factor(site), y = logtPCB)) + 
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
fox.tpcb.2 <- subset(fox.tpcb, site != c("LakeWinnebago"))
fox.log.tpcb.2 <- subset(fox.log.tpcb, site != c("LakeWinnebago"))

# Plots w/o Lake Winnebago ------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(fox.tpcb.2$tPCB)
hist(log10(fox.tpcb.2$tPCB))
# (1.2) log.tPCB
hist(fox.log.tpcb.2$logtPCB)
hist(log10(fox.log.tpcb.2$logtPCB))

# (2) Time trend plots
# (2.1) tPCB
ggplot(fox.tpcb.2, aes(y = tPCB,
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
ggplot(fox.log.tpcb.2, aes(y = logtPCB,
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
ggplot(fox.tpcb.2, aes(x = season, y = tPCB)) +
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
ggplot(fox.log.tpcb.2, aes(x = season, y = logtPCB)) +
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
ggplot(fox.tpcb.2, aes(x = factor(site), y = tPCB)) + 
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
ggplot(fox.log.tpcb.2, aes(x = factor(site), y = logtPCB)) + 
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
# Include flow data from USGS station Fox River
sitefoxN1 <- "04084445" # flow @ OX RIVER AT APPLETON, WI
sitefoxN2 <- "040851385" # water temperature @ FOX RIVER AT OIL TANK DEPOT AT GREEN BAY, WI
# Codes to retrieve data
paramflow <- "00060" # discharge, ft3/s
paramtemp <- "00010" # water temperature, C
# Retrieve USGS data
flow <- readNWISdv(sitefoxN1, paramflow,
                   min(fox.tpcb.2$date), max(fox.tpcb.2$date))
temp <- readNWISdv(sitefoxN2, paramtemp,
                   min(fox.tpcb.2$date), max(fox.tpcb.2$date))
# Add USGS data to fox.tpcb, matching dates
fox.tpcb.2$flow <- flow$X_.Primary.Stream.Flow._00060_00003[match(fox.tpcb.2$date,
                                                                flow$Date)]
fox.tpcb.2$temp <- temp$X_00010_00003[match(fox.tpcb.2$date, temp$Date)]
# Remove samples with temp = NA
fox.tpcb.2 <- na.omit(fox.tpcb.2)

# Add USGS data to fox.log.tpcb, matching dates
fox.log.tpcb.2$flow <- flow$X_.Primary.Stream.Flow._00060_00003[match(fox.log.tpcb.2$date,
                                                                    flow$Date)]
fox.log.tpcb.2$temp <- temp$X_00010_00003[match(fox.log.tpcb.2$date, temp$Date)]
# Remove samples with temp = NA
fox.log.tpcb.2 <- na.omit(fox.log.tpcb.2)

# Regressions -------------------------------------------------------------
# (1) Perform linear regression (lr)
# (1.1) tPCB vs. time
lr.fox.tpcb.t <- lm(log10(tPCB) ~ time, data = fox.tpcb.2)
# See results
summary(lr.fox.tpcb.t)
# Look at residuals
res <- resid(lr.fox.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. time
lr.fox.log.tpcb.t <- lm(logtPCB ~ time, data = fox.log.tpcb.2)
# See results
summary(lr.fox.log.tpcb.t)
# Look at residuals
res <- resid(lr.fox.log.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB vs. season
lr.fox.tpcb.s <- lm(log10(tPCB) ~ season, data = fox.tpcb.2)
# See results
summary(lr.fox.tpcb.s)
# Look at residuals
res <- resid(lr.fox.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. season
lr.fox.log.tpcb.s <- lm(logtPCB ~ season, data = fox.log.tpcb.2)
# See results
summary(lr.fox.log.tpcb.s)
# Look at residuals
res <- resid(lr.fox.log.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.5) tPCB vs. flow (fox.tpcb.2)
lr.fox.tpcb.f <- lm(log10(tPCB) ~ flow, data = fox.tpcb.2)
# See results
summary(lr.fox.tpcb.f)
# Look at residuals
res <- resid(lr.fox.tpcb.f) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.6) log.tPCB vs. flow (fox.log.tpcb.2)
lr.fox.log.tpcb.f <- lm(logtPCB ~ flow, data = fox.log.tpcb.2)
# See results
summary(lr.fox.log.tpcb.f)
# Look at residuals
res <- resid(lr.fox.log.tpcb.f) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.7) tPCB vs. water temperature (fox.tpcb.2)
lr.fox.tpcb.te <- lm(log10(tPCB) ~ temp, data = fox.tpcb.2)
# See results
summary(lr.fox.tpcb.te)
# Look at residuals
res <- resid(lr.fox.tpcb.te) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.8) log.tPCB vs. temperature (fox.log.tpcb.2)
lr.fox.log.tpcb.te <- lm(logtPCB ~ temp, data = fox.log.tpcb.2)
# See results
summary(lr.fox.log.tpcb.te)
# Look at residuals
res <- resid(lr.fox.log.tpcb.te) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) MLR
# (2.1) tPCB vs. time + season + flow + temp (fox.tpcb.2)
mlr.fox.tpcb <- lm(log10(tPCB) ~ time + season + flow + temp, data = fox.tpcb.2)
# See results
summary(mlr.fox.tpcb)
# Look at residuals
res <- resid(mlr.fox.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2.2) log.tPCB vs. time + season + flow + temp (fox.log.tpcb.2)
mlr.fox.log.tpcb <- lm(logtPCB ~ time + season + flow + temp,
                       data = fox.log.tpcb.2)
# See results
summary(mlr.fox.log.tpcb)
# Look at residuals
res <- resid(mlr.fox.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (3) Perform Linear Mixed-Effects Model (LMEM)
# (3.1) tPCB vs. time + season + flow + temp + site (fox.tpcb.2)
tpcb <- fox.tpcb.2$tPCB
time <- fox.tpcb.2$time
site <- fox.tpcb.2$site.code
season <- fox.tpcb.2$season
flow <- fox.tpcb.2$flow
tem <- fox.tpcb.2$temp

lmem.fox.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + season + flow + tem + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.fox.tpcb)
# Look at residuals
res.fox.tpcb <- resid(lmem.fox.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.fox.tpcb, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.fox.tpcb)
# Shapiro test
shapiro.test(res.fox.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.fox.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.fox.tpcb))[1, 'R2c']

# (3.2) log.tPCB vs. time + season + flow + temp + site (fox.log.tpcb.2)
log.tpcb <- fox.log.tpcb.2$logtPCB
time <- fox.log.tpcb.2$time
site <- fox.log.tpcb.2$site.code
season <- fox.log.tpcb.2$season
flow <- fox.log.tpcb.2$flow
tem <- fox.log.tpcb.2$temp

lmem.fox.log.tpcb <- lmer(log.tpcb ~ 1 + time + season + season + flow + tem + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.fox.log.tpcb)
# Look at residuals
res.fox.log.tpcb <- resid(lmem.fox.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.fox.log.tpcb, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.fox.log.tpcb)
# Shapiro test
shapiro.test(res.fox.log.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res.fox.log.tpcb, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.fox.log.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.fox.log.tpcb))[1, 'R2c']

# Modeling plots
# (1) Get predicted values tpcb
fit.values.fox.tpcb <- as.data.frame(fitted(lmem.fox.tpcb))
# Add column name
colnames(fit.values.fox.tpcb) <- c("predicted")
# Add predicted values to data.frame
fox.tpcb.2$predicted <- 10^(fit.values.fox.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(fox.tpcb.2, aes(x = tPCB, y = predicted)) +
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

ggplot(fox.tpcb.2, aes(x = tPCB, y = predicted)) +
  geom_point() +
  scale_x_log10(limits = c(10, 1e4)) +
  scale_y_log10(limits = c(10, 1e4)) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 1.3) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Plot residuals vs. predictions
plot(log10(fox.tpcb.2$predicted), res.fox.tpcb)
abline(0, 0)

# (2) Get predicted values log.tpcb
fit.values.fox.log.tpcb <- as.data.frame(fitted(lmem.fox.log.tpcb))
# Add column name
colnames(fit.values.fox.log.tpcb) <- c("predicted")
# Add predicted values to data.frame
fox.log.tpcb.2$predicted <- fit.values.fox.log.tpcb$predicted

# Plot prediction vs. observations, 1:1 line
ggplot(fox.log.tpcb.2, aes(x = logtPCB, y = predicted)) +
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

ggplot(fox.log.tpcb.2, aes(x = logtPCB, y = predicted)) +
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
plot(fox.log.tpcb.2$predicted, res.fox.log.tpcb)
abline(0, 0)
