## Water PCB concentrations data analysis per site
# Hudson River

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

# Select Hudson River data ---------------------------------------------------
hud.0 <- wdc[str_detect(wdc$SiteName, 'HudsonRiver'),]
# PCBs were discharged to the river from the General Electric
# (GE) manufacturing plants in Hudson Falls and Fort Edward, NY
# Dredging from 2009 to 2015
# https://www.epa.gov/system/files/documents/2021-08/hudson_summer2021_floodplainrifs_factsheet_final.pdf

# Data preparation --------------------------------------------------------
# Remove samples (rows) with total PCBs  = 0
hud.1 <- hud.0[!(rowSums(hud.0[, c(12:115)], na.rm = TRUE)==0),]
# Calculate total PCB
hud.tpcb <- rowSums(hud.1[, c(12:115)], na.rm = T)
# Change date format
hud.1$SampleDate <- as.Date(hud.1$SampleDate, format = "%m/%d/%y")
# Calculate sampling time
time.day <- data.frame(as.Date(hud.1$SampleDate) - min(as.Date(hud.1$SampleDate)))
# Create individual code for each site sampled
site.numb <- hud.1$SiteSampled %>% as.factor() %>% as.numeric
# Include season
yq.s <- as.yearqtr(as.yearmon(hud.1$SampleDate, "%m/%d/%Y") + 1/12)
season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                   labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
# Create data frame
hud.tpcb <- cbind(factor(hud.1$SiteSampled), hud.1$SampleDate,
                hud.1$Latitude, hud.1$Longitude, as.matrix(hud.tpcb),
                data.frame(time.day), site.numb, season.s)
# Add column manes
colnames(hud.tpcb) <- c("site", "date", "Latitude", "Longitude",
                        "tPCB", "time", "site.code", "season")

# Get coordinates per site
hud.location <- hud.tpcb[c('site', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
hud.location <- aggregate(tPCB ~ site + Latitude + Longitude,
                            data = hud.location, mean)

# (2) Calculate total log PCB
# Remove metadata
hud.log <- subset(hud.1, select = -c(ID:AroclorCongener))
# Remove Aroclor data
hud.log <- subset(hud.log, select = -c(A1016:A1260))
# Log 10 individual PCBs 
hud.log <- log10(hud.log)
# Replace -inf to NA
hud.log <- do.call(data.frame,
                   lapply(hud.log,
                          function(x) replace(x, is.infinite(x), NA)))
# Sum individual log 10 PCBs
hud.log.tpcb <- rowSums(hud.log, na.rm = T)
# Generate data.frame for analysis and plots
hud.log.tpcb <- cbind(factor(hud.1$SiteSampled), hud.1$SampleDate,
                      as.matrix(hud.log.tpcb), data.frame(time.day),
                      site.numb, season.s)
colnames(hud.log.tpcb) <- c("site", "date", "logtPCB", "time",
                            "site.code", "season")

# General plots -------------------------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(hud.tpcb$tPCB)
hist(log10(hud.tpcb$tPCB))
# (1.2) log.tPCB
hist(hud.log.tpcb$logtPCB)
hist(log10(hud.log.tpcb$logtPCB))

# (2) Time trend plots
# (2.1) tPCB
ggplot(hud.tpcb, aes(y = tPCB,
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
ggplot(hud.log.tpcb, aes(y = logtPCB,
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
ggplot(hud.tpcb, aes(x = season, y = tPCB)) +
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
ggplot(hud.log.tpcb, aes(x = season, y = logtPCB)) +
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
ggplot(hud.tpcb, aes(x = factor(site), y = tPCB)) + 
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
ggplot(hud.log.tpcb, aes(x = factor(site), y = logtPCB)) + 
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
# Remove site Bakers Falls. Upstream source
hud.tpcb.2 <- subset(hud.tpcb, site != c("BakersFalls"))
hud.log.tpcb.2 <- subset(hud.log.tpcb, site != c("BakersFalls"))

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Kalamazoo River
# Hudson River
sitehudN1 <- "01331095" # HUDSON RIVER AT STILLWATER NY
sitehudN2 <- "01359139" # HUDSON RIVER AT ALBANY NY No flow!
sitehudN3 <- "01359165" # HUDSON RIVER AT PORT OF ALBANY NY No flow!
sitehudN4 <- "01358000" # HUDSON RIVER AT GREEN ISLAND NY
sitehudN5 <- "01335755" # HUDSON RIVER AT LOCK 1 NEAR WATERFORD NY
sitehudN6 <- "01335754" # HUDSON RIVER ABOVE LOCK 1 NEAR WATERFORD NY
sitehudN7 <- "01328770" # HUDSON RIVER AT THOMSON NY
sitehudN8 <- "01327750" # HUDSON RIVER AT FORT EDWARD NY
sitehudN9 <- "01372043" # HUDSON RIVER NEAR POUGHKEEPSIE NY

# Codes to retrieve data
paramflow <- "00060" # discharge, ft3/s
paramtemp <- "00010" # water temperature, C Not available

# Flow (ft3/s)
flow.1 <- readNWISdv(sitehudN1, paramflow,
                   min(hud.tpcb$date), max(hud.tpcb$date))
flow.2 <- readNWISdv(sitehudN4, paramflow,
                     min(hud.tpcb$date), max(hud.tpcb$date))
flow.3 <- readNWISdv(sitehudN5, paramflow,
                     min(hud.tpcb$date), max(hud.tpcb$date))
flow.4 <- readNWISdv(sitehudN6, paramflow,
                     min(hud.tpcb$date), max(hud.tpcb$date))
flow.5 <- readNWISdv(sitehudN7, paramflow,
                     min(hud.tpcb$date), max(hud.tpcb$date))
flow.6 <- readNWISdv(sitehudN8, paramflow,
                     min(hud.tpcb$date), max(hud.tpcb$date))
flow.7 <- readNWISdv(sitehudN9, paramflow,
                     min(hud.tpcb$date), max(hud.tpcb$date))
temp <- readNWISdv(sitehudN2, paramtemp,
                   min(hud.tpcb$date), max(hud.tpcb$date))

# Add USGS data to fox.tpcb, matching dates
hud.tpcb$flow.1 <- flow.1$X_.Primary.Stream.Flow._00060_00003[match(hud.tpcb$date,
                                                                  flow.1$Date)]
hud.tpcb$flow.2 <- flow.2$X_.Primary.Stream.Flow._00060_00003[match(hud.tpcb$date,
                                                                     flow.2$Date)]
hud.tpcb$flow.3 <- flow.3$X_.Primary.Stream.Flow._00060_00003[match(hud.tpcb$date,
                                                                     flow.3$Date)]
hud.tpcb$flow.4 <- flow.4$X_.Primary.Stream.Flow._00060_00003[match(hud.tpcb$date,
                                                                     flow.4$Date)]
hud.tpcb$flow.5 <- flow.5$X_.Primary.Stream.Flow._00060_00003[match(hud.tpcb$date,
                                                                    flow.5$Date)]
hud.tpcb$flow.6 <- flow.6$X_.Primary.Stream.Flow._00060_00003[match(hud.tpcb$date,
                                                                    flow.6$Date)]
hud.tpcb$flow.7 <- flow.7$X_.Primary.Stream.Flow._00060_00003[match(hud.tpcb$date,
                                                                    flow.7$Date)]
hud.tpcb$temp <- temp$X_00010_00003[match(hud.tpcb$date, temp$Date)]
# Remove samples with temp = NA
hud.tpcb <- na.omit(hud.tpcb)

# Add USGS data to fox.log.tpcb, matching dates
hud.log.tpcb$flow.1 <- flow.1$X_.Primary.Stream.Flow._00060_00003[match(hud.log.tpcb$date,
                                                                      flow.1$Date)]
hud.log.tpcb$flow.2 <- flow.2$X_.Primary.Stream.Flow._00060_00003[match(hud.log.tpcb$date,
                                                                        flow.2$Date)]
hud.log.tpcb$flow.3 <- flow.3$X_.Primary.Stream.Flow._00060_00003[match(hud.log.tpcb$date,
                                                                        flow.3$Date)]
hud.log.tpcb$flow.4 <- flow.4$X_.Primary.Stream.Flow._00060_00003[match(hud.log.tpcb$date,
                                                                        flow.4$Date)]
hud.log.tpcb$flow.5 <- flow.5$X_.Primary.Stream.Flow._00060_00003[match(hud.log.tpcb$date,
                                                                        flow.5$Date)]
hud.log.tpcb$flow.6 <- flow.6$X_.Primary.Stream.Flow._00060_00003[match(hud.log.tpcb$date,
                                                                        flow.6$Date)]
hud.log.tpcb$flow.7 <- flow.7$X_.Primary.Stream.Flow._00060_00003[match(hud.log.tpcb$date,
                                                                        flow.7$Date)]
hud.log.tpcb$temp <- temp$X_00010_00003[match(hud.log.tpcb$date, temp$Date)]
# Remove samples with temp = NA
hud.log.tpcb <- na.omit(hud.log.tpcb)

# Remove samples with flow.1 = NA
kal.3 <- na.omit(kal.2)

# Regressions -------------------------------------------------------------
# (1) Perform linear regression (lr)
# (1.1) tPCB vs. time
lr.hud.tpcb.t <- lm(log10(tPCB) ~ time, data = hud.tpcb.2)
# See results
summary(lr.hud.tpcb.t)
# Look at residuals
res <- resid(lr.hud.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. time
lr.hud.log.tpcb.t <- lm(logtPCB ~ time, data = hud.log.tpcb.2)
# See results
summary(lr.hud.log.tpcb.t)
# Look at residuals
res <- resid(lr.hud.log.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB vs. season
lr.hud.tpcb.s <- lm(log10(tPCB) ~ season, data = hud.tpcb.2)
# See results
summary(lr.hud.tpcb.s)
# Look at residuals
res <- resid(lr.hud.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. season
lr.hud.log.tpcb.s <- lm(logtPCB ~ season, data = hud.log.tpcb.2)
# See results
summary(lr.hud.log.tpcb.s)
# Look at residuals
res <- resid(lr.hud.log.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.5) tPCB vs. flow (fox.tpcb.2)
lr.hud.tpcb.f <- lm(log10(tPCB) ~ flow, data = hud.tpcb.2)
# See results
summary(lr.hud.tpcb.f)
# Look at residuals
res <- resid(lr.hud.tpcb.f) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.6) log.tPCB vs. flow (fox.log.tpcb.2)
lr.hud.log.tpcb.f <- lm(logtPCB ~ flow, data = hud.log.tpcb.2)
# See results
summary(lr.hud.log.tpcb.f)
# Look at residuals
res <- resid(lr.hud.log.tpcb.f) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.7) tPCB vs. water temperature (fox.tpcb.2)
lr.hud.tpcb.te <- lm(log10(tPCB) ~ temp, data = hud.tpcb.2)
# See results
summary(lr.hud.tpcb.te)
# Look at residuals
res <- resid(lr.hud.tpcb.te) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.8) log.tPCB vs. temperature (hud.log.tpcb.2)
lr.hud.log.tpcb.te <- lm(logtPCB ~ temp, data = hud.log.tpcb.2)
# See results
summary(lr.hud.log.tpcb.te)
# Look at residuals
res <- resid(lr.hud.log.tpcb.te) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) MLR
# (2.1) tPCB vs. time + season + flow + temp (hud.tpcb.2)
mlr.hud.tpcb <- lm(log10(tPCB) ~ time + season + flow + temp,
                   data = hud.tpcb.2)
# See results
summary(mlr.hud.tpcb)
# Look at residuals
res <- resid(mlr.hud.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2.2) log.tPCB vs. time + season + flow + temp (hud.log.tpcb.2)
mlr.hud.log.tpcb <- lm(logtPCB ~ time + season + flow + temp,
                       data = hud.log.tpcb.2)
# See results
summary(mlr.hud.log.tpcb)
# Look at residuals
res <- resid(mlr.hud.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
