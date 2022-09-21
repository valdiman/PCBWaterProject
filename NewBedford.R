## Water PCB concentrations data analysis per site
# New Bedford Harbor

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

# Select nbh River data ---------------------------------------------------
nbh.0 <- wdc[str_detect(wdc$SiteName, 'NewBedford'),]

# Data preparation --------------------------------------------------------
# Remove samples (rows) with total PCBs  = 0
nbh.1 <- nbh.0[!(rowSums(nbh.0[, c(12:115)], na.rm = TRUE)==0),]
# Calculate total PCB
tpcb.nbh <- rowSums(nbh.1[, c(12:115)], na.rm = T)
# Change date format
nbh.1$SampleDate <- as.Date(nbh.1$SampleDate, format = "%m/%d/%y")
# Calculate sampling time
time.day <- data.frame(as.Date(nbh.1$SampleDate) - min(as.Date(nbh.1$SampleDate)))
# Create individual code for each site sampled
site.numb <- nbh.1$SiteSampled %>% as.factor() %>% as.numeric
# Include season
yq.s <- as.yearqtr(as.yearmon(nbh.1$SampleDate, "%m/%d/%Y") + 1/12)
season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                   labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
# Create data frame
nbh.tpcb <- cbind(factor(nbh.1$SiteSampled), nbh.1$SampleDate,
                  nbh.1$Latitude, nbh.1$Longitude, as.matrix(tpcb.nbh),
                  data.frame(time.day), site.numb, season.s)
# Add column names
colnames(nbh.tpcb) <- c("site", "date", "Latitude", "Longitude",
                        "tPCB", "time", "site.code", "season")

# Get coordinates per site to plot in Google Earth
nbh.location <- nbh.tpcb[c('site', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
nbh.location <- aggregate(tPCB ~ site + Latitude + Longitude,
                            data = nbh.location, mean)

# (2) Calculate total log PCB
# Remove metadata
nbh.log <- subset(nbh.1, select = -c(ID:AroclorCongener))
# Remove Aroclor data
nbh.log <- subset(nbh.log, select = -c(A1016:A1260))
# Log 10 individual PCBs 
nbh.log <- log10(nbh.log)
# Replace -inf to NA
nbh.log <- do.call(data.frame,
                     lapply(nbh.log,
                            function(x) replace(x, is.infinite(x), NA)))
# Sum individual log 10 PCBs
nbh.log.tpcb <- rowSums(nbh.log, na.rm = T)
# Generate data.frame for analysis and plots
nbh.log.tpcb <- cbind(factor(nbh.1$SiteSampled), nbh.1$SampleDate,
                      as.matrix(nbh.log.tpcb), data.frame(time.day),
                      site.numb, season.s)
colnames(nbh.log.tpcb) <- c("site", "date", "logtPCB", "time",
                            "site.code", "season")

# General plots -------------------------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(nbh.tpcb$tPCB)
hist(log10(nbh.tpcb$tPCB))
# (1.2) log.tPCB
hist(nbh.log.tpcb$logtPCB)
hist(log10(nbh.log.tpcb$logtPCB))

# (2) Time trend plots
# (2.1) tPCB
ggplot(nbh.tpcb, aes(y = tPCB,
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
ggplot(nbh.log.tpcb, aes(y = logtPCB,
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
ggplot(nbh.tpcb, aes(x = season, y = tPCB)) +
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
ggplot(nbh.log.tpcb, aes(x = season, y = logtPCB)) +
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
ggplot(nbh.tpcb, aes(x = factor(site), y = tPCB)) + 
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
ggplot(nbh.log.tpcb, aes(x = factor(site), y = logtPCB)) + 
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
nbh.tpcb.2 <- subset(nbh.tpcb, site != c("LakeWinnebago"))
nbh.log.tpcb.2 <- subset(nbh.log.tpcb, site != c("LakeWinnebago"))

# Plots w/o Lake Winnebago ------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(nbh.tpcb.2$tPCB)
hist(log10(nbh.tpcb.2$tPCB))
# (1.2) log.tPCB
hist(nbh.log.tpcb.2$logtPCB)
hist(log10(nbh.log.tpcb.2$logtPCB))

# (2) Time trend plots
# (2.1) tPCB
ggplot(nbh.tpcb.2, aes(y = tPCB,
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
ggplot(nbh.log.tpcb.2, aes(y = logtPCB,
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
ggplot(nbh.tpcb.2, aes(x = season, y = tPCB)) +
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
ggplot(nbh.log.tpcb.2, aes(x = season, y = logtPCB)) +
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
ggplot(nbh.tpcb.2, aes(x = factor(site), y = tPCB)) + 
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
ggplot(nbh.log.tpcb.2, aes(x = factor(site), y = logtPCB)) + 
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
# Include flow data from USGS station nbh River
sitenbhN1 <- "04084445" # flow @ OX RIVER AT APPLETON, WI
sitenbhN2 <- "040851385" # water temperature @ nbh RIVER AT OIL TANK DEPOT AT GREEN BAY, WI
# Codes to retrieve data
paramflow <- "00060" # discharge, ft3/s
paramtemp <- "00010" # water temperature, C
# Retrieve USGS data
flow <- readNWISdv(sitenbhN1, paramflow,
                   min(nbh.tpcb.2$date), max(nbh.tpcb.2$date))
temp <- readNWISdv(sitenbhN2, paramtemp,
                   min(nbh.tpcb.2$date), max(nbh.tpcb.2$date))
# Add USGS data to nbh.tpcb, matching dates
nbh.tpcb.2$flow <- flow$X_.Primary.Stream.Flow._00060_00003[match(nbh.tpcb.2$date,
                                                                flow$Date)]
nbh.tpcb.2$temp <- temp$X_00010_00003[match(nbh.tpcb.2$date, temp$Date)]
# Remove samples with temp = NA
nbh.tpcb.2 <- na.omit(nbh.tpcb.2)

# Add USGS data to nbh.log.tpcb, matching dates
nbh.log.tpcb.2$flow <- flow$X_.Primary.Stream.Flow._00060_00003[match(nbh.log.tpcb.2$date,
                                                                    flow$Date)]
nbh.log.tpcb.2$temp <- temp$X_00010_00003[match(nbh.log.tpcb.2$date, temp$Date)]
# Remove samples with temp = NA
nbh.log.tpcb.2 <- na.omit(nbh.log.tpcb.2)

# Regressions -------------------------------------------------------------
# (1) Perform linear regression (lr)
# (1.1) tPCB vs. time
lr.nbh.tpcb.t <- lm(log10(tPCB) ~ time, data = nbh.tpcb.2)
# See results
summary(lr.nbh.tpcb.t)
# Look at residuals
res <- resid(lr.nbh.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. time
lr.nbh.log.tpcb.t <- lm(logtPCB ~ time, data = nbh.log.tpcb.2)
# See results
summary(lr.nbh.log.tpcb.t)
# Look at residuals
res <- resid(lr.nbh.log.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB vs. season
lr.nbh.tpcb.s <- lm(log10(tPCB) ~ season, data = nbh.tpcb.2)
# See results
summary(lr.nbh.tpcb.s)
# Look at residuals
res <- resid(lr.nbh.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. season
lr.nbh.log.tpcb.s <- lm(logtPCB ~ season, data = nbh.log.tpcb.2)
# See results
summary(lr.nbh.log.tpcb.s)
# Look at residuals
res <- resid(lr.nbh.log.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.5) tPCB vs. flow (nbh.tpcb.2)
lr.nbh.tpcb.f <- lm(log10(tPCB) ~ flow, data = nbh.tpcb.2)
# See results
summary(lr.nbh.tpcb.f)
# Look at residuals
res <- resid(lr.nbh.tpcb.f) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.6) log.tPCB vs. flow (nbh.log.tpcb.2)
lr.nbh.log.tpcb.f <- lm(logtPCB ~ flow, data = nbh.log.tpcb.2)
# See results
summary(lr.nbh.log.tpcb.f)
# Look at residuals
res <- resid(lr.nbh.log.tpcb.f) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.7) tPCB vs. water temperature (nbh.tpcb.2)
lr.nbh.tpcb.te <- lm(log10(tPCB) ~ temp, data = nbh.tpcb.2)
# See results
summary(lr.nbh.tpcb.te)
# Look at residuals
res <- resid(lr.nbh.tpcb.te) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.8) log.tPCB vs. temperature (nbh.log.tpcb.2)
lr.nbh.log.tpcb.te <- lm(logtPCB ~ temp, data = nbh.log.tpcb.2)
# See results
summary(lr.nbh.log.tpcb.te)
# Look at residuals
res <- resid(lr.nbh.log.tpcb.te) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) MLR
# (2.1) tPCB vs. time + season + flow + temp (nbh.tpcb.2)
mlr.nbh.tpcb <- lm(log10(tPCB) ~ time + season + flow + temp, data = nbh.tpcb.2)
# See results
summary(mlr.nbh.tpcb)
# Look at residuals
res <- resid(mlr.nbh.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2.2) log.tPCB vs. time + season + flow + temp (nbh.log.tpcb.2)
mlr.nbh.log.tpcb <- lm(logtPCB ~ time + season + flow + temp,
                       data = nbh.log.tpcb.2)
# See results
summary(mlr.nbh.log.tpcb)
# Look at residuals
res <- resid(mlr.nbh.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (3) Perform Linear Mixed-Effects Model (LMEM)
# (3.1) tPCB vs. time + season + flow + temp + site (nbh.tpcb.2)
tpcb <- nbh.tpcb.2$tPCB
time <- nbh.tpcb.2$time
site <- nbh.tpcb.2$site.code
season <- nbh.tpcb.2$season
flow <- nbh.tpcb.2$flow
tem <- nbh.tpcb.2$temp

lmem.nbh.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + season + flow + tem + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.nbh.tpcb)
# Look at residuals
res.nbh.tpcb <- resid(lmem.nbh.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.nbh.tpcb, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.nbh.tpcb)
# Shapiro test
shapiro.test(res.nbh.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.nbh.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.nbh.tpcb))[1, 'R2c']

# Extract coefficient values
time.coeff <- summary(lmem.nbh.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lmem.nbh.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -log(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# (3.2) log.tPCB vs. time + season + flow + temp + site (nbh.log.tpcb.2)
log.tpcb <- nbh.log.tpcb.2$logtPCB
time <- nbh.log.tpcb.2$time
site <- nbh.log.tpcb.2$site.code
season <- nbh.log.tpcb.2$season
flow <- nbh.log.tpcb.2$flow
tem <- nbh.log.tpcb.2$temp

lmem.nbh.log.tpcb <- lmer(log.tpcb ~ 1 + time + season + season + flow + tem + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.nbh.log.tpcb)
# Look at residuals
res.nbh.log.tpcb <- resid(lmem.nbh.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.nbh.log.tpcb, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.nbh.log.tpcb)
# Shapiro test
shapiro.test(res.nbh.log.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res.nbh.log.tpcb, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.nbh.log.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.nbh.log.tpcb))[1, 'R2c']

# Modeling plots
# (1) Get predicted values tpcb
fit.values.nbh.tpcb <- as.data.frame(fitted(lmem.nbh.tpcb))
# Add column name
colnames(fit.values.nbh.tpcb) <- c("predicted")
# Add predicted values to data.frame
nbh.tpcb.2$predicted <- 10^(fit.values.nbh.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(nbh.tpcb.2, aes(x = tPCB, y = predicted)) +
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

ggplot(nbh.tpcb.2, aes(x = tPCB, y = predicted)) +
  geom_point() +
  scale_x_log10(limits = c(10, 1e4)) +
  scale_y_log10(limits = c(10, 1e4)) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 1.3) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 25, y = 10000,
           label = 'nbh River', colour = 'black', size = 4,
           fontface = 2)

# Plot residuals vs. predictions
plot(log10(nbh.tpcb.2$predicted), res.nbh.tpcb)
abline(0, 0)

# (2) Get predicted values log.tpcb
fit.values.nbh.log.tpcb <- as.data.frame(fitted(lmem.nbh.log.tpcb))
# Add column name
colnames(fit.values.nbh.log.tpcb) <- c("predicted")
# Add predicted values to data.frame
nbh.log.tpcb.2$predicted <- fit.values.nbh.log.tpcb$predicted

# Plot prediction vs. observations, 1:1 line
ggplot(nbh.log.tpcb.2, aes(x = logtPCB, y = predicted)) +
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

ggplot(nbh.log.tpcb.2, aes(x = logtPCB, y = predicted)) +
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
plot(nbh.log.tpcb.2$predicted, res.nbh.log.tpcb)
abline(0, 0)
