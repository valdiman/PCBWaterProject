## Pre water PCB concentrations data analysis

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

# Load libraries
library(ggplot2)
library(scales) # function trans_breaks
#library(gridExtra)
#library(tidyverse)
library(stringr) # str_detect
library(robustbase) # function colMedians
library(dplyr) # performs %>%
library(tibble) # adds a column
library(lme4) # performs lme
library(MuMIn) # gets Rs from lme
library(lmerTest) # gets the p-value from lme
library(zoo) # yields seasons

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc.0 <- read.csv("WaterDataCongenerAroclor08052022.csv")

# (I) Entire data ---------------------------------------------------------
# Prepare data -----------------------------------------------------------
# (1) Calculate total PCB per sample
tpcb <- rowSums(wdc.0[, c(12:115)], na.rm = T)
# Change date format
wdc.0$SampleDate <- as.Date(wdc.0$SampleDate, format = "%m/%d/%y")
# Calculate sampling time
time.day <- data.frame(as.Date(wdc.0$SampleDate) - min(as.Date(wdc.0$SampleDate)))
# Include season
yq <- as.yearqtr(as.yearmon(wdc.0$SampleDate, "%m/%d/%Y") + 1/12)
season <- factor(format(yq, "%q"), levels = 1:4,
                 labels = c("0", "s-1", "s-2", "s-3")) # winter, spring, summer, fall
# Create data frame
tpcb <- cbind(wdc.0$SiteName, wdc.0$SampleDate,
              as.matrix(tpcb), data.frame(time.day), season)
# Add names to columns
colnames(tpcb) <- c("Site", "date", "tPCB", "time", "season")

# (2) Calculate total log PCB per sample
# Remove metadata
log.pcb <- subset(wdc.0, select = -c(ID:AroclorCongener))
# Remove Aroclor data
log.pcb <- subset(log.pcb, select = -c(A1016:A1260))
# Log 10 individual PCBs 
log.pcb <- log10(log.pcb)
# Replace -inf to NA
t.log.pcb <- do.call(data.frame,
               lapply(log.pcb,
                      function(x) replace(x, is.infinite(x), NA)))
# Sum individual log 10 PCBs
t.log.pcb <- rowSums(t.log.pcb, na.rm = T)
# Generate data.frame for analysis and plots
log.tpcb <- cbind(wdc.0$SiteName, wdc.0$SampleDate,
                  as.matrix(t.log.pcb), data.frame(time.day), season)
colnames(log.tpcb) <- c("Site", "date", "logtPCB", "time", "season")

# General plots -------------------------------------------------------
# (1) Histograms
# (1.1) tpcb
hist(tpcb$tPCB)
hist(log10(tpcb$tPCB))
# (1.2) log.tPCB
hist(log.tpcb$logtPCB)
hist(log10(log.tpcb$logtPCB))

# (2) One box plot
# (2.1) tPCB
ggplot(tpcb, aes(x = "", y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  theme(aspect.ratio = 14/2) +
  xlab(expression(bold(Sigma*"PCB (n = 6593)")))+
  ylab(expression(bold("Water Concentration 1990 - 2020 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10,
                                    angle = 45, hjust = 1.8,
                                    vjust = 2)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotation_logticks(sides = "l")

# (2.2) log.tPCB
ggplot(log.tpcb, aes(x = "", y = logtPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  theme(aspect.ratio = 14/2) +
  xlab(expression(bold(Sigma*"PCB (n = 6593)")))+
  ylab(expression(bold("Water Concentration 1990 - 2020 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10,
                                    angle = 45, hjust = 1.8,
                                    vjust = 2)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotation_logticks(sides = "l")

# (3) Box plots per site
# (3.1) tPCB
ggplot(tpcb, aes(x = factor(Site), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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
ggplot(log.tpcb, aes(x = factor(Site), y = logtPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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

# (4) Time series
# (4.1) tPCB
ggplot(tpcb, aes(x = format(date,'%Y'), y = tPCB)) +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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
ggplot(log.tpcb, aes(x = format(date,'%Y'), y = logtPCB)) +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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

# (5) Seasonality
# (5.1) tPCB
ggplot(tpcb, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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

# (5.2) log.tPCB
ggplot(log.tpcb, aes(x = season, y = logtPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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

# General regressions -----------------------------------------------------
# (1) Perform linear regression (lr)
# (1.1) tPCB + 1 vs. time
lr.tpcb <- lm(log10(tPCB + 1) ~ time, data = tpcb)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. time
lr.tpcb <- lm(logtPCB ~ time, data = log.tpcb)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB + 1 vs. season
lr.tpcb <- lm(log10(tPCB + 1) ~ time + season, data = tpcb)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. season
lr.tpcb <- lm(logtPCB ~ season, data = log.tpcb)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (II) Data per site name -------------------------------------------------------
# Select site name
Blu <- wdc.0[str_detect(wdc.0$SiteName, 'BlueRiver'),]
Fox <- wdc.0[str_detect(wdc.0$SiteName, 'FoxRiver'),]
Hud <- wdc.0[str_detect(wdc.0$SiteName, 'HudsonRiver'),]
Hou <- wdc.0[str_detect(wdc.0$SiteName, 'HousatonicRiver'),]
Kal <- wdc.0[str_detect(wdc.0$SiteName, 'KalamazooRiver'),]
Mde <- wdc.0[str_detect(wdc.0$SiteName, 'MDE_TMDL'),] # Check
NB <- wdc.0[str_detect(wdc.0$SiteName, 'NewBedford'),] # Check
Por <- wdc.0[str_detect(wdc.0$SiteName, 'PortlandHarbor'),] # Check
Spo <- wdc.0[str_detect(wdc.0$SiteName, 'SpokaneRiver'),]

# (1) Calculate total PCB
d <- Kal
tpcb.s <- rowSums(d[, c(12:115)], na.rm = T)
# Change date format
d$SampleDate <- as.Date(d$SampleDate, format = "%m/%d/%y")
# Calculate sampling time
time.day <- data.frame(as.Date(d$SampleDate) - min(as.Date(d$SampleDate)))
# Create individual code for each site sampled
site.numb <- d$SiteSampled %>% as.factor() %>% as.numeric
# Include season
yq.s <- as.yearqtr(as.yearmon(d$SampleDate, "%m/%d/%Y") + 1/12)
season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                 labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
# Create data frame
tpcb.s <- cbind(factor(d$SiteSampled), d$SampleDate,
                as.matrix(tpcb.s), data.frame(time.day),
                site.numb, season.s)
# Add column manes
colnames(tpcb.s) <- c("site", "date", "tPCB", "time",
                      "site.code", "season")

# (2) Calculate total log PCB
# Remove metadata
log.pcb.s <- subset(d, select = -c(ID:AroclorCongener))
# Remove Aroclor data
log.pcb.s <- subset(log.pcb.s, select = -c(A1016:A1260))
# Log 10 individual PCBs 
log.pcb.s <- log10(log.pcb.s)
# Replace -inf to NA
t.log.pcb.s <- do.call(data.frame,
               lapply(log.pcb.s,
                      function(x) replace(x, is.infinite(x), NA)))
# Sum individual log 10 PCBs
t.log.pcb.s <- rowSums(t.log.pcb.s, na.rm = T)
# Generate data.frame for analysis and plots
log.tpcb.s <- cbind(factor(d$SiteSampled), d$SampleDate,
                    as.matrix(t.log.pcb.s), data.frame(time.day),
                    site.numb, season.s)
colnames(log.tpcb.s) <- c("site", "date", "logtPCB", "time",
                          "site.code", "season")

# Plots --------------------------------------------------------------
# (1) Histograms

tpcb.s <- read.csv("KalV2.csv")
# (1.1) tPCB
hist(tpcb.s$tPCB)
hist(log10(tpcb.s$tPCB))
# (1.2) log.tPCB
hist(log.tpcb.s$logtPCB)
hist(log10(log.tpcb.s$logtPCB))

# (2) One box plot, shows sampled sites
# (2.1) tPCB
ggplot(tpcb.s, aes(x = "", y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  geom_point(aes(colour = factor(site))) +
  theme(aspect.ratio = 14/2) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotation_logticks(sides = "l")

# (2.2) log.tPCB
ggplot(log.tpcb.s, aes(x = "", y = logtPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  geom_point(aes(colour = factor(site))) +
  theme(aspect.ratio = 14/2) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotation_logticks(sides = "l")

# (3) Time trend
# (3.1) tPCB (add +1 to avoid error with zeros)
ggplot(tpcb.s, aes(y = tPCB,
                   x = format(date,'%Y'))) +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(Sigma*"PCB concentration (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 7,
                                   angle = 45, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 10)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# (3.2) log.tPCB
ggplot(log.tpcb.s, aes(y = logtPCB,
                   x = format(date,'%Y'))) +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(Sigma*"PCB concentration (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 7,
                                   angle = 45, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 10)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# (4) Time trend plots, shows sampled sites
# (4.1) tPCB
ggplot(tpcb.s, aes(y = tPCB,
                   x = format(date,'%Y'))) +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  geom_point(aes(colour = factor(site)))

# (4.2) log.tPCB
ggplot(log.tpcb.s, aes(y = logtPCB,
                   x = format(date,'%Y'))) +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  geom_point(aes(colour = factor(site)))

# (5) Seasonality
# (5.1) tPCB
ggplot(tpcb.s, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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

# (5.2) log.tPCB
ggplot(log.tpcb.s, aes(x = season, y = logtPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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

# Regressions -------------------------------------------------------------

# (1) Perform linear regression (lr)
# (1.1) tPCB + 1 vs. time
lr.tpcb.s <- lm(log10(tPCB) ~ time, data = tpcb.s)
# See results
summary(lr.tpcb.s)
# Look at residuals
res <- resid(lr.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. time
lr.logtpcb.s <- lm(logtPCB ~ time, data = log.tpcb.s)
# See results
summary(lr.logtpcb.s)
# Look at residuals
res <- resid(lr.logtpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB + 1 vs. season
lr.tpcb <- lm(log10(tPCB + 0) ~ season, data = tpcb.s)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. season
lr.tpcb <- lm(logtPCB ~ season, data = log.tpcb.s)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.5) tPCB + 1 vs. time + season
lr.tpcb <- lm(log10(tPCB + 0) ~ time + season, data = tpcb.s)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.6) log.tPCB vs. time + season 
lr.tpcb <- lm(logtPCB ~ time + season, data = log.tpcb.s)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.7) tPCB + 1 vs flow and temp
par <- read.csv("USGSData/FoxRiverflowtemp2.csv")
par$watertemp <- as.numeric(par$watertemp)
par$watertemp2 <- as.numeric(par$watertemp2)
lr.tpcb <- lm(log10(tPCB + 0) ~ time + season + X4108660, data = tpcb.s)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) Perform Linear Mixed-Effects Model (LME)
tpcb <- tpcb.s$tPCB
time <- tpcb.s$time
site <- tpcb.s$site.code
season <- tpcb.s$season
flow <- tpcb.s$X4108660
wtemp <- par$watertemp2

# (2.1) tPCB
lmem.tpcb <- lmer(log10(tpcb) ~ 1 + time + flow + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.tpcb)
# Look at residuals
res <- resid(lmem.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res, main = "log10(C + 1)")
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2.2) log.tPCB
lmem.tpcb <- lmer(log.tpcb.s$logtPCB ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.tpcb)
# Look at residuals
res <- resid(lmem.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res, main = "log10(C + 1)")
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (III) Analysis per specific site ----------------------------------------------
# Check number of samples from the same sampling site
table(tpcb.s$site)
# Select site with 30 or > sampler per site
# From Fox: LakeWinne
# From Hud: BakersFalls, RogersIsland, Schuylerville, Stillwater,
# ThompsonIslandDam, Waterford1, Waterford2
# From Hou: DawesAveBridge, DivisionStBridge, DownstreamofLymanStBridge,
# ElmStBridge, FirstPomeroyAveBridge, FormerHousatonicStAbutment,
# HeadwatersofWoodsPond, HolmesRdBridge, HubbardAveBridge,
# LymanStBridge, NewellStBridge, NewellStParkingLotFootbridge,
# NewLenoxRdBridge, PomeroyAveBridge, SchweitzerBridge, SilverLakeOutlet,
# UpstreamofNewellStBridge
# From Kal: 26THST, 2NDST, 58THST, ALCOTTST, BLUESTARHWY, CORKST, D.AVENUE,
# E.MICHIGANAVE, FARMERST, FORMERPLAINWELLIMPOUNDMENT, GIBSONST-PORTAGECREEK,
# LAKEALLEGANDAM, M222, M89PLAINWELL, MICHIGANAVE, PLAINWELLDAM, PORTAGECREEK,
# RIVERST, ROSSMANCREEK

# Fox River
Fox <- wdc.0[str_detect(wdc.0$SiteName, 'FoxRiver'),]
Fox.1 <- Fox[str_detect(Fox$SiteSampled, 'LakeWinne'),]
# Hudson River
Hud <- wdc.0[str_detect(wdc.0$SiteName, 'HudsonRiver'),]
Hud.1 <- Hud[str_detect(Hud$SiteSampled, 'BakersFalls'),]
Hud.2 <- Hud[str_detect(Hud$SiteSampled, 'RogersIsland'),] # Check
Hud.3 <- Hud[str_detect(Hud$SiteSampled, 'Schuylerville'),] # Check, but only 2 years
Hud.4 <- Hud[str_detect(Hud$SiteSampled, 'Stillwater'),]
Hud.5 <- Hud[str_detect(Hud$SiteSampled, 'ThompsonIslandDam'),]
Hud.6 <- Hud[str_detect(Hud$SiteSampled, 'Waterford1'),] # Check
Hud.7 <- Hud[str_detect(Hud$SiteSampled, 'Waterford2'),]
# Housatonic River
Hou <- wdc.0[str_detect(wdc.0$SiteName, 'HousatonicRiver'),]
Hou.1 <- Hou[str_detect(Hou$SiteSampled, 'DawesAveBridge'),]
Hou.2 <- Hou[str_detect(Hou$SiteSampled, 'DivisionStBridge'),] # Check, 2 outliers
Hou.3 <- Hou[str_detect(Hou$SiteSampled, 'DownstreamofLymanStBridge'),]
Hou.4 <- Hou[str_detect(Hou$SiteSampled, 'ElmStBridge'),]
Hou.5 <- Hou[str_detect(Hou$SiteSampled, 'FirstPomeroyAveBridge'),]
Hou.6 <- Hou[str_detect(Hou$SiteSampled, 'FormerHousatonicStAbutment'),] # Check 
Hou.7 <- Hou[str_detect(Hou$SiteSampled, 'HeadwatersofWoodsPond'),]
Hou.8 <- Hou[str_detect(Hou$SiteSampled, 'HolmesRdBridge'),] # Check
Hou.9 <- Hou[str_detect(Hou$SiteSampled, 'HubbardAveBridge'),] # Check
Hou.10 <- Hou[str_detect(Hou$SiteSampled, 'LymanStBridge'),]
Hou.11 <- Hou[str_detect(Hou$SiteSampled, 'NewellStBridge'),]
Hou.12 <- Hou[str_detect(Hou$SiteSampled, 'NewellStParkingLotFootbridge'),]
Hou.13 <- Hou[str_detect(Hou$SiteSampled, 'NewLenoxRdBridge'),] # Check
Hou.14 <- Hou[str_detect(Hou$SiteSampled, 'PomeroyAveBridge'),] # Check
Hou.15 <- Hou[str_detect(Hou$SiteSampled, 'SchweitzerBridge'),]
Hou.16 <- Hou[str_detect(Hou$SiteSampled, 'SilverLakeOutlet'),]
Hou.17 <- Hou[str_detect(Hou$SiteSampled, 'UpstreamofNewellStBridge'),]
# Kalamazoo River
Kal <- wdc.0[str_detect(wdc.0$SiteName, 'KalamazooRiver'),]
Kal.1 <- wdc.0[str_detect(Kal$SiteSampled, '26THST'),]
Kal.2 <- wdc.0[str_detect(Kal$SiteSampled, '2NDST'),]
Kal.3 <- wdc.0[str_detect(Kal$SiteSampled, '58THST'),]
Kal.4 <- wdc.0[str_detect(Kal$SiteSampled, 'ALCOTTST'),]
Kal.5 <- wdc.0[str_detect(Kal$SiteSampled, 'BLUESTARHWY'),]
Kal.6 <- wdc.0[str_detect(Kal$SiteSampled, 'CORKST'),]
Kal.7 <- wdc.0[str_detect(Kal$SiteSampled, 'D.AVENUE'),]
Kal.8 <- wdc.0[str_detect(Kal$SiteSampled, 'E.MICHIGANAVE'),]
Kal.9 <- wdc.0[str_detect(Kal$SiteSampled, 'FARMERST'),]
Kal.10 <- wdc.0[str_detect(Kal$SiteSampled, 'FORMERPLAINWELLIMPOUNDMENT'),]
Kal.11 <- wdc.0[str_detect(Kal$SiteSampled, 'GIBSONST-PORTAGECREEK'),]
Kal.12 <- wdc.0[str_detect(Kal$SiteSampled, 'LAKEALLEGANDAM'),]
Kal.13 <- wdc.0[str_detect(Kal$SiteSampled, 'M222'),]
Kal.14 <- wdc.0[str_detect(Kal$SiteSampled, 'M89PLAINWELL'),]
Kal.15 <- wdc.0[str_detect(Kal$SiteSampled, 'MICHIGANAVE'),]
Kal.16 <- wdc.0[str_detect(Kal$SiteSampled, 'PLAINWELLDAM'),]
Kal.17 <- wdc.0[str_detect(Kal$SiteSampled, 'PORTAGECREEK'),]
Kal.18 <- wdc.0[str_detect(Kal$SiteSampled, 'RIVERST'),]
Kal.19 <- wdc.0[str_detect(Kal$SiteSampled, 'ROSSMANCREEK'),]

# (1) Calculate total PCB per sample per site
d <- Hud.2
d.1 <- rowSums(d[, c(12:115)], na.rm = T)
# Change date format
d$SampleDate <- as.Date(d$SampleDate, format = "%m/%d/%y")
# Calculate sampling time
time.day <- data.frame(as.Date(d$SampleDate) - min(as.Date(d$SampleDate)))
# Include season
yq.s.s <- as.yearqtr(as.yearmon(d$SampleDate, "%m/%d/%Y") + 1/12)
season.s.s <- factor(format(yq.s.s, "%q"), levels = 1:4,
                     labels = c("0", "s-1", "s-2", "s-3")) # winter, spring, summer, fall
# Create data frame
tpcb.s.s <- cbind(factor(d$SiteSampled), d$SampleDate,
                as.matrix(d.1), data.frame(time.day), season.s.s)
# Add column manes
colnames(tpcb.s.s) <- c("site", "date", "tPCB", "time", "season")

# (2) Calculate total log PCB per sample
# Remove metadata
d.2 <- subset(d, select = -c(ID:AroclorCongener))
# Remove Aroclor data
d.3 <- subset(d.2, select = -c(A1016:A1260))
# Log 10 individual PCBs 
d.4 <- log10(d.3)
# Replace -inf to NA
d.4 <- do.call(data.frame, lapply(d.4,
               function(x) replace(x, is.infinite(x), NA)))
# Sum individual log 10 PCBs
d.5 <- rowSums(d.4, na.rm = T)
# Generate data.frame for analysis and plots
log.tpcb.s.s <- cbind(factor(d$SiteSampled), d$SampleDate,
                    as.matrix(d.5), data.frame(time.day), season.s.s)
colnames(log.tpcb.s.s) <- c("site", "date", "logtPCB", "time", "season")

# Plots -------------------------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(tpcb.s.s$tPCB)
hist(log10(tpcb.s.s$tPCB))
# (1.2) log.tPCB
hist(log.tpcb.s.s$logtPCB)
hist(log10(log.tpcb.s.s$logtPCB))

# (2) Time trend plots
# (2.1) tPCB
ggplot(tpcb.s.s, aes(y = tPCB,
                   x = format(date,'%Y'))) +
  geom_point() +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15)

# (2.2) log.tPCB
ggplot(log.tpcb.s.s, aes(y = logtPCB,
                   x = format(date,'%Y'))) +
  geom_point() +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15)

# (3) Seasonality
# (3.1) tPCB
ggplot(tpcb.s.s, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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
ggplot(log.tpcb.s.s, aes(x = season, y = logtPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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

# Regressions -------------------------------------------------------------
# (1) Perform linear regression (lr)
# (1.1) tPCB + 1 vs. time
lr.tpcb.s.s <- lm(log10(tPCB + 0) ~ time, data = tpcb.s.s)
# See results
summary(lr.tpcb.s.s)
# Look at residuals
res <- resid(lr.tpcb.s.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. time
lr.logtpcb.s.s <- lm(logtPCB ~ time, data = log.tpcb.s.s)
# See results
summary(lr.logtpcb.s.s)
# Look at residuals
res <- resid(lr.logtpcb.s.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB + 1 vs. season
lr.tpcb <- lm(log10(tPCB + 0) ~ season, data = tpcb.s.s)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. season
lr.tpcb <- lm(logtPCB ~ season, data = log.tpcb.s.s)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.5) tPCB + 1 vs. time + season
lr.tpcb <- lm(log10(tPCB + 0) ~ time + season, data = tpcb.s.s)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.6) log.tPCB vs. time + season 
lr.tpcb <- lm(logtPCB ~ time + season, data = log.tpcb.s.s)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.7) tPCB + 1 vs. time. season, flow
# Read flow data
flow <- read.csv("USGSData/flow.csv")


lr.tpcb <- lm(log10(tPCB + 0) ~ time + season + flow, data = tpcb.s.s)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
