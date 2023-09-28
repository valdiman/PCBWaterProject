## Water PCB concentrations data analysis
## Fox River 2005 - 2018
## This code version has all the preliminary analyses,
## including sum of log10 PCB and linear and multiple linear regression

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

# Select Fox River data ---------------------------------------------------
fox.0 <- wdc[str_detect(wdc$LocationName, 'Fox River'),]
# Lake Winnebago is a background site.
# Data preparation --------------------------------------------------------
{
  # Remove samples (rows) with total PCBs  = 0
  fox.1 <- fox.0[!(rowSums(fox.0[, c(14:117)], na.rm = TRUE)==0),]
  # Calculate total PCB
  tpcb.fox <- rowSums(fox.1[, c(14:117)], na.rm = T)
  # Calculate total log PCB
  # Remove metadata
  fox.log <- subset(fox.1, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  fox.log <- subset(fox.log, select = -c(A1016:A1260))
  # Log10 individual PCBs 
  fox.log <- log10(fox.log)
  # Replace -inf to NA
  fox.log <- do.call(data.frame,
                   lapply(fox.log,
                          function(x) replace(x, is.infinite(x), NA)))
  # Sum individual log 10 PCBs
  fox.log.tpcb <- rowSums(fox.log, na.rm = T)
  # Change date format
  fox.1$SampleDate <- as.Date(fox.1$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(fox.1$SampleDate) - min(as.Date(fox.1$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- fox.1$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(fox.1$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  fox.tpcb <- cbind(factor(fox.1$SiteID), fox.1$SampleDate,
                    fox.1$Latitude, fox.1$Longitude, as.matrix(tpcb.fox),
                    as.matrix(fox.log.tpcb), data.frame(time.day),
                    site.numb, season.s)
  # Add column names
  colnames(fox.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "logtPCB", "time", "site.code", "season")
}

# Get coordinates per site to plot in Google Earth
fox.location <- fox.tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
fox.location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                            data = fox.location, mean)

# General plots -------------------------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(fox.tpcb$tPCB)
hist(log10(fox.tpcb$tPCB)) # Better approach
# (1.2) log.tPCB
hist(fox.tpcb$logtPCB)

# (2) Time trend plots
# (2.1) tPCB
ggplot(fox.tpcb, aes(y = tPCB,
                     x = format(date,'%Y'))) +
  geom_point(shape = 21, size = 2, fill = "#66ccff") +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9)) +
  annotate("text", x = 5.8, y = 10^5, label = "Fox River",
           size = 3)
  

# (2.2) log.tPCB
ggplot(fox.tpcb, aes(y = logtPCB,
                         x = format(date,'%Y'))) +
  geom_point(shape = 21, size = 2, fill = "#66ccff") +
  xlab("") +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9)) +
  annotate("text", x = 5.8, y = 90, label = "Fox River",
           size = 3)

# (3) Seasonality
# (3.1) tPCB
ggplot(fox.tpcb, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# (3.2) log.tPCB
ggplot(fox.tpcb, aes(x = season, y = logtPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# (4) Sites
# (4.1) tPCB
ggplot(fox.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# (4.2) log.tPCB
ggplot(fox.tpcb, aes(x = factor(SiteID), y = logtPCB)) + 
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# Remove site -------------------------------------------------------------
# Remove site Lake Winnebago (background site)
fox.tpcb.2 <- subset(fox.tpcb, SiteID != c("WCPCB-FOX001"))

# Plots w/o Lake Winnebago ------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(fox.tpcb.2$tPCB)
hist(log10(fox.tpcb.2$tPCB))
# (1.2) log.tPCB
hist(fox.tpcb.2$logtPCB)

# (2) Time trend plots
# (2.1) tPCB
ggplot(fox.tpcb.2, aes(y = tPCB,
                     x = format(date,'%Y'))) +
  geom_point(shape = 21, fill = "#66ccff") +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9)) +
  annotate("text", x = 5.8, y = 10^5, label = "Fox River",
           size = 3)

# (2.2) log.tPCB
ggplot(fox.tpcb.2, aes(y = logtPCB,
                         x = format(date,'%Y'))) +
  geom_point(shape = 21, fill = "#66ccff") +
  xlab("") +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9)) +
  annotate("text", x = 5.8, y = 90, label = "Fox River",
           size = 3)

# (3) Seasonality
# (3.1) tPCB
ggplot(fox.tpcb.2, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# (3.2) log.tPCB
ggplot(fox.tpcb.2, aes(x = season, y = logtPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# (4) Sites
# (4.1) tPCB
ggplot(fox.tpcb.2, aes(x = factor(SiteID), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# (4.2) log.tPCB
ggplot(fox.tpcb.2, aes(x = factor(SiteID), y = logtPCB)) + 
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# Include USGS flow data --------------------------------------------------
{
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
  # Add USGS data to fox.tpcb, matching dates, conversion to m3/s
  fox.tpcb.2$flow <- 0.03*flow$X_.Primary.Stream.Flow._00060_00003[match(fox.tpcb.2$date,
                                                                         flow$Date)]
  fox.tpcb.2$temp <- 273.15 + temp$X_00010_00003[match(fox.tpcb.2$date,
                                                       temp$Date)]
  # Remove samples with temp = NA
  fox.tpcb.2 <- na.omit(fox.tpcb.2)
}

# tPCB Regressions --------------------------------------------------------
# (1) Perform linear regression (lr)
# (1.1) tPCB vs. time
lr.fox.tpcb.t <- lm(log10(tPCB) ~ time, data = fox.tpcb.2)
# See results
summary(lr.fox.tpcb.t)
# Look at residuals
{
  res <- resid(lr.fox.tpcb.t) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. time
lr.fox.log.tpcb.t <- lm(logtPCB ~ time, data = fox.tpcb.2)
# See results
summary(lr.fox.log.tpcb.t)
# Look at residuals
{
  res <- resid(lr.fox.log.tpcb.t) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB vs. season
lr.fox.tpcb.s <- lm(log10(tPCB) ~ season, data = fox.tpcb.2)
# See results
summary(lr.fox.tpcb.s)
# Look at residuals
{
  res <- resid(lr.fox.tpcb.s) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. season
lr.fox.log.tpcb.s <- lm(logtPCB ~ season, data = fox.tpcb.2)
# See results
summary(lr.fox.log.tpcb.s)
# Look at residuals
{
  res <- resid(lr.fox.log.tpcb.s) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.5) tPCB vs. flow
lr.fox.tpcb.f <- lm(log10(tPCB) ~ flow, data = fox.tpcb.2)
# See results
summary(lr.fox.tpcb.f)
# Look at residuals
{
  res <- resid(lr.fox.tpcb.f) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.6) log.tPCB vs. flow
lr.fox.log.tpcb.f <- lm(logtPCB ~ flow, data = fox.tpcb.2)
# See results
summary(lr.fox.log.tpcb.f)
# Look at residuals
{
  res <- resid(lr.fox.log.tpcb.f) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.7) tPCB vs. water temperature
lr.fox.tpcb.te <- lm(log10(tPCB) ~ temp, data = fox.tpcb.2)
# See results
summary(lr.fox.tpcb.te)
# Look at residuals
{
  res <- resid(lr.fox.tpcb.te) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.8) log.tPCB vs. temperature
lr.fox.log.tpcb.te <- lm(logtPCB ~ temp, data = fox.tpcb.2)
# See results
summary(lr.fox.log.tpcb.te)
# Look at residuals
{
  res <- resid(lr.fox.log.tpcb.te) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) MLR
# (2.1) tPCB vs. time + season + flow + temp
mlr.fox.tpcb <- lm(log10(tPCB) ~ time + season + flow + temp, data = fox.tpcb.2)
# See results
summary(mlr.fox.tpcb)
# Look at residuals
{
  res <- resid(mlr.fox.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
# Predictions
fit.mlr.values.fox.tpcb <- as.data.frame(predict(mlr.fox.tpcb))

# (2.2) log.tPCB vs. time + season + flow + temp. Best approach
mlr.fox.log.tpcb <- lm(logtPCB ~ time + season + flow + temp,
                       data = fox.tpcb.2)
# See results
summary(mlr.fox.log.tpcb)
# Look at residuals
{
  res <- resid(mlr.fox.log.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (3) Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- fox.tpcb.2$tPCB
log.tpcb <- fox.tpcb.2$logtPCB
time <- fox.tpcb.2$time
site <- fox.tpcb.2$site.code
season <- fox.tpcb.2$season
flow <- fox.tpcb.2$flow
tem <- fox.tpcb.2$temp
# (3.1) tPCB vs. time + season + flow + temp + site
lmem.fox.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + flow + tem + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.fox.tpcb)
# Look at residuals
{
  res.fox.tpcb <- resid(lmem.fox.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.fox.tpcb, main = "log10(C)")
  qqnorm(res.fox.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.fox.tpcb)
}
# Shapiro test
shapiro.test(res.fox.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
# Random effect site Std Dev
RandonEffectSiteStdDev <- as.data.frame(VarCorr(lmem.fox.tpcb))[1,'sdcor']
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.fox.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.fox.tpcb))[1, 'R2c']

# Extract coefficient values
time.coeff <- summary(lmem.fox.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lmem.fox.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -ln(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# (3.2) log.tPCB vs. time + season + flow + temp + site
lmem.fox.log.tpcb <- lmer(log.tpcb ~ 1 + time + season + flow + tem + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.fox.log.tpcb)
# Look at residuals
{
  res.fox.log.tpcb <- resid(lmem.fox.log.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.fox.log.tpcb,
         main = expression(paste("Normal Q-Q Plot", " (", Sigma,
                                 "log"[10]*"PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.fox.log.tpcb)
}
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
fit.lme.values.fox.tpcb <- as.data.frame(fitted(lmem.fox.tpcb))
# Add column name
colnames(fit.lme.values.fox.tpcb) <- c("predicted")
# Add predicted values to data.frame
fox.tpcb.2$predicted <- 10^(fit.lme.values.fox.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(fox.tpcb.2, aes(x = tPCB, y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "#66ccff") +
  scale_y_log10(limits = c(10, 10^5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10, 10^5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", linewidth = 1.3) +
  geom_abline(intercept = 0.5, slope = 1, col = "blue", linewidth = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = -0.5, slope = 1, col = "blue", linewidth = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 100, y = 10^4.8,
           label = expression("Fox River (R"^2*"= 0.78)"),
           size = 4, fontface = 2)

# Plot residuals vs. predictions
{
  plot(log10(fox.tpcb.2$predicted), res.fox.tpcb,
       points(log10(fox.tpcb.2$predicted), res.fox.tpcb, pch = 16, 
              col = "#66ccff"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  abline(0, 0)
  abline(h = c(-1, 1), col = "grey")
  abline(v = seq(2, 3.5, 0.5), col = "grey")
  }

# (2) Get predicted values log.tpcb
fit.lme.values.fox.log.tpcb <- as.data.frame(fitted(lmem.fox.log.tpcb))
# Add column name
colnames(fit.lme.values.fox.log.tpcb) <- c("predicted")
# Add predicted values to data.frame
fox.tpcb.2$predictedlog <- fit.lme.values.fox.log.tpcb$predicted

# Plot prediction vs. observations, 1:1 line
ggplot(fox.tpcb.2, aes(x = logtPCB, y = predictedlog)) +
  geom_point(shape = 21, size = 2, fill = "#66ccff") +
  scale_x_continuous(limits = c(-2, 100)) +
  scale_y_continuous(limits = c(-2, 100)) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 9),
        axis.title.x = element_text(face = "bold", size = 9)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  geom_abline(intercept = 0, slope = 1, col = "red", linewidth = 1.3) +
  geom_abline(intercept = 2, slope = 1, col = "blue", linewidth = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = -2, slope = 1, col = "blue", linewidth = 0.8) # 2:1 line (factor of 2)

# Plot residuals vs. predictions
plot(fox.tpcb.2$predictedlog, res.fox.log.tpcb)
abline(0, 0)

# Plot time series with mlr and lme predictions
# Create a data frame to storage data
{
  time.serie.tpcb <- as.data.frame(matrix(nrow = length(fox.tpcb.2[,1]),
                                          ncol = 4))
  # Add name to columns
  colnames(time.serie.tpcb) <- c('date', 'tPCB', 'mlrtPCB', 'lmetPCB')
  # Add data
  time.serie.tpcb$date <- fox.tpcb.2$date
  time.serie.tpcb$tPCB <- fox.tpcb.2$tPCB
  time.serie.tpcb$mlrtPCB <- 10^(fit.mlr.values.fox.tpcb)
  time.serie.tpcb$lmetPCB <- 10^(fit.lme.values.fox.tpcb)
  # Change again the names
  colnames(time.serie.tpcb[,3]) <- c("mlrtPCB")
  colnames(time.serie.tpcb[,4]) <- c("lmetPCB")
  # Change data.frame format to be plotted
  time.serie.tpcb.2 <- melt(time.serie.tpcb, id.vars = c("date"))
}
# Plot
ggplot(time.serie.tpcb.2, aes(x = date, y = value, group = variable)) +
  geom_point(aes(shape = variable, color = variable, size = variable,
                 fill = variable)) +
  scale_shape_manual(values = c(21, 3, 4)) +
  scale_color_manual(values = c('black','#8856a7', '#b2df8a')) +
  scale_size_manual(values = c(2, 1, 1)) +
  scale_fill_manual(values = c("#1b98e0", '#8856a7', '#b2df8a')) +
  scale_x_date(labels = date_format("%Y-%m")) +
  scale_y_log10(limits = c(10, 10000)) +
  xlab("") +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9)) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# Individual PCB Analysis -------------------------------------------------
# Use fox.1 (no 0s samples)
# Prepare data.frame
{
  fox.pcb <- subset(fox.1, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  fox.pcb <- subset(fox.pcb, select = -c(A1016:A1260))
  # Log10 individual PCBs 
  fox.pcb <- log10(fox.pcb)
  # Replace -inf to NA
  fox.pcb <- do.call(data.frame,
                     lapply(fox.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Add site ID
  fox.pcb$SiteID <- fox.1$SiteID
  # Change date format
  fox.pcb$SampleDate <- as.Date(fox.1$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  fox.pcb$time <- as.Date(fox.pcb$SampleDate) - min(as.Date(fox.pcb$SampleDate))
  # Create individual code for each site sampled
  fox.pcb$site.numb <- fox.pcb$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(fox.pcb$SampleDate, "%m/%d/%Y") + 1/12)
  fox.pcb$season <- factor(format(yq.s, "%q"), levels = 1:4,
                           labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Remove site Lake Winnebago (background site)
  fox.pcb.2 <- subset(fox.pcb, SiteID != c("WCPCB-FOX001"))
  # Include flow data from USGS station Fox River
  sitefoxN1 <- "04084445" # flow @ OX RIVER AT APPLETON, WI
  sitefoxN2 <- "040851385" # water temperature @ FOX RIVER AT OIL TANK DEPOT AT GREEN BAY, WI
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow <- readNWISdv(sitefoxN1, paramflow,
                     min(fox.pcb.2$SampleDate), max(fox.pcb.2$SampleDate))
  temp <- readNWISdv(sitefoxN2, paramtemp,
                     min(fox.pcb.2$SampleDate), max(fox.pcb.2$SampleDate))
  # Add USGS data to fox.pcb, matching dates, conversion to m3/s
  fox.pcb.2$flow <- 0.03*flow$X_.Primary.Stream.Flow._00060_00003[match(fox.pcb.2$SampleDate,
                                                                        flow$Date)]
  fox.pcb.2$temp <- 273.15 + temp$X_00010_00003[match(fox.pcb.2$SampleDate, temp$Date)]
  # Remove samples with temperature = NA
  fox.pcb.2 <- fox.pcb.2[!is.na(fox.pcb.2$temp), ]
  # Remove individual PCB that have 45% or less values
  fox.pcb.2[,colSums(is.na(fox.pcb.2)) > nrow(fox.pcb.2) - 45] <- NULL
  # Remove metadata
  fox.pcb.3 <- subset(fox.pcb.2, select = -c(SiteID:temp))
}

# Get covariates
time <- fox.pcb.2$time
flow <- fox.pcb.2$flow
temper <- fox.pcb.2$temp
season <- fox.pcb.2$season
site <- fox.pcb.2$site.numb

# MLR for individual PCBs  ------------------------------------------------
# Create matrix to storage results/coefficients
mlr.pcb <- matrix(nrow = length(fox.pcb.3), ncol = 22)

for(i in 1:length(fox.pcb.3)) {
  fit.mlr <- lm(fox.pcb.3[,i] ~ time + flow + temper + season)
  mlr.pcb[i,1] <- summary(fit.mlr)$coef[1,"Estimate"] # intercept
  mlr.pcb[i,2] <- summary(fit.mlr)$coef[1,"Std. Error"] # intercept error
  mlr.pcb[i,3] <- summary(fit.mlr)$coef[1,"Pr(>|t|)"] # intercept p-value
  mlr.pcb[i,4] <- summary(fit.mlr)$coef[2,"Estimate"] # time
  mlr.pcb[i,5] <- summary(fit.mlr)$coef[2,"Std. Error"] # time error
  mlr.pcb[i,6] <- summary(fit.mlr)$coef[2,"Pr(>|t|)"] # time p-value
  mlr.pcb[i,7] <- summary(fit.mlr)$coef[3,"Estimate"] # flow
  mlr.pcb[i,8] <- summary(fit.mlr)$coef[3,"Std. Error"] # flow error
  mlr.pcb[i,9] <- summary(fit.mlr)$coef[3,"Pr(>|t|)"] # flow p-value
  mlr.pcb[i,10] <- summary(fit.mlr)$coef[4,"Estimate"] # temperature
  mlr.pcb[i,11] <- summary(fit.mlr)$coef[4,"Std. Error"] # temperature error
  mlr.pcb[i,12] <- summary(fit.mlr)$coef[4,"Pr(>|t|)"] # temperature p-value
  mlr.pcb[i,13] <- summary(fit.mlr)$coef[5,"Estimate"] # season 2
  mlr.pcb[i,14] <- summary(fit.mlr)$coef[5,"Std. Error"] # season 2 error
  mlr.pcb[i,15] <- summary(fit.mlr)$coef[5,"Pr(>|t|)"] # season 2 p-value
  mlr.pcb[i,16] <- summary(fit.mlr)$coef[6,"Estimate"] # season 3
  mlr.pcb[i,17] <- summary(fit.mlr)$coef[6,"Std. Error"] # season 3 error
  mlr.pcb[i,18] <- summary(fit.mlr)$coef[6,"Pr(>|t|)"] # season 3 p-value
  mlr.pcb[i,19] <- -log(2)/mlr.pcb[i,4]/365 # t0.5
  mlr.pcb[i,20] <- abs(-log(2)/mlr.pcb[i,4]/365)*mlr.pcb[i,5]/abs(mlr.pcb[i,4]) # t0.5 error
  mlr.pcb[i,21] <- summary(fit.mlr)$adj.r.squared # R2 adj
  mlr.pcb[i,22] <- shapiro.test(resid(fit.mlr))$p.value
}

# Add column names
colnames(mlr.pcb) <- c("intercept", "intercep.error", "intercept.pv",
                       "time", "time.error", "time.pv", "flow", "flow.error",
                       "flow.pv", "temperature", "temperature.error", "temperature.pv",
                       "season2", "season2.error", "season2.pv", "season3",
                       "season3.error", "season3.pv", "t0.5", "t0.5.error",
                       "R2.adj", "Normality")
# Just 3 significant figures
mlr.pcb <- formatC(signif(mlr.pcb, digits = 3))
# Add congener names
congeners <- colnames(fox.pcb.3)
mlr.pcb <- cbind(congeners, mlr.pcb)
# Export results
write.csv(mlr.pcb, file = "Output/Data/csv/mlrFoxPCB.csv")

# LME for individual PCBs -------------------------------------------------
# Create matrix to store results
lme.pcb <- matrix(nrow = length(fox.pcb.3[1,]), ncol = 24)

# Perform LME
for (i in 1:length(fox.pcb.3[1,])) {
    fit <- lmer(fox.pcb.3[,i] ~ 1 + time + flow + temper + season + (1|site),
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
    lme.pcb[i,10] <- fixef(fit)[4] # temperature
    lme.pcb[i,11] <- summary(fit)$coef[4,"Std. Error"] # temperature error
    lme.pcb[i,12] <- summary(fit)$coef[4,"Pr(>|t|)"] # temperature p-value
    lme.pcb[i,13] <- fixef(fit)[5] # season 2
    lme.pcb[i,14] <- summary(fit)$coef[5,"Std. Error"] # season 2 error
    lme.pcb[i,15] <- summary(fit)$coef[5,"Pr(>|t|)"] # season 2 p-value
    lme.pcb[i,16] <- fixef(fit)[6] # season 3
    lme.pcb[i,17] <- summary(fit)$coef[6,"Std. Error"] # season 3 error
    lme.pcb[i,18] <- summary(fit)$coef[6,"Pr(>|t|)"] # season 3 p-value
    lme.pcb[i,19] <- -log(2)/lme.pcb[i,4]/365 # t0.5
    lme.pcb[i,20] <- abs(-log(2)/lme.pcb[i,4]/365)*lme.pcb[i,5]/abs(lme.pcb[i,4]) # t0.5 error
    lme.pcb[i,21] <- as.data.frame(VarCorr(fit))[1,'sdcor']
    lme.pcb[i,22] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
    lme.pcb[i,23] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
    lme.pcb[i,24] <- shapiro.test(resid(fit))$p.value
}

# Just 3 significant figures
lme.pcb <- formatC(signif(lme.pcb, digits = 3))
# Add PCB congener names
lme.pcb <- cbind(congeners, lme.pcb)
# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "flow", "flow.error", "flow.pv", "temperature",
                       "temperature.error", "temperature.pv", "season2",
                       "season2.error", "season2, pv", "season3",
                       "season3.error", "season3.pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality")

# Export results
write.csv(lme.pcb, file = "Output/Data/csv/LmeFoxPCB.csv")

# Get predicted values for selected PCBs
# tPCB vs. time + season + flow + temp
# (1) mlr
mlr.fox.pcbi <- lm(fox.pcb.3$PCB17 ~ time + season + flow + temp,
                   data = fox.tpcb.2)
# See results
summary(mlr.fox.pcbi)
# Look at residuals
{
  res.mlr.fox.pcbi <- resid(mlr.fox.pcbi) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
# Predictions
fit.mlr.values.fox.pcbi <- as.data.frame(predict(mlr.fox.pcbi))

# (2) lme
lme.fox.pcbi <- lmer(fox.pcb.3$PCB17 ~ 1 + time + flow + temper + season +
                       (1|site), REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                  check.nobs.vs.rankZ = "ignore",
                                  check.nobs.vs.nRE="ignore"))

# See results
summary(lme.fox.pcbi)
# Look at residuals
{
  res.lme.fox.pcbi <- resid(lme.fox.pcbi) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.lme.fox.pcbi, main = expression(paste("Normal Q-Q Plot PCB 17")))
  # Add a straight diagonal line to the plot
  qqline(res.lme.fox.pcbi)
}
# Shapiro test
shapiro.test(res.lme.fox.pcbi)
# One-sample Kolmogorov-Smirnov test
ks.test(res.lme.fox.pcbi, 'pnorm')

# (1) Get predicted values pcbi
date.pcbi <- format(fox.pcb.2$SampleDate, "%Y-%m-%d")
obs <- fox.pcb.3$PCB17
fox.pcbi <- cbind(date.pcbi, obs)
# Remove NA value from observations
fox.pcbi <- na.omit(fox.pcbi)
fit.mlr.values.pcbi <- as.data.frame(predict(mlr.fox.pcbi))
fit.lme.values.pcbi <- as.data.frame(fitted(lme.fox.pcbi))
fox.pcbi <- cbind(fox.pcbi, fit.mlr.values.pcbi,
                             fit.lme.values.pcbi)
colnames(fox.pcbi) <- c("date", "obs", "mlr", 'lme')
fox.pcbi.2$date <- as.Date(fox.pcbi.2$date)
fox.pcbi$obs <- as.numeric(fox.pcbi$obs)

# Plot residuals vs. predictions
# mlr
{
  plot(fox.pcbi$mlr, res.mlr.fox.pcbi,
       points(fox.pcbi$mlr, res.mlr.fox.pcbi,
              pch = 16, col = "#66ccff"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted concentration PCB 17 (pg/L)")),
       ylab = "Residual (mlr)")
  abline(0, 0)
  abline(h = seq(-2, 2, 1), col = "grey")
  abline(v = seq(0.8, 2, 0.2), col = "grey")
}
# lme
{
  plot(fox.pcbi$lme, res.lme.fox.pcbi,
       points(fox.pcbi$lme, res.lme.fox.pcbi, pch = 16, 
              col = "#66ccff"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted concentration PCB 17 (pg/L)")),
       ylab = "Residual (lme)")
  abline(0, 0)
  abline(h = seq(-2, 2, 1), col = "grey")
  abline(v = seq(0.5, 2, 0.5), col = "grey")
}

# Modeling plots
# Change data.frame format to be plotted
fox.pcbi.2 <- melt(fox.pcbi, id.vars = c("date"))
# Plot
ggplot(fox.pcbi.2, aes(x = date, y = 10^(value), group = variable)) +
  geom_point(aes(shape = variable, color = variable, size = variable,
                 fill = variable)) +
  scale_shape_manual(values = c(21, 3, 4)) +
  scale_color_manual(values = c('black','#8856a7', '#b2df8a')) +
  scale_size_manual(values = c(2, 1, 1)) +
  scale_fill_manual(values = c("#1b98e0", '#8856a7', '#b2df8a')) +
  scale_x_date(labels = date_format("%Y-%m")) +
  scale_y_log10(limits = c(1, 1000)) +
  xlab("") +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste("PCB 17 (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

