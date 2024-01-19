## Water PCB concentrations analysis.
## Data were obtained from EPA and contractors from PCB Superfund
## sites in USA. Using log10 of the sum of PCB.

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
install.packages("reshape")
install.packages("sf")

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
  library(reshape)
  library(sf)
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# Total Concentration Analysis --------------------------------------------
# Data preparation
{
  # Change date format
  SampleDate <- as.Date(wdc$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Create individual code for each site sampled
  site.numb <- wdc$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(wdc$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  tpcb <- cbind(factor(wdc$SiteID), SampleDate,
                wdc$Latitude, wdc$Longitude, wdc$tPCB,
                data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                      "tPCB", "time", "site.code", "season")
}

# Total PCB Regressions ---------------------------------------------------
# Get variables
tPCB <- tpcb$tPCB
time <- tpcb$time
site <- tpcb$site.code
season <- tpcb$season
# (1) Perform linear regression (lr)
# tPCB vs. time
lr.tpcb.t <- lm(log10(tPCB) ~ time)
# See results
summary(lr.tpcb.t)
# Plot residuals. Create a Q-Q plot and save it.
{
  # Create a new PNG graphics device
  png("Output/Plots/Global/qq_plotlrtPCBV02.png", width = 800, height = 600)
  res <- resid(lr.tpcb.t) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res)
  # Close the PNG device
  dev.off()
}

# Modeling plots
# (1) Get predicted values tpcb
fit.values.lr.tpcb <- as.data.frame(fitted(lr.tpcb.t))
# Add column name
colnames(fit.values.lr.tpcb) <- c("lr.predicted")
# Add predicted values to data.frame
tpcb$lr.predicted <- 10^(fit.values.lr.tpcb$lr.predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(tpcb, aes(x = tPCB, y = lr.predicted)) +
  geom_point(shape = 21, size = 2, fill = "#66ccff") +
  scale_y_log10(limits = c(0.1, 10^8),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.1, 10^8),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lr concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 1) +
  geom_abline(intercept = 0.30103, slope = 1, col = "blue",
              linewidth = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = -0.30103, slope = 1, col = "blue",
              linewidth = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# (2) tPCB vs. season
lr.tpcb.s <- lm(log10(tPCB) ~ season)
# See results
summary(lr.tpcb.s)
# Look at residuals
{
  res <- resid(lr.tpcb.s) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}

# (3) MLR
mlr.tpcb <- lm(log10(tPCB) ~ time + season)
# See results
summary(mlr.tpcb)
# Look at residuals
{
  res <- resid(mlr.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res)
  # Add a straight diagonal line to the plot
  qqline(res)
}

# (4) Perform Linear Mixed-Effects Model (lme)
lme.tpcb <- lmer(log10(tPCB) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.tpcb)
#Create a Q-Q plot and save it.
{
  # Create a new PNG graphics device
  png("Output/Plots/Global/qq_plotlmetPCBV02.png", width = 800, height = 600)
  res <- resid(lme.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res)
  # Close the PNG device
  dev.off()
}
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lme.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lme.tpcb))[1, 'R2c']
# Extract coefficient values
time.coeff <- summary(lme.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lme.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -log(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# Modeling plots
# (1) Get predicted values tpcb
fit.values.tpcb <- as.data.frame(fitted(lme.tpcb))
# Add column name
colnames(fit.values.tpcb) <- c("lme.predicted")
# Add predicted values to data.frame
tpcb$lmepredicted <- 10^(fit.values.tpcb$lme.predicted)

# Plot prediction vs. observations, 1:1 line
tPCBObsPred <- ggplot(tpcb, aes(x = tPCB, y = lmepredicted)) +
  geom_point(shape = 21, size = 2, fill = "white") +
  scale_y_log10(limits = c(0.1, 10^8),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.1, 10^8),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = 0.30103, slope = 1, col = "blue",
              linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = -0.30103, slope = 1, col = "blue",
              linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  annotation_logticks(sides = "bl")

# See plot
print(tPCBObsPred)

# Save plot in folder
ggsave("Output/Plots/Global/tPCBObsPredV04.png", plot = tPCBObsPred,
       width = 8, height = 6, dpi = 300)

# Plot residuals vs. predictions
  {
    # Open a PNG graphics device
    png("Output/Plots/Global/res_plotlmetPCBV02.png", width = 800, height = 600)
    # Create your plot
    plot(tpcb$lmepredicted, resid(lme.tpcb),
         points(tpcb$lmepredicted, resid(lme.tpcb), pch = 16, col = "white"),
         ylim = c(-4, 4),
         xlim = c(1, 10^6.1),
         xlab = expression(paste("Predicted lme concentration ",
                                 Sigma, "PCB (pg/L)")),
         ylab = "Residual")
    # Add lines to the plot
    abline(0, 0)
    abline(h = seq(-4, 4, 1), col = "grey")
    abline(v = seq(0, 1200000, 200000), col = "grey")
    # Close the PNG graphics device
    dev.off()
  }

# Individual PCB Regressions ----------------------------------------------
# PCB5.8
pcb5.8 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB5.8,
               data.frame(time.day), site.numb, season.s)
# Add column names
colnames(pcb5.8) <- c("SiteID", "date", "PCB5.8", "time",
                     "site.code", "season")
# Remove 0s and NA values
pcb5.8 <- pcb5.8[complete.cases(pcb5.8$PCB5.8) & pcb5.8$PCB5.8 != 0, ]

# Get variables
PCBi <- pcb5.8$PCB5.8
time <- pcb5.8$time
site <- pcb5.8$site.code
season <- pcb5.8$season

# (1) Perform linear regression (lr)
# (1.1) PCB vs. time
lr.pcb5.8.t <- lm(log10(PCBi) ~ time)
# See results
summary(lr.pcb5.8.t)

# (1.2) PCB vs. season
lr.pcb5.8.s <- lm(log10(PCBi) ~ season)
# See results
summary(lr.pcb5.8.s)

# (1.3) MLR
mlr.pcb5.8 <- lm(log10(PCBi) ~ time + season)
# See results
summary(mlr.pcb5.8)

# (1.4) Perform Linear Mixed-Effects Model (lme)
lme.pcb5.8 <- lmer(log10(PCBi) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.pcb5.8)

# Shapiro test
shapiro.test(resid(lme.pcb5.8)) # p-value <<< 0.5

# PCB11
pcb11 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB11,
               data.frame(time.day), site.numb, season.s)
# Add column names
colnames(pcb11) <- c("SiteID", "date", "PCB11", "time",
                     "site.code", "season")
# Remove 0s and NA values
pcb11 <- pcb11[complete.cases(pcb11$PCB11) & pcb11$PCB11 != 0, ]

# Get variables
PCBi <- pcb11$PCB11
time <- pcb11$time
site <- pcb11$site.code
season <- pcb11$season

# (1) Perform linear regression (lr)
# (1.1) PCB vs. time
lr.pcb11.t <- lm(log10(PCBi) ~ time)
# See results
summary(lr.pcb11.t)

# (1.2) PCB vs. season
lr.pcb11.s <- lm(log10(PCBi) ~ season)
# See results
summary(lr.pcb11.s)

# (1.3) MLR
mlr.pcb11 <- lm(log10(PCBi) ~ time + season)
# See results
summary(mlr.pcb11)

# (1.4) Perform Linear Mixed-Effects Model (lme)
lme.pcb11 <- lmer(log10(PCBi) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.pcb11)

# Shapiro test
shapiro.test(resid(lme.pcb11)) # p-value <<< 0.5

# PCB18.30
pcb18.30 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB18.30,
               data.frame(time.day), site.numb, season.s)
# Add column names
colnames(pcb18.30) <- c("SiteID", "date", "PCB18.30", "time",
                     "site.code", "season")
# Remove 0s and NA values
pcb18.30 <- pcb18.30[complete.cases(pcb18.30$PCB18.30) & pcb18.30$PCB18.30 != 0, ]

# Get variables
PCBi <- pcb18.30$PCB18.30
time <- pcb18.30$time
site <- pcb18.30$site.code
season <- pcb18.30$season

# (1) Perform linear regression (lr)
# (1.1) PCB vs. time
lr.pcb18.30.t <- lm(log10(PCBi) ~ time)
# See results
summary(lr.pcb18.30.t)

# (1.2) PCB vs. season
lr.pcb18.30.s <- lm(log10(PCBi) ~ season)
# See results
summary(lr.pcb18.30.s)

# (1.3) MLR
mlr.pcb18.30 <- lm(log10(PCBi) ~ time + season)
# See results
summary(mlr.pcb18.30)

# (1.4) Perform Linear Mixed-Effects Model (lme)
lme.pcb18.30 <- lmer(log10(PCBi) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.pcb18.30)

# Shapiro test
shapiro.test(resid(lme.pcb18.30)) # p-value <<< 0.5

# PCB20.21.28.31.33.50.53
pcb20 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB20.21.28.31.33.50.53,
               data.frame(time.day), site.numb, season.s)
# Add column names
colnames(pcb20) <- c("SiteID", "date", "PCB20", "time",
                     "site.code", "season")
# Remove 0s and NA values
pcb20 <- pcb20[complete.cases(pcb20$PCB20) & pcb20$PCB20 != 0, ]

# Get variables
PCBi <- pcb20$PCB20
time <- pcb20$time
site <- pcb20$site.code
season <- pcb20$season

# (1) Perform linear regression (lr)
# (1.1) PCB vs. time
lr.pcb20.t <- lm(log10(PCBi) ~ time)
# See results
summary(lr.pcb20.t)

# (1.2) PCB vs. season
lr.pcb20.s <- lm(log10(PCBi) ~ season)
# See results
summary(lr.pcb11.s)

# (1.3) MLR
mlr.pcb20 <- lm(log10(PCBi) ~ time + season)
# See results
summary(mlr.pcb20)

# (1.4) Perform Linear Mixed-Effects Model (lme)
lme.pcb20 <- lmer(log10(PCBi) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.pcb20)

# Shapiro test
shapiro.test(resid(lme.pcb20)) # p-value <<< 0.5

# PCB44+47+65
pcb44 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB44.47.65,
               data.frame(time.day), site.numb, season.s)
# Add column names
colnames(pcb44) <- c("SiteID", "date", "PCB44", "time",
                     "site.code", "season")
# Remove 0s and NA values
pcb44 <- pcb44[complete.cases(pcb44$PCB44) & pcb44$PCB44 != 0, ]

# Get variables
PCBi <- pcb44$PCB44
time <- pcb44$time
site <- pcb44$site.code
season <- pcb44$season

# (1) Perform linear regression (lr)
# (1.1) PCB vs. time
lr.pcb44.t <- lm(log10(PCBi) ~ time)
# See results
summary(lr.pcb44.t)

# (1.2) PCB vs. season
lr.pcb44.s <- lm(log10(PCBi) ~ season)
# See results
summary(lr.pcb44.s)

# (1.3) MLR
mlr.pcb44 <- lm(log10(PCBi) ~ time + season)
# See results
summary(mlr.pcb44)

# (1.4) Perform Linear Mixed-Effects Model (lme)
lme.pcb44 <- lmer(log10(PCBi) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.pcb44)

# Shapiro test
shapiro.test(resid(lme.pcb44)) # p-value <<< 0.5

# PCB 67
pcb67 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB67,
               data.frame(time.day), site.numb, season.s)
# Add column names
colnames(pcb67) <- c("SiteID", "date", "PCB67", "time",
                     "site.code", "season")
# Remove 0s and NA values
pcb67 <- pcb67[complete.cases(pcb67$PCB67) & pcb67$PCB67 != 0, ]

# Get variables
PCBi <- pcb67$PCB67
time <- pcb67$time
site <- pcb67$site.code
season <- pcb67$season

# (1) Perform linear regression (lr)
# (1.1) PCB vs. time
lr.pcb67.t <- lm(log10(PCBi) ~ time)
# See results
summary(lr.pcb67.t)

# (1.2) PCB vs. season
lr.pcb67.s <- lm(log10(PCBi) ~ season)
# See results
summary(lr.pcb67.s)

# (1.3) MLR
mlr.pcb67 <- lm(log10(PCBi) ~ time + season)
# See results
summary(mlr.pcb67)

# (1.4) Perform Linear Mixed-Effects Model (lme)
lme.pcb67 <- lmer(log10(PCBi) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.pcb67)

# Shapiro test
shapiro.test(resid(lme.pcb67)) # p-value <<< 0.5

# PCB 106+118
pcb106.118 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB106.118,
                    data.frame(time.day), site.numb, season.s)
# Add column names
colnames(pcb106.118) <- c("SiteID", "date", "PCB106.118", "time",
                          "site.code", "season")
# Remove 0s and NA values
pcb106.118 <- pcb106.118[complete.cases(pcb106.118$PCB106.118) & pcb106.118$PCB106.118 != 0, ]

# Get variables
PCBi <- pcb106.118$PCB106.118
time <- pcb106.118$time
site <- pcb106.118$site.code
season <- pcb106.118$season

# (1) Perform linear regression (lr)
# (1.1) PCB vs. time
lr.pcb106.118.t <- lm(log10(PCBi) ~ time)
# See results
summary(lr.pcb106.118.t)

# (1.2) PCB vs. season
lr.pcb106.118.s <- lm(log10(PCBi) ~ season)
# See results
summary(lr.pcb106.118.s)

# (1.3) MLR
mlr.pcb106.118 <- lm(log10(PCBi) ~ time + season)
# See results
summary(mlr.pcb106.118)

# (1.4) Perform Linear Mixed-Effects Model (lme)
lme.pcb106.118 <- lmer(log10(PCBi) ~ 1 + time + season + (1|site),
                       REML = FALSE,
                       control = lmerControl(check.nobs.vs.nlev = "ignore",
                                             check.nobs.vs.rankZ = "ignore",
                                             check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.pcb106.118)

# Shapiro test
shapiro.test(resid(lme.pcb106.118)) # p-value <<< 0.5
