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
install.packages("dataRetrieval")
install.packages("reshape")
install.packages("sf")
install.packages("sfheaders")

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
  library(sf)
  library(sfheaders) # Create file to be used in Google Earth
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

# Regressions -------------------------------------------------------------
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
lmem.tpcb <- lmer(log10(tPCB) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lmem.tpcb)
#Create a Q-Q plot and save it.
{
  # Create a new PNG graphics device
  png("Output/Plots/Global/qq_plotlmetPCBV02.png", width = 800, height = 600)
  res <- resid(lmem.tpcb) # get list of residuals
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
R2.nre <- as.data.frame(r.squaredGLMM(lmem.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.tpcb))[1, 'R2c']
# Extract coefficient values
time.coeff <- summary(lmem.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lmem.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -log(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# Modeling plots
# (1) Get predicted values tpcb
fit.values.tpcb <- as.data.frame(fitted(lmem.tpcb))
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
    plot(tpcb$lmepredicted, resid(lmem.tpcb),
         points(tpcb$lmepredicted, resid(lmem.tpcb), pch = 16, col = "white"),
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

#Until here!

# Spatial Plots and Analysis ----------------------------------------------
# tPCB
# List of sites you want to include in the plot
sites_to_include <- c("Housatonic River", "New Bedford Harbor", "Passaic River",
                      "Hudson River", "Kalamazoo River", "Fox River",
                      "Portland Harbor", "Lake Michigan Mass Balance",
                      "Spokane River", "Chesapeake Bay")

# Filter the data to include only the specified sites
filtered_data <- wdc %>%
  filter(LocationName %in% sites_to_include)

# Total PCBs
tpcb.site <- ggplot(filtered_data, aes(x = factor(LocationName),
                              y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 20/15) +
  ylab(expression(bold("Water Concentration " *Sigma*"PCB (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 12,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "white") +
  geom_boxplot(lwd = 0.5, width = 0.7, outlier.shape = NA, alpha = 0)

print(tpcb.site)

# Save plot in folder
ggsave("Output/Plots/Global/tPCBSiteV02.png", plot = tpcb.site,
       width = 5, height = 10, dpi = 300)

# Individual congeners
# Filter out rows with NA and 0 values in the 'PCBi' column
filtered_datai <- wdc %>%
  filter(LocationName %in% sites_to_include, !is.na(PCB11),
         !(PCB11 == 0))

# Create the ggplot
pcbi.site <- ggplot(filtered_datai, aes(x = factor(LocationName),
                                        y = PCB11)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 20/15) +
  ylab(expression(bold("PCB 11 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 14)) +
  theme(axis.text.x = element_text(face = "bold", size = 12,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "white") +
  geom_boxplot(lwd = 0.5, width = 0.7, outlier.shape = NA, alpha = 0)

print(pcbi.site)

# Save plot in folder
ggsave("Output/Plots/Global/PCB11Site.png", plot = pcbi.site,
       width = 5, height = 10, dpi = 300)

filtered_datai <- wdc %>%
  filter(LocationName %in% sites_to_include, !is.na(PCB20.21.28.31.33.50.53),
         !(PCB20.21.28.31.33.50.53 == 0))

# Create the ggplot
pcbi.site <- ggplot(filtered_datai, aes(x = factor(LocationName),
                                        y = PCB20.21.28.31.33.50.53)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 20/15) +
  ylab(expression(bold("PCB 20+21+28+31+33+50+53 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 14,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "white") +
  geom_boxplot(lwd = 0.5, width = 0.7, outlier.shape = NA, alpha = 0)

print(pcbi.site)

# Save plot in folder
ggsave("Output/Plots/Global/PCB20Site.png", plot = pcbi.site,
       width = 5, height = 10, dpi = 300)

# PCB44+47+65
filtered_datai <- wdc %>%
  filter(LocationName %in% sites_to_include, !is.na(PCB44.47.65),
         !(PCB44.47.65 == 0))

# Create the ggplot
pcbi.site <- ggplot(filtered_datai, aes(x = factor(LocationName),
                                        y = PCB44.47.65)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 20/15) +
  ylab(expression(bold("PCB 44+47+65 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 14,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "white") +
  geom_boxplot(lwd = 0.5, width = 0.7, outlier.shape = NA, alpha = 0)

print(pcbi.site)

# Save plot in folder
ggsave("Output/Plots/Global/PCB44Site.png", plot = pcbi.site,
       width = 5, height = 10, dpi = 300)

# PCB 67
filtered_datai <- wdc %>%
  filter(LocationName %in% sites_to_include, !is.na(PCB67),
         !(PCB67 == 0))

# Create the ggplot
pcbi.site <- ggplot(filtered_datai, aes(x = factor(LocationName),
                                        y = PCB67)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 20/15) +
  ylab(expression(bold("PCB 67 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 14,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "white") +
  geom_boxplot(lwd = 0.5, width = 0.7, outlier.shape = NA, alpha = 0)

print(pcbi.site)

# Save plot in folder
ggsave("Output/Plots/Global/PCB67Site.png", plot = pcbi.site,
       width = 5, height = 10, dpi = 300)

