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
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# General information -----------------------------------------------------
# Number of locations and number of site per location
location_count <- wdc %>%
  group_by(LocationName) %>%
  summarise(count = n())

print(location_count)
# Median amount of samples per location
median(location_count$count)

# Number of locations per states
state_count <- wdc %>%
  group_by(StateSampled) %>%
  summarise(count = n())

print(state_count)

# Number of sites and number of replicates per site
site_count <- wdc %>%
  group_by(SiteID) %>%
  summarise(count = n())

print(site_count)

# Media of number of site available
median(site_count$count)

# Number of replicates per sites from the same day
site_repli_count <- wdc %>%
  group_by(SampleID) %>%
  summarise(count = n())

print(site_repli_count)

# Media of number of replicates per sites from the same day
median(site_repli_count$count)

# Find the SampleID with the highest count
max_count_sample <- site_repli_count %>%
  filter(count == max(count)) %>%
  pull(SampleID)

filtered_wdc <- wdc %>%
  filter(SampleID == max_count_sample)

# Extract Site Name
extracted_string <- filtered_wdc %>%
  select(SiteName) %>%
  pull()

# Display Site Name
print(extracted_string[1])

# Aroclor summary ---------------------------------------------------------
# Number of samples analyzed using Aroclor method
count_Aroclor <- sum(wdc$AroclorCongener == "Aroclor")
total_samples <- length(wdc[,1])
percent_aroclor <- count_Aroclor/total_samples*100

# Calculate sample % for each Aroclor mixtures
aroclors <- c('A1016', 'A1221', 'A1232', 'A1242', 'A1248',
              'A1254', 'A1260')

# Calculate the number of non-NA values (frequency of numbers) in each Aroclor
frequency_aroclors <- lapply(wdc[aroclors], function(column) {
  length(na.omit(column))
})

# Percentage of Aroclor mixtures in relation to all Aroclors
for (i in seq_along(aroclors)) {
  column_name <- aroclors[i]
  print(paste(column_name))
  print(frequency_aroclors[[i]]/count_Aroclor*100)
}

# Congener frequency ------------------------------------------------------
{
  # Remove metadata
  wdc.nmd <- subset(wdc, select = -c(SampleID:AroclorCongener)) #nmd = no meta data
  # Remove Aroclor data
  wdc.nmd <- subset(wdc.nmd, select = -c(A1016:A1260))
  # (2) Only consider congener data
  wdc.cong <- subset(wdc, AroclorCongener == "Congener")
  # Remove metadata
  wdc.cong.1 <- subset(wdc.cong, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data and tPCB
  wdc.cong.1 <- subset(wdc.cong.1, select = -c(A1016:tPCB))
}

# Create a frequency detection plot
{
  wdc.cong.freq <- colSums(! is.na(wdc.cong.1) & (wdc.cong.1 !=0))/nrow(wdc.cong.1)
  wdc.cong.freq <- data.frame(wdc.cong.freq)
  colnames(wdc.cong.freq) <- c("PCB.frequency")
  congener <- row.names(wdc.cong.freq)
  wdc.cong.freq <- cbind(congener, wdc.cong.freq$PCB.frequency)
  colnames(wdc.cong.freq) <- c("congener", "PCB.frequency")
  wdc.cong.freq <- data.frame(wdc.cong.freq)
  wdc.cong.freq$congener <- as.character(wdc.cong.freq$congener)
  wdc.cong.freq$congener <- gsub('\\.', '+', wdc.cong.freq$congener) # replace dot for +
  wdc.cong.freq$PCB.frequency <- as.numeric(as.character(wdc.cong.freq$PCB.frequency))
  wdc.cong.freq$congener <- factor(wdc.cong.freq$congener,
                                   levels = rev(wdc.cong.freq$congener)) # change the order to be plotted.
}

# Summary statistic of frequency of detection
summary(wdc.cong.freq$PCB.frequency)

# Frequency detection plot
plot.cong.freq <- ggplot(wdc.cong.freq, aes(x = 100*PCB.frequency, y = congener)) +
  geom_bar(stat = "identity", fill = "#66ccff", color = "black") +
  geom_vline(xintercept = 100*mean(wdc.cong.freq$PCB.frequency),
             color = "red") +
  ylab("") +
  theme_bw() +
  xlim(c(0,100)) +
  theme(aspect.ratio = 20/5) +
  xlab(expression(bold("Frequency detection (%)"))) +
  theme(axis.text.x = element_text(face = "bold", size = 8),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.text.y = element_text(face = "bold", size = 7))

print(plot.cong.freq)  # Print the plot

# Save map in folder
ggsave("Output/Plots/Global/FreqPCBV01.png", plot = plot.cong.freq,
       width = 5, height = 10, dpi = 300)

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

# Get coordinates per site to plot in Google Earth
location <- tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = location, mean)

# Summary statistic of total PCB (congeners + Aroclor) in pg/L not including 0s
summary(wdc$tPCB)

# Highest sample location and time
max_sample <- wdc[which.max(wdc$tPCB), ]
max_sample <- max_sample[c("LocationName", "SampleDate", "SiteName")]
print(max_sample)

# Global plots ------------------------------------------------------------
# Histogram
hist(tpcb$tPCB)
hist(log10(tpcb$tPCB))

## Total PCBs in 1 box plot
## include 64 and 640 pg/L from EPA
plot.box.tPCB <- ggplot(tpcb, aes(x = "", y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "white") +
  geom_boxplot(lwd = 0.8, width = 0.7, outlier.shape = NA, alpha = 0) +
  theme_classic() +
  theme(aspect.ratio = 14/2) +
  xlab(expression(bold(Sigma*"PCB")))+
  ylab(expression(bold("Water Concentration 1979 - 2020 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 16),
        axis.title.y = element_text(face = "bold", size = 16)) +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, vjust = 5)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_hline(yintercept = 640, color = "#9999CC",
             linewidth = 0.8) +
  geom_hline(yintercept = 64, color = "#CC6666",
             linewidth = 0.8)

# Print or save the plot
print(plot.box.tPCB)

# Save map in folder
ggsave("Output/Plots/Global/tPCBBoxPlotV02.png", plot = plot.box.tPCB,
       width = 5, height = 10, dpi = 300)

# Calculate % samples above both EPA thresholds
EPA640 <- sum(tpcb$tPCB > 640)/nrow(tpcb)*100
EPA64 <- sum(tpcb$tPCB > 64)/nrow(tpcb)*100

# Individual congeners
summary(wdc.cong.1, na.rm = T, zero = T)
# Get the max value for each congener
cong.max <-as.numeric(sub('.*:', '',
                          summary(wdc.cong.1, na.rm = T,
                                  zero = T)[6,]))
# Add congener
cong.max <- cbind(congener, data.frame(cong.max))

# Obtain the median for each individual congener
cong.median <- as.numeric(sub('.*:',
                              '', summary(wdc.cong.1, na.rm = T,
                                          zero = T)[3,]))
# Add congener
cong.median <- cbind(congener, data.frame(cong.median))

# Min
print(min(cong.median$cong.median))
#Max
print(max(cong.median$cong.median))
# Mean
print(mean(cong.median$cong.median))

# Individual PCB boxplot
PCBi_boxplot <- ggplot(stack(wdc.cong.1), aes(x = ind, y = values)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot(width = 0.6, shape = 21, outlier.fill = "white",
               fill = "white", outlier.shape = 21) +
  scale_x_discrete(labels = wdc.cong.freq$congener) + # use to change the "." to "+"
  theme_bw() +
  theme(aspect.ratio = 25/135) +
  xlab(expression("")) +
  ylab(expression(bold("PCB congener concentration (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8,
                                   color = "black"),
        axis.title.y = element_text(face = "bold", size = 8,
                                    color = "black")) +
  theme(axis.text.x = element_text(face = "bold", size = 6,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.6, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# Print the plot
print(PCBi_boxplot)

# Save map in folder
ggsave("Output/Plots/Global/PCBiBoxPlotV02.png", plot = PCBi_boxplot,
       width = 10, height = 5, dpi = 300)

# Regression analysis and plots---------------------------------------------
# Plots
# (1) Time trend plots
plot.time.tPCB <- ggplot(tpcb, aes(y = tPCB,
                                   x = format(date,'%Y'))) +
  geom_point(shape = 21, size = 2.5, fill = "white") +
  xlab("") +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB 1979 - 2020 (pg/L)"))))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  annotation_logticks(sides = "l") +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black")) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 11))

# Print the plot
print(plot.time.tPCB)

# Save plot in folder
ggsave("Output/Plots/Global/tPCBTimeV02.png", plot = plot.time.tPCB,
       width = 10, height = 3, dpi = 300)

# (2) Seasonality
ggplot(tpcb, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB 1979 - 2020 (pg/L)"))))) +
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
  png("Output/Plots/Global/qq_plotlrtPCBV01.png", width = 800, height = 600)
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
  geom_abline(intercept = 0, slope = 1, col = "black", size = 1) +
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
  png("Output/Plots/Global/qq_plotlmetPCBV01.png", width = 800, height = 600)
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

print(tPCBObsPred)  # Print the plot

# Save plot in folder
ggsave("Output/Plots/Global/tPCBObsPredV03.png", plot = tPCBObsPred,
       width = 8, height = 6, dpi = 300)

# Plot residuals vs. predictions
  {
    # Open a PNG graphics device
    png("Output/Plots/Global/res_plotlmetPCBV01.png", width = 800, height = 600)
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
# States
sites <- c("CA", "CT", "DC*", "DE", "ID", "IN", "MA", "MD", "MI", "MO",
           "MT", "NM", "NY", "OH", "ON*", "OR", "PA", "TX", "VA", "WA", "WI")

# Total PCBs
tpcb.state <- ggplot(wdc, aes(x = factor(StateSampled, levels = sites),
                y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB 1979 - 2019 (pg/L)"))))) +
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
  geom_boxplot(lwd = 0.5, width = 0.7, outlier.shape = NA, alpha = 0) +
  geom_hline(yintercept = 640, color = "#9999CC",
             size = 0.8) + # U.S. EPA Water Quality Criterion for Human Health from fish consumption, associated with an incremental cancer risk of 10−5
  geom_hline(yintercept = 64, color = "#CC6666",
             size = 0.8) # associated with an incremental cancer risk of 10−6.

print(tpcb.state)

# Save map in folder
ggsave("Output/Plots/Global/tPCBStateV01.png", plot = tpcb.state,
       width = 10, height = 5, dpi = 300)

# Selected StateSampled and individual PCB congener
wdc.pcbi <- subset(wdc, select = c(StateSampled, PCB5.8))
# Remove samples with 0s
wdc.pcbi <- wdc.pcbi[!(wdc.pcbi[2] == 0), ]
# Remove samples this NA
wdc.pcbi <- wdc.pcbi[!is.na(wdc.pcbi[2]),]
# Plot
ggplot(wdc.pcbi, aes(x = factor(StateSampled, levels = sites),
                   y = PCB5.8)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concentration",
                            paste("PCB 4+10 1979 - 2019 (pg/L)"))))) +
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
  geom_boxplot(lwd = 0.5, width = 0.7, outlier.shape = NA, alpha = 0)


