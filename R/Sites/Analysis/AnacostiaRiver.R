## Water PCB concentrations data analysis per site
## Anacostia River
## data only tPCB

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
install.packages("tidyr")
install.packages('patchwork')
install.packages("scales")
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
  library(tidyr) # function gather
  library(patchwork) # combine plots
  library(sf) # Create file to be used in Google Earth
  library(sfheaders) # Create file to be used in Google Earth
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# Select anrsatonic River data ---------------------------------------------------
anr <- wdc[str_detect(wdc$LocationName, 'Anacostia River'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  anr$SampleDate <- as.Date(anr$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(anr$SampleDate) - min(as.Date(anr$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- anr$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(anr$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  anr.tpcb <- cbind(factor(anr$SiteID), anr$SampleDate,
                    anr$Latitude, anr$Longitude, as.matrix(anr$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(anr.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Get coordinates per site to plot in Google Earth
location <- anr.tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = location, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/AnacostiaRiverLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# General plots -------------------------------------------------------------------
# (1) Histograms
hist(anr.tpcb$tPCB)
hist(log10(anr.tpcb$tPCB))

# (2) Time trend plots
ANRTime <- ggplot(anr.tpcb, aes(y = tPCB, x = format(date, '%Y-%m'))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  xlab("") +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000, 10000),  # Specify the desired breaks
    labels = label_comma()(c(1, 10, 100, 1000, 10000))  # Specify the desired labels
  ) +
  theme_classic() +
  ylab(expression(bold(Sigma*"PCB (pg/L)"))) +
  theme(
    axis.text.y = element_text(face = "bold", size = 22),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.x = element_text(size = 24, angle = 60, hjust = 1),
    axis.title.x = element_text(face = "bold", size = 17),
    plot.margin = margin(0, 0, 0, 0, unit = "cm"))

# Pint plot
print(ANRTime)

# Save plot in folder
ggsave("Output/Plots/Sites/Temporal/AnacostiaRiverTime.png",
       plot = ANRTime, width = 18, height = 8, dpi = 500)

# (3) Seasonality
ggplot(anr.tpcb, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotate("text", x = 4, y = 10^5.6, label = "anrsotonic River",
           size = 3)

# (4) Sites
ggplot(anr.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concentration",
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
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotate("text", x = 15, y = 10^4.5, label = "Anacostia River",
           size = 3)

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- anr.tpcb$tPCB
time <- anr.tpcb$time
site <- anr.tpcb$site.code
season <- anr.tpcb$season

# tPCB vs. time + flow + season + site
lme.anr.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lme.anr.tpcb)
# Look at residuals
{
  res.anr.tpcb <- resid(lme.anr.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/AnacostiaRiverQ-QtPCB.pdf")
  qqnorm(res.anr.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.anr.tpcb)
  dev.off()
}
# Shapiro test
shapiro.test(resid(lme.anr.tpcb)) # p-value = 0.04

# Create matrix to store results
{
  lme.tpcb <- matrix(nrow = 1, ncol = 21)
  lme.tpcb[1] <- fixef(lme.anr.tpcb)[1] # intercept
  lme.tpcb[2] <- summary(lme.anr.tpcb)$coef[1,"Std. Error"] # intercept error
  lme.tpcb[3] <- summary(lme.anr.tpcb)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb[4] <- fixef(lme.anr.tpcb)[2] # time
  lme.tpcb[5] <- summary(lme.anr.tpcb)$coef[2,"Std. Error"] # time error
  lme.tpcb[6] <- summary(lme.anr.tpcb)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb[7] <- fixef(lme.anr.tpcb)[3] # season 1
  lme.tpcb[8] <- summary(lme.anr.tpcb)$coef[3,"Std. Error"] # season 1 error
  lme.tpcb[9] <- summary(lme.anr.tpcb)$coef[3,"Pr(>|t|)"] # season 1 p-value
  lme.tpcb[10] <- fixef(lme.anr.tpcb)[4] # season 2
  lme.tpcb[11] <- summary(lme.anr.tpcb)$coef[4,"Std. Error"] # season 2 error
  lme.tpcb[12] <- summary(lme.anr.tpcb)$coef[4,"Pr(>|t|)"] # season 2 p-value
  lme.tpcb[13] <- fixef(lme.anr.tpcb)[5] # season 3
  lme.tpcb[14] <- summary(lme.anr.tpcb)$coef[5,"Std. Error"] # season 3 error
  lme.tpcb[15] <- summary(lme.anr.tpcb)$coef[5,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb[16] <- -log(2)/lme.tpcb[4]/365 # t0.5
  lme.tpcb[17] <- abs(-log(2)/lme.tpcb[4]/365)*lme.tpcb[5]/abs(lme.tpcb[4]) # t0.5 error
  lme.tpcb[18] <- as.data.frame(VarCorr(lme.anr.tpcb))[1,'sdcor']
  lme.tpcb[19] <- as.data.frame(r.squaredGLMM(lme.anr.tpcb))[1, 'R2m']
  lme.tpcb[20] <- as.data.frame(r.squaredGLMM(lme.anr.tpcb))[1, 'R2c']
  lme.tpcb[21] <- shapiro.test(resid(lme.anr.tpcb))$p.value
}

# Just 3 significant figures
lme.tpcb <- formatC(signif(lme.tpcb, digits = 3))
# Add column names
colnames(lme.tpcb) <- c("Intercept", "Intercept.error",
                        "Intercept.pv", "time", "time.error", "time.pv",
                        "season1", "season1.error", "season1.pv",
                        "season2", "season2.error", "season2.pv",
                        "season3", "season3.error", "season3.pv", "t05",
                        "t05.error", "RandonEffectSiteStdDev", "R2nR", "R2R",
                        "Normality")

# Export results
write.csv(lme.tpcb,
          file = "Output/Data/Sites/csv/AnacostiaRiver/AnacostiaRiverLmetPCB.csv")

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.anr.tpcb <- as.data.frame(fitted(lme.anr.tpcb))
# Add column name
colnames(fit.lme.values.anr.tpcb) <- c("predicted")
# Add predicted values to data.frame
anr.tpcb$predicted <- 10^(fit.lme.values.anr.tpcb$predicted)
# Create overall plot prediction vs. observations
predic.obs <- data.frame(tPCB = anr.tpcb$tPCB, predicted = anr.tpcb$predicted)
predic.obs <- data.frame(Location = anr$LocationName[1], predic.obs)
# Save new data
write.csv(predic.obs,
          "Output/Data/Sites/csv/AnacostiaRiver/AnacostiaRiverObsPredtPCB.csv")

# Plot prediction vs. observations, 1:1 line
p <- ggplot(anr.tpcb, aes(x = tPCB, y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(1, 10^5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# See plot
print(p)

# Save plot
ggsave("Output/Plots/Sites/ObsPred/AnacostiaRiver/AnacostiaRiverObsPredtPCB.png",
       plot = p, width = 8, height = 8, dpi = 500)

# Plot residuals vs. predictions
{
  # Open a PNG graphics device
  png("Output/Plots/Sites/Residual/res_plotlmeAnacostiaRivertPCB.png",
      width = 800, height = 600)
  # Create your plot
  plot(anr.tpcb$predicted, resid(lme.anr.tpcb),
       points(anr.tpcb$predicted, resid(lme.anr.tpcb), pch = 16, 
              col = "white"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  # Add lines to the plot
  abline(0, 0)
  abline(h = seq(-2, 2, 1), col = "grey")
  abline(v = seq(0, 10000, 2000), col = "grey")
  # Close the PNG graphics device
  dev.off()
}

# Estimate a factor of 2 between observations and predictions
anr.tpcb$factor2 <- anr.tpcb$tPCB/anr.tpcb$predicted
factor2.tpcb <- nrow(anr.tpcb[anr.tpcb$factor2 > 0.5 & anr.tpcb$factor2 < 2,
                              ])/length(anr.tpcb[,1])*100

# Convert the vector to a data frame
factor2.tpcb <- data.frame(Factor_2 = factor2.tpcb)

# Export results
write.csv(factor2.tpcb,
          file = "Output/Data/Sites/csv/AnacostiaRiver/AnacostiaRiverFactor2tPCB.csv")

