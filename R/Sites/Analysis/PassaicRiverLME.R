## Water PCB concentrations data analysis per site
## Passaic River

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

# Select Passaic River data ---------------------------------------------------
pass <- wdc[str_detect(wdc$LocationName, 'Passaic River'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  pass$SampleDate <- as.Date(pass$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(pass$SampleDate) - min(as.Date(pass$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- pass$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(pass$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  pass.tpcb <- cbind(factor(pass$SiteID), pass$SampleDate,
                    pass$Latitude, pass$Longitude, as.matrix(pass$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(pass.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Get coordinates per site to plot in Google Earth
location <- pass.tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = location, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/PassaicRiverLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# General plots -------------------------------------------------------------------
# (1) Histograms
hist(pass.tpcb$tPCB)
hist(log10(pass.tpcb$tPCB))

# (2) Time trend plots
PASTime <- ggplot(pass.tpcb, aes(y = tPCB, x = format(date, '%Y'))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  xlab("") +
  scale_y_log10(
    breaks = c(10, 100, 1000, 10000, 100000),  # Specify the desired breaks
    labels = label_comma()(c(10, 100, 1000, 10000, 100000))  # Specify the desired labels
  ) +
  theme_classic() +
  ylab(expression(bold(Sigma*"PCB (pg/L)"))) +
  theme(
    axis.text.y = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 20, angle = 60, hjust = 1),
    axis.title.x = element_text(face = "bold", size = 17),
    plot.margin = margin(0.1, 0, 0, 0, unit = "cm"))

# Print plot
print(PASTime)

# Save plot in folder
ggsave("Output/Plots/Sites/Temporal/PassaicRiverTime.png",
       plot = PASTime, width = 6, height = 5, dpi = 500)

# (3) Seasonality
ggplot(pass.tpcb, aes(x = season, y = tPCB)) +
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
              shape = 21, fill = "white") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotate("text", x = 1, y = 20, label = "Passaic River",
           size = 3)

# (4) Sites
ggplot(pass.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
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
              shape = 21, fill = "white") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotate("text", x = 5, y = 20, label = "Passaic River",
           size = 3)

# Include USGS flow and temperature data --------------------------------------------------
{
  # Include flow data from USGS station Passaic River
  sitepassN1 <- "01381900" # No temp
  sitepassN2 <- "01379500" # No temp
  sitepassN3 <- "01389005" # No flow
  sitepassN4 <- "01389010" # No temp
  sitepassN5 <- "01389500" # No temp
  sitepassN6 <- "01389890" # No temp

  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow.1 <- readNWISdv(sitepassN1, paramflow,
                     min(pass.tpcb$date), max(pass.tpcb$date))
  flow.2 <- readNWISdv(sitepassN2, paramflow,
                     min(pass.tpcb$date), max(pass.tpcb$date))
  flow.3 <- readNWISdv(sitepassN4, paramflow,
                     min(pass.tpcb$date), max(pass.tpcb$date))
  flow.4 <- readNWISdv(sitepassN5, paramflow,
                     min(pass.tpcb$date), max(pass.tpcb$date))
  flow.5 <- readNWISdv(sitepassN6, paramflow,
                     min(pass.tpcb$date), max(pass.tpcb$date))
  temp <- readNWISdv(sitepassN3, paramtemp,
                     min(pass.tpcb$date), max(pass.tpcb$date))
  
  # Add USGS data to pass.tpcb.2, matching dates, conversion to m3/s
  pass.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(pass.tpcb$date,
                                                      flow.1$Date)]
  pass.tpcb$flow.2 <- 0.03*flow.1$X_00060_00003[match(pass.tpcb$date,
                                                      flow.2$Date)]
  pass.tpcb$flow.3 <- 0.03*flow.1$X_00060_00003[match(pass.tpcb$date,
                                                      flow.3$Date)]
  pass.tpcb$flow.4 <- 0.03*flow.1$X_00060_00003[match(pass.tpcb$date,
                                                      flow.4$Date)]
  pass.tpcb$flow.5 <- 0.03*flow.1$X_00060_00003[match(pass.tpcb$date,
                                                      flow.5$Date)]
  pass.tpcb$temp <- 273.15 + temp$X_.from.middle.intake_00010_00003[match(pass.tpcb$date,
                                                       temp$Date)]
}

# Remove site -------------------------------------------------------------
# Remove site located in the ocean.Possible typo in original coordinates.
pass.tpcb.1 <- subset(pass.tpcb, SiteID != c("WCPCB-PASS022"))

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- pass.tpcb$tPCB
time <- pass.tpcb$time
site <- pass.tpcb$site.code
season <- pass.tpcb$season
flow <- pass.tpcb$flow.5
tem <- pass.tpcb$temp
# tPCB vs. time + season + flow + temp + site
lme.pass.tpcb <- lmer(log10(tpcb) ~ 1 + time + flow + tem + season + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"))

# See results
summary(lme.pass.tpcb)
# Look at residuals
{
  res.pass.tpcb <- resid(lme.pass.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/PassaicRiverQ-QtPCB.pdf")
  qqnorm(res.pass.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.pass.tpcb)
  dev.off()
}
# Shapiro test
shapiro.test(resid(lme.pass.tpcb)) # Lme doesn't work.

# Individual PCB Analysis -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  pass.pcb <- subset(pass, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  pass.pcb <- subset(pass.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  pass.pcb <- log10(pass.pcb)
  # Replace -inf to NA
  pass.pcb <- do.call(data.frame,
                     lapply(pass.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  pass.pcb.1 <- pass.pcb[,
                       -which(colSums(is.na(pass.pcb))/nrow(pass.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(pass$SiteID)
  # Change date format
  SampleDate <- as.Date(pass$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Create individual code for each site sampled
  site.numb <- pass$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(pass$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to pass.pcb
  pass.pcb.1 <- cbind(pass.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     site.numb, season.s)
  # Include flow data from USGS station Passaic River
  sitepassN1 <- "01381900" # No temp
  sitepassN2 <- "01379500" # No temp
  sitepassN3 <- "01389005" # No flow
  sitepassN4 <- "01389010" # No temp
  sitepassN5 <- "01389500" # No temp
  sitepassN6 <- "01389890" # No temp
  
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow.1 <- readNWISdv(sitepassN1, paramflow,
                       min(pass.pcb.1$SampleDate), max(pass.pcb.1$SampleDate))
  flow.2 <- readNWISdv(sitepassN2, paramflow,
                       min(pass.pcb.1$SampleDate), max(pass.pcb.1$SampleDate))
  flow.3 <- readNWISdv(sitepassN4, paramflow,
                       min(pass.pcb.1$SampleDate), max(pass.pcb.1$SampleDate))
  flow.4 <- readNWISdv(sitepassN5, paramflow,
                       min(pass.pcb.1$SampleDate), max(pass.pcb.1$SampleDate))
  flow.5 <- readNWISdv(sitepassN6, paramflow,
                       min(pass.pcb.1$SampleDate), max(pass.pcb.1$SampleDate))
  temp <- readNWISdv(sitepassN3, paramtemp,
                     min(pass.pcb.1$SampleDate), max(pass.pcb.1$SampleDate))
  
  # Add USGS data to pass.tpcb.1, matching dates, conversion to m3/s
  pass.pcb.1$flow.1 <- 0.03*flow.1$X_00060_00003[match(pass.pcb.1$SampleDate,
                                                      flow.1$Date)]
  pass.pcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(pass.pcb.1$SampleDate,
                                                      flow.2$Date)]
  pass.pcb.1$flow.3 <- 0.03*flow.3$X_00060_00003[match(pass.pcb.1$SampleDate,
                                                      flow.3$Date)]
  pass.pcb.1$flow.4 <- 0.03*flow.4$X_00060_00003[match(pass.pcb.1$SampleDate,
                                                      flow.4$Date)]
  pass.pcb.1$flow.5 <- 0.03*flow.5$X_00060_00003[match(pass.pcb.1$SampleDate,
                                                      flow.5$Date)]
  pass.pcb.1$temp <- 273.15 + temp$X_.from.middle.intake_00010_00003[match(pass.pcb.1$SampleDate,
                                                                          temp$Date)]
  # Remove metadata
  pass.pcb.2 <- subset(pass.pcb.1, select = -c(SiteID:temp))
}

# LME for individual PCBs -------------------------------------------------
# Get covariates
time <- pass.pcb.1$time
season <- pass.pcb.1$season
site <- pass.pcb.1$site.numb
flow <- pass.pcb.1$flow.3 # flow.3 yiels 6 PCB congeners
tem <- pass.pcb.1$temp

# Create matrix to store results
lme.pcb <- matrix(nrow = length(pass.pcb.2[1,]), ncol = 27)

# Perform LME
for (i in 1:length(pass.pcb.2[1,])) {
  fit <- lmer(pass.pcb.2[,i] ~ 1 + time + flow + tem + season + (1|site),
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
  lme.pcb[i,10] <- fixef(fit)[4] # temp
  lme.pcb[i,11] <- summary(fit)$coef[4,"Std. Error"] # temp error
  lme.pcb[i,12] <- summary(fit)$coef[4,"Pr(>|t|)"] # temp p-value
  lme.pcb[i,13] <- fixef(fit)[4] # # season 1
  lme.pcb[i,14] <- summary(fit)$coef[5,"Std. Error"] # season 1 error
  lme.pcb[i,15] <- summary(fit)$coef[5,"Pr(>|t|)"] # # season 1 p-value
  lme.pcb[i,16] <- fixef(fit)[6] # season 2
  lme.pcb[i,17] <- summary(fit)$coef[6,"Std. Error"] # season 2 error
  lme.pcb[i,18] <- summary(fit)$coef[6,"Pr(>|t|)"] # season 2 p-value
  #lme.pcb[i,19] <- fixef(fit)[7] # season 3
  #lme.pcb[i,20] <- summary(fit)$coef[7,"Std. Error"] # season 3 error
  #lme.pcb[i,21] <- summary(fit)$coef[7,"Pr(>|t|)"] # season 3 p-value
  lme.pcb[i,22] <- -log(2)/lme.pcb[i,4]/365 # t0.5
  lme.pcb[i,23] <- abs(-log(2)/lme.pcb[i,4]/365)*lme.pcb[i,5]/abs(lme.pcb[i,4]) # t0.5 error
  lme.pcb[i,24] <- as.data.frame(VarCorr(fit))[1,'sdcor']
  lme.pcb[i,25] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
  lme.pcb[i,26] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
  lme.pcb[i,27] <- shapiro.test(resid(fit))$p.value
}

# Just 3 significant figures
lme.pcb <- formatC(signif(lme.pcb, digits = 3))
# Add congener names
congeners <- colnames(pass.pcb.2)
lme.pcb <- as.data.frame(cbind(congeners, lme.pcb))
# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "flow", "flow.error", "flow.pv", "temp", "temp.error",
                       "temp.pv", "season1", "season1.error", "season1.pv",
                       "season2", "season2.error", "season2, pv", "season3",
                       "season3.error", "season3.pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality")
# Remove congeners with no normal distribution
# Shapiro test p-value < 0.05
lme.pcb$Normality <- as.numeric(lme.pcb$Normality)
# Get the congeners that are not showing normality
lme.pcb.out <- lme.pcb[lme.pcb$Normality < 0.05, ]
lme.pcb <- lme.pcb[lme.pcb$Normality > 0.05, ]

# Export results
write.csv(lme.pcb, file = "Output/Data/Sites/csv/PassaicRiver/PassaicLmePCB.csv")

# Generate predictions
# Select congeners that are not showing normality to be remove from pass.pcb.2
df <- data.frame(names_to_remove = lme.pcb.out$Congeners)
# Get column indices to remove
cols_to_remove <- which(names(pass.pcb.2) %in% df$names_to_remove)
# Remove columns from che.pcb.2 with congeners that don't show normality
pass.pcb.3 <- pass.pcb.2[, -cols_to_remove]

# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(pass.pcb.3[,1]),
                      ncol = length(pass.pcb.3[1,]))

for (i in 1:length(pass.pcb.3[1,])) {
  fit <- lmer(pass.pcb.3[,i] ~ 1 + time + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  lme.fit.pcb[,i] <- fitted(fit)
}

# Estimate a factor of 2 between observations and predictions
factor2 <- 10^(pass.pcb.3)/10^(lme.fit.pcb)
factor2.pcb <- sum(factor2 > 0.5 & factor2 < 2,
                   na.rm = TRUE)/(sum(!is.na(factor2)))*100

# Convert the vector to a data frame
factor2.pcb <- data.frame(Factor_2 = factor2.pcb)

# Export results
write.csv(factor2.pcb,
          file = "Output/Data/Sites/csv/PassaicRiver/PassaicFactor2PCB.csv")

# Individual PCB congener plots -------------------------------------------
# (1)
# Plot 1:1 for all congeners
# Transform lme.fit.pcb to data.frame
lme.fit.pcb <- as.data.frame(lme.fit.pcb)
# Add congener names to lme.fit.pcb columns
colnames(lme.fit.pcb) <- colnames(pass.pcb.3)
# Add code number to first column
df1 <- cbind(code = row.names(pass.pcb.3), pass.pcb.3)
df2 <- cbind(code = row.names(lme.fit.pcb), lme.fit.pcb)

for (i in 2:length(df1)) {
  col_name <- if (i == 1) {
    ""  # leave the name empty for the first plot
  } else {
    names(df1)[i] # use the column name for other plots
  }
  
  # create plot for each pair of columns
  p <- ggplot(data = data.frame(x = df1$code, y1 = 10^(df1[, i]), y2 = 10^(df2[, i])),
              aes(x = y1, y = y2)) +
    geom_point(shape = 21, size = 3, fill = "white") +
    scale_y_log10(limits = c(0.01, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(limits = c(0.01, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
    ylab(expression(bold("Predicted lme concentration PCBi (pg/L)"))) +
    theme_bw() +
    theme(aspect.ratio = 15/15) +
    annotation_logticks(sides = "bl") +
    geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
    geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
    geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
    annotate('text', x = 10^1, y = 10^4, label = gsub("\\.", "+", names(df1)[i]),
             size = 3, fontface = 2)
  # save plot
  ggsave(paste0("Output/Plots/Sites/ObsPred/PassaicRiver/", col_name, ".png"),
         plot = p, width = 6, height = 6, dpi = 500)
}

# (2)
# All plots in one page
# Create a list to store all the plots
plot_list <- list()

# Loop over the columns of df1 and df2
for (i in 2:length(df1)) {
  col_name <- paste(names(df1)[i], sep = "")  # use the column name for plot title
  # Create plot for each pair of columns and add to plot_list
  p <- ggplot(data = data.frame(x = df1$code, y1 = 10^(df1[, i]), y2 = 10^(df2[, i])),
              aes(x = y1, y = y2)) +
    geom_point(shape = 21, size = 3, fill = "white") +
    scale_y_log10(limits = c(0.01, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(limits = c(0.01, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
    ylab(expression(bold("Predicted lme concentration PCBi (pg/L)"))) +
    theme_bw() +
    theme(aspect.ratio = 15/15, 
          axis.title = element_text(size = 8)) +
    annotation_logticks(sides = "bl") +
    annotate('text', x = 25, y = 10^4, label = gsub("\\.", "+", col_name),
             size = 2.5, fontface = 2) +
    geom_abline(intercept = 0, slope = 1, col = "white", linewidth = 0.7) +
    geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
    geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7)
  
  plot_list[[i-1]] <- p  # add plot to list
}
# Combine all the plots using patchwork
combined_plot <- wrap_plots(plotlist = plot_list, ncol = 4)
# Save the combined plot
ggsave("Output/Plots/Sites/ObsPred/PassaicRiver/combined_plot.png", combined_plot,
       width = 15, height = 15, dpi = 500)

# (3)
# Create a list to store all the cleaned data frames
cleaned_df_list <- list()
# Loop over the columns of df1 and df2
for (i in 2:length(df1)) {
  # Create a new data frame by binding the columns of df1 and df2 for each pair of columns
  df_pair <- cbind(df1[,1], df1[,i], df2[,i])
  colnames(df_pair) <- c("code", "observed", "predicted")
  # Remove the rows with missing values
  cleaned_df_pair <- na.omit(df_pair)
  # Add the cleaned data frame to the list
  cleaned_df_list[[i-1]] <- cleaned_df_pair
}

{
  # Modify data to be plotted
  # Combine all the cleaned data frames using rbind
  combined_cleaned_df <- do.call(rbind, cleaned_df_list)
  # Convert the matrix to a data frame
  combined_cleaned_df <- as.data.frame(combined_cleaned_df)
  # Convert the code column to a factor
  combined_cleaned_df$code <- as.factor(combined_cleaned_df$code)
  # Convert the observed and predicted columns to numeric
  combined_cleaned_df[,2:3] <- apply(combined_cleaned_df[,2:3], 2, as.numeric)
}

# Export results for plotting
# Add column LocationName
combined_cleaned_df$LocationName <- "Passaic River"
write.csv(combined_cleaned_df,
          file = "Output/Data/Sites/csv/PassaicRiver/PassaicObsPredPCB.csv")

# Plot all the pairs together
p <- ggplot(combined_cleaned_df, aes(x = 10^(observed), y = 10^(predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(0.01, 10^4), 
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.01, 10^4), 
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration PCBi (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 15/15, 
        axis.title = element_text(size = 10)) +
  annotation_logticks(sides = "bl") +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  annotate("text", x = 0.1, y = 10^3.6,
           label = expression(atop("Passaic River",
                                   paste("6 PCB congeners (n = 77 pairs)"))),
           size = 4, fontface = 2)
# See plot
print(p)
# Save plot
ggsave("Output/Plots/Sites/ObsPred/PassaicRiver/PassaicObsPredPCB.png",
       plot = p, width = 8, height = 8, dpi = 500)

         