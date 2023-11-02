## Water PCB concentrations data analysis per site
## Lake Michigan Mass Balance & Great Lakes

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

# Select LMMB and Great Lakes data ---------------------------------------------------
grl <- wdc[str_detect(wdc$LocationName, 'Lake Michigan Mass Balance|Great Lakes'), ]

# Just open lake data, remove data from tributaries
grl <- grl[!grepl("^Tributary", grl$SiteName), ]

# Data preparation --------------------------------------------------------
{
  # Change date format
  grl$SampleDate <- as.Date(grl$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(grl$SampleDate) - min(as.Date(grl$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- grl$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(grl$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  grl.tpcb <- cbind(factor(grl$SiteID), grl$SampleDate,
                    grl$Latitude, grl$Longitude, as.matrix(grl$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(grl.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Get coordinates per site to plot in Google Earth
location <- grl.tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = location, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/GreatLakesLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# General plots -------------------------------------------------------------------
# (1) Histograms
hist(grl.tpcb$tPCB)
hist(log10(grl.tpcb$tPCB))

# (2) Time trend plots
GLTime <- ggplot(grl.tpcb, aes(y = tPCB, x = format(date, '%Y'))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  xlab("") +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000, 10000, 100000),  # Specify the desired breaks
    labels = label_comma()(c(1, 10, 100, 1000, 10000, 100000))  # Specify the desired labels
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
print(GLTime)

# Save plot in folder
ggsave("Output/Plots/Sites/Temporal/plotGreatLakesTime.png",
       plot = GLTime, width = 6, height = 5, dpi = 500)

# (3) Seasonality
ggplot(grl.tpcb, aes(x = season, y = tPCB)) +
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
  annotate("text", x = 1, y = 20, label = "Great Lakes",
           size = 3)

# (4) Sites
ggplot(grl.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
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
  annotate("text", x = 5, y = 20, label = "Great Lakes",
           size = 3)

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- grl.tpcb$tPCB
time <- grl.tpcb$time
site <- grl.tpcb$site.code
season <- grl.tpcb$season
# tPCB vs. time + season + flow + temp + site
lme.grl.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.grl.tpcb)
# Look at residuals
{
  res.grl.tpcb <- resid(lme.grl.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/GreatLakesQ-QtPCB.pdf")
  qqnorm(res.grl.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.grl.tpcb)
  dev.off()
}
# Shapiro test
shapiro.test(resid(lme.grl.tpcb))

# Remove samples less than 10 and more than 30,000 pg/L (n = 9)
grl.tpcb.1 <- subset(grl.tpcb, tPCB > 10 & tPCB < 30000)

# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- grl.tpcb.1$tPCB
time <- grl.tpcb.1$time
site <- grl.tpcb.1$site.code
season <- grl.tpcb.1$season
# tPCB vs. time + season + flow + temp + site
lme.grl.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.grl.tpcb)
# Look at residuals
{
  res.grl.tpcb <- resid(lme.grl.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/GreatLakesQ-QtPCBV02.pdf")
  qqnorm(res.grl.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.grl.tpcb)
  dev.off()
}
# Shapiro test
shapiro.test(resid(lme.grl.tpcb))

# Create matrix to store results
{
  lme.tpcb <- matrix(nrow = 1, ncol = 18)
  lme.tpcb[1] <- fixef(lme.grl.tpcb)[1] # intercept
  lme.tpcb[2] <- summary(lme.grl.tpcb)$coef[1,"Std. Error"] # intercept error
  lme.tpcb[3] <- summary(lme.grl.tpcb)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb[4] <- fixef(lme.grl.tpcb)[2] # time
  lme.tpcb[5] <- summary(lme.grl.tpcb)$coef[2,"Std. Error"] # time error
  lme.tpcb[6] <- summary(lme.grl.tpcb)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb[7] <- fixef(lme.grl.tpcb)[3] # season 2
  lme.tpcb[8] <- summary(lme.grl.tpcb)$coef[3,"Std. Error"] # season 2 error
  lme.tpcb[9] <- summary(lme.grl.tpcb)$coef[3,"Pr(>|t|)"] # season 2 p-value
  lme.tpcb[10] <- fixef(lme.grl.tpcb)[4] # season 3
  lme.tpcb[11] <- summary(lme.grl.tpcb)$coef[4,"Std. Error"] # season 3 error
  lme.tpcb[12] <- summary(lme.grl.tpcb)$coef[4,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb[13] <- -log(2)/lme.tpcb[4]/365 # t0.5
  lme.tpcb[14] <- abs(-log(2)/lme.tpcb[4]/365)*lme.tpcb[5]/abs(lme.tpcb[4]) # t0.5 error
  lme.tpcb[15] <- as.data.frame(VarCorr(lme.grl.tpcb))[1,'sdcor']
  lme.tpcb[16] <- as.data.frame(r.squaredGLMM(lme.grl.tpcb))[1, 'R2m']
  lme.tpcb[17] <- as.data.frame(r.squaredGLMM(lme.grl.tpcb))[1, 'R2c']
  lme.tpcb[18] <- shapiro.test(resid(lme.grl.tpcb))$p.value
}

# Just 3 significant figures
lme.tpcb <- formatC(signif(lme.tpcb, digits = 3))
# Add column names
colnames(lme.tpcb) <- c("Intercept", "Intercept.error",
                        "Intercept.pv", "time", "time.error", "time.pv",
                        "season2", "season2.error", "season2.pv", "season3",
                        "season3.error", "season3.pv", "t05", "t05.error",
                        "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality")

# Export results
write.csv(lme.tpcb, file = "Output/Data/Sites/csv/GreatLakes/GreatLakesLmetPCB.csv")

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.grl.tpcb <- as.data.frame(fitted(lme.grl.tpcb))
# Add column name
colnames(fit.lme.values.grl.tpcb) <- c("predicted")
# Add predicted values to data.frame
grl.tpcb.1$predicted <- 10^(fit.lme.values.grl.tpcb$predicted)
# Create overall plot prediction vs. observations
predic.obs <- data.frame(tPCB = grl.tpcb.1$tPCB, predicted = grl.tpcb.1$predicted)
predic.obs <- data.frame(Location = grl$LocationName[1], predic.obs)
# Save new data
write.csv(predic.obs, "Output/Data/Sites/csv/GreatLakes/GreatLakesPredic_Obser.csv")

# Plot prediction vs. observations, 1:1 line
p <- ggplot(grl.tpcb.1, aes(x = tPCB, y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(1, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 10, y = 10^3.5,
           label = expression(atop(" Great Lakes (R"^2*"= 0.76)",
                                   paste("t"[1/2]*" = 45 ± 13 (yr)"))),
           size = 4, fontface = 2)
# See plot
print(p)
# Save plot
ggsave("Output/Plots/Sites/ObsPred/GreatLakes/GreatLakesObsPredtPCB.png",
       plot = p, width = 8, height = 8, dpi = 500)

# Plot residuals vs. predictions
{
  # Open a PNG graphics device
  png("Output/Plots/Sites/Residual/GreatLakesResidualtPCB.png", width = 800,
      height = 600)
  # Create your plot
  plot(grl.tpcb.1$predicted, resid(lme.grl.tpcb),
       points(grl.tpcb.1$predicted, resid(lme.grl.tpcb), pch = 16, 
              col = "white"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  # Add lines to the plot
  abline(0, 0)
  abline(h = seq(-2, 2, 1), col = "grey")
  abline(v = seq(0, 1000, 200), col = "grey")
  # Close the PNG graphics device
  dev.off()
}

# Estimate a factor of 2 between observations and predictions
grl.tpcb.1$factor2 <- grl.tpcb.1$tPCB/grl.tpcb.1$predicted
factor2.tpcb <- nrow(grl.tpcb.1[grl.tpcb.1$factor2 > 0.5 & grl.tpcb.1$factor2 < 2,
                                ])/length(grl.tpcb.1[,1])*100

# Individual PCB Analysis -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  grl.pcb <- subset(grl, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  grl.pcb <- subset(grl.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  grl.pcb <- log10(grl.pcb)
  # Replace -inf to NA
  grl.pcb <- do.call(data.frame,
                     lapply(grl.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  grl.pcb.1 <- grl.pcb[,
                       -which(colSums(is.na(grl.pcb))/nrow(grl.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(grl$SiteID)
  # Change date format
  SampleDate <- as.Date(grl$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Create individual code for each site sampled
  site.numb <- grl$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(grl$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to fox.pcb.1
  grl.pcb.1 <- cbind(grl.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     site.numb, season.s)
  # Remove metadata
  grl.pcb.2 <- subset(grl.pcb.1, select = -c(SiteID:season.s))
}

# LME for individual PCBs -------------------------------------------------
# Get covariates
time <- grl.pcb.1$time
season <- grl.pcb.1$season
site <- grl.pcb.1$site.numb

# Create matrix to store results
lme.pcb <- matrix(nrow = length(grl.pcb.2[1,]), ncol = 18)

# Perform LME
for (i in 1:length(grl.pcb.2[1,])) {
  fit <- lmer(grl.pcb.2[,i] ~ 1 + time + season + (1|site),
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
  lme.pcb[i,7] <- fixef(fit)[3] # # season 2
  lme.pcb[i,8] <- summary(fit)$coef[3,"Std. Error"] # season 2 error
  lme.pcb[i,9] <- summary(fit)$coef[3,"Pr(>|t|)"] # # season 2 p-value
  lme.pcb[i,10] <- fixef(fit)[4] # season 3
  lme.pcb[i,11] <- summary(fit)$coef[4,"Std. Error"] # season 3 error
  lme.pcb[i,12] <- summary(fit)$coef[4,"Pr(>|t|)"] # season 3 p-value
  lme.pcb[i,13] <- -log(2)/lme.pcb[i,4]/365 # t0.5
  lme.pcb[i,14] <- abs(-log(2)/lme.pcb[i,4]/365)*lme.pcb[i,5]/abs(lme.pcb[i,4]) # t0.5 error
  lme.pcb[i,15] <- as.data.frame(VarCorr(fit))[1,'sdcor']
  lme.pcb[i,16] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
  lme.pcb[i,17] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
  lme.pcb[i,18] <- shapiro.test(resid(fit))$p.value
}

# Just 3 significant figures
lme.pcb <- formatC(signif(lme.pcb, digits = 3))
# Add congener names
congeners <- colnames(grl.pcb.2)
lme.pcb <- as.data.frame(cbind(congeners, lme.pcb))
# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "season2", "season2.error", "season2, pv", "season3",
                       "season3.error", "season3.pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality")

# Remove congeners with no normal distribution
# Shapiro test p-value < 0.05
lme.pcb$Normality <- as.numeric(lme.pcb$Normality)
# Get the congeners that are not showing normality
lme.pcb.out <- lme.pcb[lme.pcb$Normality < 0.045, ]
lme.pcb <- lme.pcb[lme.pcb$Normality > 0.045, ]

# Only one congener and just above 0.05

# Export results
write.csv(lme.pcb, file = "Output/Data/Sites/csv/GreatLakes/GreatLakesLmePCB.csv")

# Generate predictions
# Select congeners that are not showing normality to be remove from pass.pcb.2
df <- data.frame(names_to_remove = lme.pcb.out$Congeners)
# Get column indices to remove
cols_to_remove <- which(names(grl.pcb.2) %in% df$names_to_remove)
# Remove columns from che.pcb.2 with congeners that don't show normality
grl.pcb.3 <- grl.pcb.2[, -cols_to_remove]

# Create vector to store results
lme.fit.pcb <- rep(NA, 333)

fit <- lmer(grl.pcb.3 ~ 1 + time + season + (1|site),
            REML = FALSE,
            control = lmerControl(check.nobs.vs.nlev = "ignore",
                                  check.nobs.vs.rankZ = "ignore",
                                  check.nobs.vs.nRE="ignore"),
            na.action = na.exclude)
lme.fit.pcb <- fitted(fit)

# Estimate a factor of 2 between observations and predictions
factor2 <- 10^(grl.pcb.3)/10^(lme.fit.pcb)
factor2.pcb <- sum(factor2 > 0.5 & factor2 < 2,
                   na.rm = TRUE)/(sum(!is.na(factor2)))*100

# Individual PCB congener plots -------------------------------------------
# (1)
# Plot 1:1 for all congeners
# Transform lme.fit.pcb to data.frame
lme.fit.pcb <- as.data.frame(lme.fit.pcb)
# Add congener names to lme.fit.pcb columns
colnames(lme.fit.pcb) <- colnames(grl.pcb.3)
# Add code number to first column
df1 <- cbind(code = row.names(grl.pcb.3), grl.pcb.3)
df2 <- cbind(code = row.names(lme.fit.pcb), lme.fit.pcb)

# Assuming df1 is a numeric vector, and df2 is a data frame
# Convert df1 to a data frame
df1 <- data.frame(code = 1:333, grl.pcb.3 = df1)

# Combine df1 and df2 into a single data frame
combined_df <- cbind(df1, df2)

# Remove rows with NAs from combined_df
combined_df <- na.omit(combined_df)

# Make column names unique
colnames(combined_df) <- make.unique(colnames(combined_df))

# Update combined_df
# Identify the column indices you want to delete
columns_to_delete <- c("code", "code.1")

# Delete the specified columns
combined_df <- combined_df[, !(names(combined_df) %in% columns_to_delete)]

# Define new column names
new_names <- c("observed", "predicted")

# Set the new column names
colnames(combined_df) <- new_names

# Add column LocationName
combined_df$LocationName <- "Great Lakes"

# Create the plot
col_name <- "PCB89"  # Set the name to "PCB89"
p <- ggplot(data = combined_df, aes(x = 10^(observed), y = 10^(predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(0.01, 10^2), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.01, 10^2), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration PCBi (pg/L)")))+
  ylab(expression(bold("Predicted lme concentration PCBi (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2) +
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  annotate('text', x = 0.1, y = 10^2, label = gsub("\\.", "+", col_name),
           size = 3, fontface = 2)

print(p)

# Save the plot with the name "PCB89"
ggsave("Output/Plots/Sites/ObsPred/Greatlakes/PCB89.png", plot = p,
       width = 6, height = 6, dpi = 500)

# Export results for plotting
write.csv(combined_df,
          file = "Output/Data/Sites/csv/GreatLakes/ObsPredGreatLakesPCB.csv")




