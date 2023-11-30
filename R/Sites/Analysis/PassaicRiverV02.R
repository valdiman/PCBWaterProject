## Water PCB concentrations data analysis per site
## Passaic River
## Source: Diamond Alkali. 80-120 Lister Avenue, Newark, NJ.

# Install packages
install.packages("randomForest")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("robustbase")
install.packages("dplyr")
install.packages("tibble")
install.packages("Matrix")
install.packages("zoo")
install.packages("dataRetrieval")
install.packages("reshape")
install.packages("tidyr")
install.packages('patchwork')
install.packages("scales")
install.packages("sf")
install.packages("units")
install.packages("sfheaders")

# Load libraries
{
  library(ggplot2)
  library(scales) # function trans_breaks
  library(stringr) # str_detect
  library(robustbase) # function colMedians
  library(dplyr) # performs %>%
  library(tibble) # adds a column
  library(zoo) # yields seasons
  library(dataRetrieval) # read data from USGS
  library(reshape)
  library(tidyr) # function gather
  library(patchwork) # combine plots
  library(sf) # Create file to be used in Google Earth
  library(units)
  library(randomForest)
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

# Source locations
source1 <- c(Latitude = 40.740238, Longitude = -74.133805)  # Former Diamond Alkali

# Create sf points
source1_sf <- st_point(c(source1["Latitude"], source1["Longitude"]))

# Initialize a list to store point geometries
point_list <- list()

# Loop to create point geometries
for (i in 1:nrow(nbh.tpcb)) {
  point <- st_point(c(nbh.tpcb[i, "Latitude"], nbh.tpcb[i, "Longitude"]))
  point_list[[i]] <- point
}

# Create an sf object from the list of point geometries
nbh.tpcb_sf <- st_sf(geometry = st_sfc(point_list), crs = 4326)

# Ensure that both objects have the same CRS
st_crs(nbh.tpcb_sf) <- st_crs(source1_sf)

# Calculate distances
distances <- st_distance(nbh.tpcb_sf, source1_sf) * 100

# Add distances to hud.tpcb data.frame
nbh.tpcb$Distance_to_source1 <- distances




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

# Random Forest Model -----------------------------------------------------
# Using all the data, hud.tpcb.1
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(pass.tpcb), 0.8 * nrow(pass.tpcb))
train_data <- pass.tpcb[train_indices, ]
test_data <- pass.tpcb[-train_indices, ]

# Fit the Model (1)
rf_model.1 <- randomForest(log10(tPCB) ~ time + site.code + season +
                             flow.1 + flow.2 + flow.3 + flow.4 + flow.5 +
                             temp, data = train_data)

# Make Predictions
predictions.1 <- predict(rf_model.1, newdata = test_data)

# Evaluate Model Performance
mse.1 <- mean((predictions.1 - log10(test_data$tPCB))^2)
rmse.1 <- sqrt(mse.1)
mae.1 <- mean(abs(predictions.1 - log10(test_data$tPCB)))
r_squared.1 <- 1 - (sum((log10(test_data$tPCB) - predictions.1)^2)/sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2))

# Feature Importance
importance.1 <- importance(rf_model.1)
# Plot features
barplot(importance.1[, 1], names.arg = rownames(importance.1),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting
plot_data.1 <- data.frame(
  Actual = log10(test_data$tPCB),
  Predicted = predictions.1
)

# Create the scatter plot using ggplot2
ggplot(plot_data.1, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10^2, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^2, 10^7),
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
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Estimate a factor of 2 between observations and predictions
# Create a data frame with observed and predicted values
compare_df.1 <- data.frame(observed = test_data$tPCB,
                         predicted = 10^predictions.1)

# Estimate a factor of 2 between observations and predictions
compare_df.1$factor2 <- compare_df.1$observed/compare_df.1$predicted

# Calculate the percentage of observations within the factor of 2
factor2_percentage.1 <- nrow(compare_df.1[compare_df.1$factor2 > 0.5 & compare_df.1$factor2 < 2, ])/nrow(compare_df.1)*100

# Fit the Model (2)
rf_model.2 <- randomForest(log10(tPCB) ~ time + site.code + season +
                             flow.1 + flow.2 + flow.3 + flow.4 + flow.5 +
                             temp, data = train_data)

# Make Predictions
predictions.2 <- predict(rf_model.2, newdata = test_data)

# Evaluate Model Performance
mse.2 <- mean((predictions.2 - log10(test_data$tPCB))^2)
rmse.2 <- sqrt(mse.2)
mae.2 <- mean(abs(predictions.2 - log10(test_data$tPCB)))
r_squared.2 <- 1 - (sum((log10(test_data$tPCB) - predictions.2)^2)/sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2))

# Feature Importance
importance.2 <- importance(rf_model.2)
# Plot features
barplot(importance.2[, 1], names.arg = rownames(importance.2),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting
plot_data.2 <- data.frame(
  Actual = log10(test_data$tPCB),
  Predicted = predictions.2
)

# Create the scatter plot using ggplot2
ggplot(plot_data.2, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10^2, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^2, 10^7),
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
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Estimate a factor of 2 between observations and predictions
# Create a data frame with observed and predicted values
compare_df.2 <- data.frame(observed = test_data$tPCB,
                           predicted = 10^predictions.2)

# Estimate a factor of 2 between observations and predictions
compare_df.2$factor2 <- compare_df.2$observed/compare_df.2$predicted

# Calculate the percentage of observations within the factor of 2
factor2_percentage.2 <- nrow(compare_df.2[compare_df.2$factor2 > 0.5 & compare_df.2$factor2 < 2, ])/nrow(compare_df.2)*100

