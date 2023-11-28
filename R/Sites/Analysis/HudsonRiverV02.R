## Water PCB concentrations data analysis per site
## Hudson River
## Model 4 shows better performance

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

# Select Hudson River data ---------------------------------------------------
hud <- wdc[str_detect(wdc$LocationName, 'Hudson River'),]
# PCBs were discharged to the river from the General Electric
# (GE) manufacturing plants in Hudson Falls and Fort Edward, NY
# Dredging from 2009 to 2015
# https://www.epa.gov/system/files/documents/2021-08/hudson_summer2021_floodplainrifs_factsheet_final.pdf

# Data preparation --------------------------------------------------------
{
  # Change date format
  hud$SampleDate <- as.Date(hud$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(hud$SampleDate) - min(as.Date(hud$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- hud$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hud$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  hud.tpcb <- cbind(factor(hud$SiteID), hud$SampleDate,
                    hud$Latitude, hud$Longitude, as.matrix(hud$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(hud.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Source locations
source1 <- c(Latitude = 43.295369, Longitude = -73.590631)  # GE Hudson Falls Plant
source2 <- c(Latitude = 43.28639, Longitude = -73.588380)  # GE Fort Edward Plant

# Create sf points
source1_sf <- st_point(c(source1["Latitude"], source1["Longitude"]))
source2_sf <- st_point(c(source2["Latitude"], source2["Longitude"]))

# Initialize a list to store point geometries
point_list <- list()

# Loop to create point geometries
for (i in 1:nrow(hud.tpcb)) {
  point <- st_point(c(hud.tpcb[i, "Latitude"], hud.tpcb[i, "Longitude"]))
  point_list[[i]] <- point
}

# Create an sf object from the list of point geometries
hud.tpcb_sf <- st_sf(geometry = st_sfc(point_list), crs = 4326)

# Ensure that both objects have the same CRS
st_crs(hud.tpcb_sf) <- st_crs(source1_sf)
st_crs(hud.tpcb_sf) <- st_crs(source2_sf)

# Calculate distances
distances1 <- st_distance(hud.tpcb_sf, source1_sf) * 100
distances2 <- st_distance(hud.tpcb_sf, source2_sf) * 100

# Add distances to hud.tpcb data.frame
hud.tpcb$Distance_to_source1 <- distances1
hud.tpcb$Distance_to_source2 <- distances2

# Remove site -------------------------------------------------------------
## Remove site Bakers Falls. Upstream source
## North Bakers Falls = WCPCB-HUD006 and
## South Bakers Falls = WCPCB-HUD006.
hud.tpcb.1 <- subset(hud.tpcb, SiteID != c("WCPCB-HUD006"))
hud.tpcb.1 <- subset(hud.tpcb.1, SiteID != c("WCPCB-HUD010"))

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Hudson River
{
  sitehudN1 <- "01331095" # HUDSON RIVER AT STILLWATER NY No temp!
  sitehudN2 <- "01335754" # HUDSON RIVER ABOVE LOCK 1 NEAR WATERFORD NY, no temp!
  sitehudN3 <- "01328770" # HUDSON RIVER AT THOMSON NY, no temp!
  sitehudN4 <- "01327750" # HUDSON RIVER AT FORT EDWARD NY, no temp!
  sitehudN5 <- "01359139" # HUDSON RIVER AT ALBANY NY No flow!
  
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  
  # Flow (ft3/s)
  flow.1 <- readNWISdv(sitehudN1, paramflow,
                       min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  flow.2 <- readNWISdv(sitehudN2, paramflow,
                       min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  flow.3 <- readNWISdv(sitehudN3, paramflow,
                       min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  flow.4 <- readNWISdv(sitehudN4, paramflow,
                       min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  # Water temperature in Celsius
  temp <- readNWISdv(sitehudN5, paramtemp,
                     min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  
  # Add USGS data to hud.tpcb.2, matching dates
  hud.tpcb.1$flow.1 <- 0.03*flow.1$X_00060_00003[match(hud.tpcb.1$date,
                                                       flow.1$Date)]
  hud.tpcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(hud.tpcb.1$date,
                                                       flow.2$Date)]
  hud.tpcb.1$flow.3 <- 0.03*flow.3$X_00060_00003[match(hud.tpcb.1$date,
                                                       flow.3$Date)]
  hud.tpcb.1$flow.4 <- 0.03*flow.4$X_00060_00003[match(hud.tpcb.1$date,
                                                       flow.4$Date)]
  hud.tpcb.1$temp <- 273.15 + temp$X_00010_00003[match(hud.tpcb.1$date,
                                                       temp$Date)]
  # Remove samples with temp = NA
  hud.tpcb.2 <- na.omit(hud.tpcb.1)
}

# Random Forest Model -----------------------------------------------------
# Using all the data, hud.tpcb.1
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(hud.tpcb.1), 0.8 * nrow(hud.tpcb.1))
train_data <- hud.tpcb.1[train_indices, ]
test_data <- hud.tpcb.1[-train_indices, ]

# Handle missing values from flow.3
train_data <- na.omit(train_data)
test_data <- na.omit(test_data)

# Fit the Model (1)
rf_model.1 <- randomForest(log10(tPCB) ~ time + site.code + season + flow.1 + flow.2
                        + flow.4 + temp + Distance_to_source1 +
                           Distance_to_source2, data = train_data)

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
  scale_y_log10(limits = c(10^2, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^2, 10^6),
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

# Fit the Model (2)
rf_model.2 <- randomForest(log10(tPCB) ~ time + site.code + flow.2
                         + temp + Distance_to_source2,
                         data = train_data)

# Make Predictions
predictions.2 <- predict(rf_model.2, newdata = test_data)

# Evaluate Model Performance
mse.2 <- mean((predictions.2 - log10(test_data$tPCB))^2)
rmse.2 <- sqrt(mse.2)
mae.2 <- mean(abs(predictions.2 - log10(test_data$tPCB)))
r_squared.2 <- 1 - (sum((log10(test_data$tPCB) - predictions.2)^2)/sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2))

# Feature Importance
importance.2 <- importance(rf_model.2)
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
  scale_y_log10(limits = c(10^2, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^2, 10^6),
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

# Fit the Model (3)
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(hud.tpcb.2), 0.8 * nrow(hud.tpcb.2))
train_data <- hud.tpcb.2[train_indices, ]
test_data <- hud.tpcb.2[-train_indices, ]

# Handle missing values from flow.3
train_data <- na.omit(train_data)
test_data <- na.omit(test_data)

# Fit the Model
rf_model.3 <- randomForest(log10(tPCB) ~ time + site.code + season + flow.1 + flow.2
                           + flow.3 + flow.4 + temp + Distance_to_source1 +
                             Distance_to_source2, data = train_data)

# Make Predictions
predictions.3 <- predict(rf_model.3, newdata = test_data)

# Evaluate Model Performance
mse.3 <- mean((predictions.3 - log10(test_data$tPCB))^2)
rmse.3 <- sqrt(mse.3)
mae.3 <- mean(abs(predictions.3 - log10(test_data$tPCB)))
r_squared.3 <- 1 - (sum((log10(test_data$tPCB) - predictions.3)^2)/sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2))

# Feature Importance
importance.3 <- importance(rf_model.3)
barplot(importance.3[, 1], names.arg = rownames(importance.3),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting
plot_data.3 <- data.frame(
  Actual = log10(test_data$tPCB),
  Predicted = predictions.3
)

# Create the scatter plot using ggplot2
ggplot(plot_data.3, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10^2, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^2, 10^6),
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

# Fit the Model (4)
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(hud.tpcb.2), 0.8 * nrow(hud.tpcb.2))
train_data <- hud.tpcb.2[train_indices, ]
test_data <- hud.tpcb.2[-train_indices, ]

# Handle missing values from flow.3
train_data <- na.omit(train_data)
test_data <- na.omit(test_data)

# Fit the Model (4)
rf_model.4 <- randomForest(log10(tPCB) ~ time + site.code + season +
                             flow.3 + temp + Distance_to_source1,
                           data = train_data)

# Make Predictions
predictions.4 <- predict(rf_model.4, newdata = test_data)

# Evaluate Model Performance
mse.4 <- mean((predictions.4 - log10(test_data$tPCB))^2)
rmse.4 <- sqrt(mse.4)
mae.4 <- mean(abs(predictions.4 - log10(test_data$tPCB)))
r_squared.4 <- 1 - (sum((log10(test_data$tPCB) - predictions.4)^2)/sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2))

# Feature Importance
importance.4 <- importance(rf_model.4)
barplot(importance.4[, 1], names.arg = rownames(importance.4),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting
plot_data.4 <- data.frame(
  Actual = log10(test_data$tPCB),
  Predicted = predictions.4
)

# Create the scatter plot using ggplot2
ggplot(plot_data.4, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10^2, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^2, 10^6),
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

# Fit the Model (5)
rf_model.5 <- randomForest(log10(tPCB) ~ time + site.code + flow.3 +
                             temp + Distance_to_source1,
                           data = train_data)

# Make Predictions
predictions.5 <- predict(rf_model.5, newdata = test_data)

# Evaluate Model Performance
mse.5 <- mean((predictions.5 - log10(test_data$tPCB))^2)
rmse.5 <- sqrt(mse.5)
mae.5 <- mean(abs(predictions.5 - log10(test_data$tPCB)))
r_squared.5 <- 1 - (sum((log10(test_data$tPCB) - predictions.5)^2)/sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2))

# Feature Importance
importance.5 <- importance(rf_model.5)
barplot(importance.5[, 1], names.arg = rownames(importance.5),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting
plot_data.5 <- data.frame(
  Actual = log10(test_data$tPCB),
  Predicted = predictions.5
)

# Create the scatter plot using ggplot2
ggplot(plot_data.5, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10^2, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^2, 10^6),
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

