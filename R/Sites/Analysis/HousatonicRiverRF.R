## Water PCB concentrations data analysis per site
## Housatonic River
## Aroclors 1254 and 1260, no congener analysis
## GE facility map @ https://semspub.epa.gov/work/01/574882.pdf
## Random forest model

# Install packages
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
install.packages('ranger')
install.packages('caret')

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
  library(ranger) # Random Forest functions
  library(caret) # For cross-validation
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# Select Housatonic River data ---------------------------------------------------
hou <- wdc[str_detect(wdc$LocationName, 'Housatonic River'),]

# Add distance to source --------------------------------------------------
{
  # Define source coordinates
  source1 <- c(Latitude = 42.456479, Longitude = -73.217587) # GE Pittsfield
  # Create an sf point for the source
  source1_sf <- st_sfc(st_point(c(source1["Longitude"], source1["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(source1_sf) <- 4326
  # Transform source1_sf to UTM Zone 18N (EPSG:32618)
  source1_sf_utm <- st_transform(source1_sf, 32618)
  # Convert the data frame to an sf object
  sf_hou <- st_as_sf(hou, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_hou <- st_set_crs(sf_hou, 4326)
  # Transform to UTM Zone 18N
  sf_hou_utm <- st_transform(sf_hou, 32618)
  # Calculate distances in meters from each location to the source
  distances_meters <- st_distance(sf_hou_utm, source1_sf_utm)
  # Convert distances to kilometers
  distances_km <- units::set_units(distances_meters, "km")
  # Extract numeric values and assign to the DistanceToSource column
  hou$DistanceToSource <- as.numeric(distances_km[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  hou$SampleDate <- as.Date(hou$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(hou$SampleDate),
                                  min(as.Date(hou$SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hou$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  hou.tpcb <- cbind(factor(hou$SiteID), hou$SampleDate, as.matrix(hou$tPCB),
                    data.frame(time.day), season.s, hou$DistanceToSource)
  # Add column names
  colnames(hou.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                          "DistanceToSource")
}

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Housatonic River
{
  siteHouN1 <- "01197000" # EAST BRANCH HOUSATONIC RIVER AT COLTSVILLE, MA
  siteHouN2 <- "01197500" # HOUSATONIC RIVER NEAR GREAT BARRINGTON, MA
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  # Retrieve USGS data
  flow.1 <- readNWISdv(siteHouN1, paramflow,
                       min(hou.tpcb$date), max(hou.tpcb$date))
  flow.2 <- readNWISdv(siteHouN2, paramflow,
                       min(hou.tpcb$date), max(hou.tpcb$date))
  # Add USGS data to hou.tpcb, matching dates (m3/s, 0.03 conversion factor)
  hou.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(hou.tpcb$date,
                                                     flow.1$Date)]
  hou.tpcb$flow.2 <- 0.03*flow.2$X_00060_00003[match(hou.tpcb$date,
                                                     flow.2$Date)]
}

# Random Forest Model tPCB ------------------------------------------------
# Remove columns not used here
# Use flow.1
hou.tpcb <- select(hou.tpcb, -c(date, flow.2))

# Train-Test Split
set.seed(123)

train_indices <- sample(1:nrow(hou.tpcb), 0.8 * nrow(hou.tpcb))
train_data <- hou.tpcb[train_indices, ]
test_data <- hou.tpcb[-train_indices, ]

# Define hyperparameter grid
param_grid <- expand.grid(
  mtry = seq(1, ncol(train_data) - 1),  # Adjust mtry values based on your data
  splitrule = c("gini", "extratrees"),
  min.node.size = c(5, 10, 20)
)

# Prepare training control
ctrl <- trainControl(method = "cv", number = 5, search = "grid")

# Perform grid search with cross-validation using ranger
rf_model <- train(
  log10(tPCB) ~ time + SiteID + season + flow.1 + DistanceToSource,
  data = train_data,
  method = "ranger",
  importance = 'permutation',
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best mtry
best_mtry <- rf_model$bestTune$mtry

final_rf_model <- ranger(
  formula = log10(tPCB) ~ time + SiteID + season + flow.1 + DistanceToSource,
  data = train_data,
  num.trees = 1000, # need to manualy modify this parameter
  mtry = best_mtry,
  importance = 'permutation',
  seed = 123
)

# Get predictions on the test data
predictions <- predict(final_rf_model, data = test_data)$predictions

# Evaluate model performance
mse <- mean((predictions - log10(test_data$tPCB))^2)
rmse <- sqrt(mse)

# Calculate R-squared
ss_res <- sum((log10(test_data$tPCB) - predictions)^2)
ss_tot <- sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2)
r_squared <- 1 - (ss_res / ss_tot)

# Print RMSE and R-squared
print(paste("RMSE:", rmse))
print(paste("R-squared:", r_squared))

# Estimate a factor of 2 between observations and predictions
# Create a data frame with observed and predicted values
comparison <- data.frame(observed = test_data$tPCB,
                         predicted = 10^predictions)

# Estimate a factor of 2 between observations and predictions
comparison$factor2 <- comparison$observed/comparison$predicted

# Calculate the percentage of observations within the factor of 2
factor2_percentage <- nrow(comparison[comparison$factor2 > 0.5 & comparison$factor2 < 2
                                      , ])/nrow(comparison)*100

# Print Factor2
print(paste("Factor2:", factor2_percentage))

# Create the data frame directly
performance_RF <- data.frame(Heading = c("RMSE", "R2", "Factor2"),
                             Value = c(rmse, r_squared,
                                       factor2_percentage))

# Export results
write.csv(performance_RF,
          file = "Output/Data/Sites/csv/HousatonicRiver/HousatonicRiverRFtPCB.csv",
          row.names = FALSE)

# Create a data frame for plotting
plot_data <- data.frame(
  Location = rep("Housatonic River", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions
)

# Export results
write.csv(plot_data,
          file = "Output/Data/Sites/csv/HousatonicRiver/HousatonicRiverRFObsPredtPCB.csv",
          row.names = FALSE)

# Create the scatter plot using ggplot2
plotRF <- ggplot(plot_data, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10^3, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^3, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue",
              linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue",
              linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Print the plot
print(plotRF)

# Save plot in folder
ggsave("Output/Plots/Sites/ObsPred/HousatonicRiver/HousatonicRiverRFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)
