## Water PCB concentrations data analysis per site
## Richardson Hill Road Landfill
## Arolclor method
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

# Select Richardson Hill Road Landfill data ---------------------------------------------------
rhr <- wdc[str_detect(wdc$LocationName, 'Richardson Hill Road Landfill'),]

# Located southern location & calculate distance to other locations -------
{
  # Identify the northern sample based on the maximum Longitude
  index_of_northern_sample <- which.min(rhr$Longitude)
  # Extract coordinates for the northern sample
  northern_sample_latitude <- rhr$Latitude[index_of_northern_sample]
  northern_sample_longitude <- rhr$Longitude[index_of_northern_sample]
  # Define source coordinates for the northern sample
  northern_source <- c(Latitude = northern_sample_latitude,
                       Longitude = northern_sample_longitude)
  # Create an sf point for the northern source
  northern_source_sf <- st_sfc(st_point(c(northern_source["Longitude"],
                                          northern_source["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(northern_source_sf) <- 4326
  # Transform northern_source_sf to UTM Zone 18N (EPSG:32610)
  northern_source_sf_utm <- st_transform(northern_source_sf, 32618)
  # Convert the data frame to an sf object for the northern sample (rhr)
  sf_rhr <- st_as_sf(rhr, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_rhr <- st_set_crs(sf_rhr, 4326)
  # Transform to UTM Zone 18N
  sf_rhr_utm <- st_transform(sf_rhr, 32618)
  # Calculate distances in meters from each location to northern source
  distances_meters_rhr <- st_distance(sf_rhr_utm, northern_source_sf_utm)
  # Convert distances to kilometers
  distances_km_rhr <- units::set_units(distances_meters_rhr, "km")
  # Extract numeric values and assign to the DistanceToSouthernSource column
  rhr$DistanceToSouthernLocation <- as.numeric(distances_km_rhr[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  rhr$SampleDate <- as.Date(rhr$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(rhr$SampleDate),
                                  min(as.Date(rhr$SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(rhr$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  rhr.tpcb <- cbind(factor(rhr$SiteID), as.matrix(rhr$tPCB), data.frame(time.day),
                    season.s, rhr$DistanceToSouthernLocation)
  # Add column names
  colnames(rhr.tpcb) <- c("SiteID", "tPCB", "time", "season",
                          "DistanceToSouthernLocation")
}

# Random Forest Model tPCB ------------------------------------------------
# Set seed for reproducibility
set.seed(123)

# Train-Test Split
train_indices <- sample(1:nrow(rhr.tpcb), 0.8 * nrow(rhr.tpcb))
train_data <- rhr.tpcb[train_indices, ]
test_data <- rhr.tpcb[-train_indices, ]

# Define hyperparameter grid
param_grid <- expand.grid(
  mtry = seq(1, ncol(train_data) - 1),  # Adjust mtry values based on your data
  splitrule = c("gini", "extratrees"),
  min.node.size = c(2, 3, 5)
)

# Prepare training control
ctrl <- trainControl(method = "cv", number = 5, search = "grid")

# Perform grid search with cross-validation using ranger
rf_model <- train(
  log10(tPCB) ~ time + SiteID + season + DistanceToSouthernLocation,
  data = train_data,
  method = "ranger",
  importance = 'permutation',
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best mtry
best_mtry <- rf_model$bestTune$mtry

final_rf_model <- ranger(
  formula = log10(tPCB) ~ time + SiteID + season + DistanceToSouthernLocation,
  data = train_data,
  num.trees = 100, # need to manually modify this parameter
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
          file = "Output/Data/Sites/csv/Richardson/RichardsonRFtPCB.csv",
          row.names = FALSE)

# Create a data frame for plotting
plot_data <- data.frame(
  Location = rep("Richardson Hill Road Landfill", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions
)

# Export results
write.csv(plot_data,
          file = "Output/Data/Sites/csv/Richardson/RichardsonRFObsPredtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(plot_data, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 1, fill = "white") +
  scale_y_log10(limits = c(100, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(100, 10^7),
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
ggsave("Output/Plots/Sites/ObsPred/Richardson/RichardsonRFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# It doesn't work. R2 < 0.0