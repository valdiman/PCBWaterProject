## Water PCB concentrations data analysis per site
## Lake Michigan Mass Balance & Great Lakes
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
install.packages("gbm3")
install.packages("caret")

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
  library(gbm3) # Random Forest functions
  library(caret) # Random Forest optimization
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# Select LMMB data ---------------------------------------------------
grl <- wdc[str_detect(wdc$LocationName, 'Lake Michigan Mass Balance'), ]

# Just get Lake Michigan data, remove data from tributaries
grl <- grl[!grepl("^Tributary", grl$SiteName), ]

# Calculate central location ---------------------------------------------
{
  # Calculate the mean latitude and longitude
  center_lat <- mean(grl$Latitude)
  center_lon <- mean(grl$Longitude)
  # Create a data frame for the center
  center_df <- data.frame(
    SiteID = "Center",  # You can set a unique identifier for the center
    Latitude = center_lat,
    Longitude = center_lon
  )
  # Convert the center data frame to an sf object
  sf_center <- st_as_sf(center_df, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326) for the center
  sf_center <- st_set_crs(sf_center, 4326)
  # Transform to UTM Zone 16 for the center
  sf_center_utm <- st_transform(sf_center, 32616)
  # Convert the data frame to an sf object
  sf_grl <- st_as_sf(grl, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_grl <- st_set_crs(sf_grl, 4326)
  # Transform to UTM Zone 16
  sf_grl <- st_transform(sf_grl, 32616)
  # Calculate distances in meters from each location to the center
  distances_meters <- st_distance(sf_grl, sf_center_utm)
  # Convert distances to kilometers
  distances_km <- units::set_units(distances_meters, "km")
  # Extract numeric values and assign to the DistanceToCentroid column
  grl$DistanceToCentroid <- as.numeric(distances_km[, 1])
}

# tPCB data preparation ---------------------------------------------------
{
  # Change date format
  grl$SampleDate <- as.Date(grl$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(grl$SampleDate),
                                  min(as.Date(grl$SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(grl$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  grl.tpcb <- cbind(factor(grl$SiteID), grl$SampleDate, as.matrix(grl$tPCB),
                    data.frame(time.day), season.s, grl$DistanceToCentroid)
  # Add column names
  colnames(grl.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                          "DistanceToCentroid")
}

# Add water temperature data ----------------------------------------------
# See code: R/ExtractingData/LMMB/WaterTemp.R
{
  # Read water temperature
  wtp <- read.csv("Output/Data/Sites/csv/GreatLakes/WaterTemp/LakeMichiganWT.csv")
  # Convert date columns to Date format
  wtp$Date <- as.Date(wtp$Date)
  # Add water temperature to grl.tpcb
  grl.tpcb$temp <- wtp$WTMP_K[match(grl.tpcb$date, wtp$Date)]
  # Remove samples with temp = NA
  grl.tpcb <- na.omit(grl.tpcb)
}

# Random Forest Model tPCB ------------------------------------------------
# Train-Test Split
set.seed(123)

# Train-test split
train_indices <- sample(1:nrow(grl.tpcb), 0.8 * nrow(grl.tpcb))
train_data <- grl.tpcb[train_indices, ]
test_data <- grl.tpcb[-train_indices, ]

# Define hyperparameter grid for tuning
param_grid <- expand.grid(
  n.trees = c(500, 1000, 3000, 5000),
  interaction.depth = c(5, 10, 15),
  shrinkage = c(0.001, 0.01, 0.1),
  n.minobsinnode = c(10, 20, 30)  # Add this parameter
)

# Perform grid search with cross-validation
ctrl <- trainControl(method = "cv", number = 5)
gbm_model <- train(
  log10(tPCB) ~ time + SiteID + season + temp,
  data = train_data,
  method = "gbm",
  distribution = "gaussian",
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best hyperparameters
best_n_trees <- gbm_model$bestTune$n.trees
best_depth <- gbm_model$bestTune$interaction.depth
best_shrinkage <- gbm_model$bestTune$shrinkage

# Fit the GBM model with the best hyperparameters
final_gbm_model <- gbm(
  log10(tPCB) ~ time + SiteID + season + temp,
  data = train_data,
  distribution = "gaussian",
  n.trees = best_n_trees,
  interaction.depth = best_depth,
  shrinkage = best_shrinkage
)

# Make predictions on test data
predictions <- predict(object = final_gbm_model, newdata = test_data,
                       n.trees = best_n_trees)

# Evaluate model performance
mse <- mean((predictions - log10(test_data$tPCB))^2)
rmse <- sqrt(mse)
ss_res <- sum((log10(test_data$tPCB) - predictions)^2)
ss_tot <- sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2)
r_squared <- 1 - (ss_res / ss_tot)

# Print results
cat("Best n.trees:", best_n_trees, "\n")
cat("Best interaction depth:", best_depth, "\n")
cat("Best shrinkage:", best_shrinkage, "\n")

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

# Create the data frame directly
performance_RF <- data.frame(Heading = c("RMSE", "R2", "Factor2"),
                             Value = c(rmse, r_squared,
                                       factor2_percentage))
# Print the modified data frame
print(performance_RF)

# Export results
write.csv(performance_RF,
          file = "Output/Data/Sites/csv/GreatLakes/GreatLakesRFtPCBV02.csv",
          row.names = FALSE)

# Create a data frame for plotting
plot_data <- data.frame(
  Location = rep("Great Lakes", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions
)

# Export results
write.csv(plot_data,
          file = "Output/Data/Sites/csv/GreatLakes/GreatLakesRFObsPredtPCBV02.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(plot_data, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10, 10^4),
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
ggsave("Output/Plots/Sites/ObsPred/GreatLakes/GreatLakesRFtPCBV02.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Random Forest doesn't work for tPCB, R2<0!

