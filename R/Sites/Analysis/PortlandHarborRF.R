## Water PCB concentrations data analysis per site
## Portland Harbor
## Random Forest model

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

# Select Portland Harbor data ---------------------------------------------------
por <- wdc[str_detect(wdc$LocationName, 'Portland Harbor'),]

# Located northern location & calculate distance to other locations -------
{
  # Identify the northern sample based on the maximum Longitude
  index_of_northern_sample <- which.max(por$Longitude)
  # Extract coordinates for the northern sample
  northern_sample_latitude <- por$Latitude[index_of_northern_sample]
  northern_sample_longitude <- por$Longitude[index_of_northern_sample]
  # Define source coordinates for the northern sample
  northern_source <- c(Latitude = northern_sample_latitude,
                       Longitude = northern_sample_longitude)
  # Create an sf point for the northern source
  northern_source_sf <- st_sfc(st_point(c(northern_source["Longitude"],
                                          northern_source["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(northern_source_sf) <- 4326
  # Transform northern_source_sf to UTM Zone 10N (EPSG:32610)
  northern_source_sf_utm <- st_transform(northern_source_sf, 32610)
  # Convert the data frame to an sf object for the northern sample (spo)
  sf_por <- st_as_sf(por, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_por <- st_set_crs(sf_por, 4326)
  # Transform to UTM Zone 10N
  sf_por_utm <- st_transform(sf_por, 32610)
  # Calculate distances in meters from each location to northern source
  distances_meters_por <- st_distance(sf_por_utm, northern_source_sf_utm)
  # Convert distances to kilometers
  distances_km_por <- units::set_units(distances_meters_por, "km")
  # Extract numeric values and assign to the DistanceToNorthernSource column
  por$DistanceToNorthernLocation <- as.numeric(distances_km_por[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  por$SampleDate <- as.Date(por$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(por$SampleDate),
                                  min(as.Date(por$SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(por$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  por.tpcb <- cbind(factor(por$SiteID), por$SampleDate, as.matrix(por$tPCB),
                    data.frame(time.day), season.s, por$DistanceToNorthernLocation)
  # Add column names
  colnames(por.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                          "DistanceToNorthernLocation")
  # Include USGC station Portland Harbor flow and water temperature
  sitePorN1 <- "14211720" # WILLAMETTE RIVER AT PORTLAND, OR
  sitePorN2 <- "14211820" # COLUMBIA SLOUGH AT PORTLAND, OR
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  # Flow (ft3/s)
  flow.1 <- readNWISdv(sitePorN1, paramflow,
                       min(por.tpcb$date), max(por.tpcb$date))
  flow.2 <- readNWISdv(sitePorN2, paramflow,
                       min(por.tpcb$date), max(por.tpcb$date))
  temp <- readNWISdv(sitePorN1, paramtemp,
                       min(por.tpcb$date), max(por.tpcb$date))
  # Add USGS data to por.tpcb, matching dates
  por.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(por.tpcb$date, flow.1$Date)]
  por.tpcb$flow.2 <- 0.03*flow.2$X_00060_00003[match(por.tpcb$date, flow.2$Date)]
  por.tpcb$temp <- 273.15 + temp$X_00010_00003[match(por.tpcb$date, temp$Date)]
  # Remove samples with temp = NA
  por.tpcb.1 <- na.omit(por.tpcb)
}

# Random Forest Model -----------------------------------------------------
# Remove columns not used here
# Use flow.2 and DistanceSource1
por.tpcb.1 <- select(por.tpcb.1, -c(date, flow.1))

# Set seed for reproducibility
set.seed(123)

# Train-Test Split
train_indices <- sample(1:nrow(por.tpcb.1), 0.8 * nrow(por.tpcb.1))
train_data <- por.tpcb.1[train_indices, ]
test_data <- por.tpcb.1[-train_indices, ]

# Define hyperparameter grid
param_grid <- expand.grid(
  mtry = seq(1, ncol(train_data) - 1),  # Adjust mtry values based on your data
  splitrule = c("gini", "extratrees"),
  min.node.size = c(3, 4, 5)
)

# Prepare training control
ctrl <- trainControl(method = "cv", number = 5, search = "grid")

# Perform grid search with cross-validation using ranger
# Better fit with flow.2
rf_model <- train(
  log10(tPCB) ~ time + SiteID + season + flow.2 + temp + DistanceToNorthernLocation,
  data = train_data,
  method = "ranger",
  importance = 'permutation',
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best mtry
best_mtry <- rf_model$bestTune$mtry

final_rf_model <- ranger(
  formula = log10(tPCB) ~ time + SiteID + season + flow.2 + temp +
    DistanceToNorthernLocation,
  data = train_data,
  num.trees = 5000, # need to manually modify this parameter
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
          file = "Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFtPCB.csv",
          row.names = FALSE)

# Create a data frame for plotting and exporting
plot_data <- data.frame(
  Location = rep("Portland Harbor", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions
)

# Export results
write.csv(plot_data,
          file = "Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFObsPredtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(plot_data, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 1, fill = "white") +
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
ggsave("Output/Plots/Sites/ObsPred/PortlandHarbor/PortlandHarborRFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Individual PCB Analysis -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  por.pcb <- subset(por, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  por.pcb <- subset(por.pcb, select = -c(A1016:DistanceToNorthernLocation))
  # Log10 individual PCBs 
  por.pcb <- log10(por.pcb)
  # Replace -inf to NA
  por.pcb <- do.call(data.frame,
                     lapply(por.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  por.pcb.1 <- por.pcb[,
                       -which(colSums(is.na(por.pcb))/nrow(por.pcb) > 0.7)]
  # Change date format
  SampleDate <- as.Date(por$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(SampleDate),
                                  min(as.Date(SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(por$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add distance to northern location sample
  DistanceToNorthernLocation <- por$DistanceToNorthernLocation
  # Add date and time to por.pcb.1
  por.pcb.1 <- cbind(por.pcb.1, as.factor(por$SiteID), SampleDate,
                     data.frame(time.day), season.s, DistanceToNorthernLocation)
  # Include flow data from USGS station Portland Harbor
  sitePorN1 <- "14211720" # WILLAMETTE RIVER AT PORTLAND, OR
  sitePorN2 <- "14211820" # COLUMBIA SLOUGH AT PORTLAND, OR
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  # Flow (ft3/s)
  flow.2 <- readNWISdv(sitePorN2, paramflow,
                       min(por.pcb.1$SampleDate), max(por.pcb.1$SampleDate))
  temp <- readNWISdv(sitePorN1, paramtemp,
                     min(por.pcb.1$SampleDate), max(por.pcb.1$SampleDate))
  # Add USGS data to por.tpcb, matching dates
  por.pcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(por.pcb.1$SampleDate, flow.2$Date)]
  por.pcb.1$temp <- 273.15 + temp$X_00010_00003[match(por.pcb.1$SampleDate, temp$Date)]
  # Remove samples with temp = NA
  por.pcb.2 <- por.pcb.1[!is.na(por.pcb.1$temp), ]
  # Remove metadata not use in the random forest
  por.pcb.2 <- por.pcb.2[, !(names(por.pcb.2) %in% c("SampleDate"))]
}

# Set the seed for reproducibility
set.seed(123)

# Find the numeric columns (columns starting with "PCB")
pcb_numeric_columns <- grep("^PCB", colnames(por.pcb.2), value = TRUE)

# Find the corresponding character columns
char_columns <- setdiff(colnames(por.pcb.2), pcb_numeric_columns)

# Initialize the results matrix
rf_results <- data.frame(
  Congener = pcb_numeric_columns,
  RMSE = rep(NA, length(pcb_numeric_columns)),
  R_squared = rep(NA, length(pcb_numeric_columns)),
  Factor2_Percentage = rep(NA, length(pcb_numeric_columns))
)

# Create an empty data frame to store all predicted and actual data
all_results <- data.frame()

# Define parameter grid. More values can be included.
num_trees_grid <- c(50, 100, 200)
mtry_grid <- c(2, 3)
min_node_size_grid <- c(3, 4, 5)

# Initialize variables to store best parameters and performance
best_params <- c(Inf, Inf, Inf)  # Initial best performance (lower is better)
best_performance <- c(Inf, -Inf, -Inf)  # Initial best performance (higher is better)

# Perform grid search
for (num_trees in num_trees_grid) {
  for (mtry in mtry_grid) {
    for (min_node_size in min_node_size_grid) {
      
      # Initialize performance metrics
      avg_mse <- 0
      avg_r_squared <- 0
      avg_factor2_percentage <- 0
      
      # Iterate over each numeric column
      for (i in seq_along(pcb_numeric_columns)) {
        # Combine numeric and character data
        combined_data <- cbind(por.pcb.2[, pcb_numeric_columns[i], drop = FALSE],
                               por.pcb.2[, char_columns])
        
        # Exclude rows with missing values
        combined_data <- na.omit(combined_data)
        
        # Sample indices for training
        train_indices <- sample(1:nrow(combined_data), 0.8 * nrow(combined_data))
        
        # Create separate training and testing sets
        train_data <- combined_data[train_indices, ]
        test_data <- combined_data[-train_indices, ]
        
        # Train the ranger model with specified hyperparameters
        ranger_model <- ranger(
          dependent.variable.name = pcb_numeric_columns[i],
          data = train_data,
          num.trees = num_trees,
          mtry = mtry,
          min.node.size = min_node_size,
          seed = 123
        )
        
        # Predict on the test set
        predictions <- predict(ranger_model, data = test_data)$predictions
        
        # Calculate evaluation metrics
        mse <- mean((predictions - test_data[, pcb_numeric_columns[i]])^2)
        r_squared <- 1 - sum((test_data[, pcb_numeric_columns[i]] - predictions)^2) / sum((test_data[, pcb_numeric_columns[i]] - mean(test_data[, pcb_numeric_columns[i]]))^2)
        compare_df <- data.frame(observed = test_data[, pcb_numeric_columns[i]], predicted = predictions)
        compare_df$factor2 <- compare_df$observed / compare_df$predicted
        factor2_percentage <- sum(compare_df$factor2 > 0.5 & compare_df$factor2 < 2) / nrow(compare_df) * 100
        
        # Update average performance metrics
        avg_mse <- avg_mse + mse
        avg_r_squared <- avg_r_squared + r_squared
        avg_factor2_percentage <- avg_factor2_percentage + factor2_percentage
        
        # Append to the all_results dataframe
        col_results <- data.frame(
          Location = rep("Portland Harbor", nrow(test_data)),
          Congener = rep(pcb_numeric_columns[i], nrow(test_data)),
          Actual = test_data[, pcb_numeric_columns[i]],
          Predicted = predictions,
          R_squared = r_squared
        )
        all_results <- rbind(all_results, col_results)
        
        # Update rf_results with the average performance metrics for the current Congener
        rf_results$RMSE[i] <- sqrt(mse)
        rf_results$R_squared[i] <- r_squared
        rf_results$Factor2_Percentage[i] <- factor2_percentage
      }
      
      # Average performance metrics across all numeric columns
      avg_mse <- avg_mse / length(pcb_numeric_columns)
      avg_r_squared <- avg_r_squared / length(pcb_numeric_columns)
      avg_factor2_percentage <- avg_factor2_percentage / length(pcb_numeric_columns)
      
      # Update best parameters and performance if better
      if (avg_mse < best_performance[1] && avg_r_squared > best_performance[2] && avg_factor2_percentage > best_performance[3]) {
        best_params <- c(num_trees, mtry, min_node_size)
        best_performance <- c(avg_mse, avg_r_squared, avg_factor2_percentage)
      }
    }
  }
}

# Output best parameters and performance
print("Best Parameters:")
print(best_params)

# Remove congeners w/R2 < 0
rf_results <- rf_results %>%
  filter(R_squared >= 0)

# Add location name
rf_results <- cbind(Location = rep("Portland Harbor", nrow(rf_results)),
                    rf_results)

# Remove rows in all_results where R_squared < 0
all_results <- all_results %>%
  filter(R_squared >= 0)

# Remove the "R_squared" column from all_results
all_results <- all_results %>% select(-R_squared)

# Export results
write.csv(rf_results,
          file = "Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFPCB.csv",
          row.names = FALSE)

# Export combined results
write.csv(all_results,
          file = "Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFObsPredPCB.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(all_results, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 1, fill = "white") +
  scale_y_log10(limits = c(0.0001, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.0001, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration PCBi (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Print the plot
print(plotRFPCBi)

# Save plot in folder
ggsave("Output/Plots/Sites/ObsPred/PortlandHarbor/PortlandHarborRFPCB.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)

