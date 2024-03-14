## Water PCB concentrations data analysis per site
## Passaic River
## Source: Diamond Alkali. 80-120 Lister Avenue, Newark, NJ.

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

# Select Passaic River data ---------------------------------------------------
pas <- wdc[str_detect(wdc$LocationName, 'Passaic River'),]

# Add distance to source --------------------------------------------------
{
  # Define source coordinates
  source <- c(Latitude = 40.740238, Longitude = -74.133805)  # Former Diamond Alkali
  # Create an sf point for the source
  source_sf <- st_sfc(st_point(c(source["Longitude"], source["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(source_sf) <- 4326
  # Transform source1_sf to UTM Zone 18N (EPSG:32618)
  source_sf_utm <- st_transform(source_sf, 32618)
  # Convert the data frame to an sf object
  sf_pas <- st_as_sf(pas, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_pas <- st_set_crs(sf_pas, 4326)
  # Transform to UTM Zone 18N
  sf_pas_utm <- st_transform(sf_pas, 32618)
  # Calculate distances in meters from each location to source
  distances_meters <- st_distance(sf_pas_utm, source_sf_utm)
  # Convert distances to kilometers
  distances_km <- units::set_units(distances_meters, "km")
  # Extract numeric values and assign to the DistanceToSource column
  pas$DistanceToSource <- as.numeric(distances_km[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  pas$SampleDate <- as.Date(pas$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(pas$SampleDate),
                                  min(as.Date(pas$SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(pas$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  pas.tpcb <- cbind(factor(pas$SiteID), pas$SampleDate, as.matrix(pas$tPCB),
                    data.frame(time.day), season.s, pas$DistanceToSource)
  # Add column names
  colnames(pas.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                           "DistanceToSource")
}

# Include USGS flow and temperature data --------------------------------------------------
{
  # Include flow data from USGS station Passaic River
  sitepasN1 <- "01381900" # No temp
  sitepasN2 <- "01379500" # No temp
  sitepasN3 <- "01389005" # No flow
  sitepasN4 <- "01389010" # No temp
  sitepasN5 <- "01389500" # No temp
  sitepasN6 <- "01389890" # No temp

  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow.1 <- readNWISdv(sitepasN1, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  flow.2 <- readNWISdv(sitepasN2, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  flow.3 <- readNWISdv(sitepasN4, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  flow.4 <- readNWISdv(sitepasN5, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  flow.5 <- readNWISdv(sitepasN6, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  temp <- readNWISdv(sitepasN3, paramtemp,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  
  # Add USGS data to pass.tpcb.2, matching dates, conversion to m3/s
  pas.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.1$Date)]
  pas.tpcb$flow.2 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.2$Date)]
  pas.tpcb$flow.3 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.3$Date)]
  pas.tpcb$flow.4 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.4$Date)]
  pas.tpcb$flow.5 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.5$Date)]
  pas.tpcb$temp <- 273.15 + temp$X_.from.middle.intake_00010_00003[match(pas.tpcb$date,
                                                       temp$Date)]
}

# Remove site -------------------------------------------------------------
# Remove site located in the ocean. Possible typo in original coordinates.
pas.tpcb.1 <- subset(pas.tpcb, SiteID != c("WCPCB-PAS022"))

# Random Forest Model tPCB ------------------------------------------------
# Remove columns not used here
# Using flow.1
pas.tpcb.1 <- select(pas.tpcb.1, -c(date, flow.2, flow.3, flow.4, flow.5))
# Train-Test Split
set.seed(123)

# Train-test split
train_indices <- sample(1:nrow(pas.tpcb.1), 0.8 * nrow(pas.tpcb.1))
train_data <- pas.tpcb.1[train_indices, ]
test_data <- pas.tpcb.1[-train_indices, ]

# Define hyperparameter grid
param_grid <- expand.grid(
  mtry = seq(1, ncol(train_data) - 1),  # Adjust mtry values based on your data
  splitrule = c("gini", "extratrees"),
  min.node.size = c(3, 4, 5)
)

# Prepare training control
ctrl <- trainControl(method = "cv", number = 5, search = "grid")

# Perform grid search with cross-validation using ranger
rf_model <- train(
  log10(tPCB) ~ time + SiteID + season + flow.1 + temp + DistanceToSource,
  data = train_data,
  method = "ranger",
  importance = 'permutation',
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best mtry
best_mtry <- rf_model$bestTune$mtry

final_rf_model <- ranger(
  formula = log10(tPCB) ~ time + SiteID + season + flow.1 + temp +
    DistanceToSource,
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
          file = "Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFtPCB.csv",
          row.names = FALSE)

# Create a data frame for plotting and exporting
plot_data <- data.frame(
  Location = rep("Passaic River", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions
)

# Export results
write.csv(plot_data,
          file = "Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFObsPredtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(plot_data, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 1, fill = "white") +
  scale_y_log10(limits = c(1, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^6),
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
ggsave("Output/Plots/Sites/ObsPred/PassaicRiver/PassaicRiverRFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Random Forest Model individual PCBs -------------------------------------
{
  # Remove metadata
  pas.pcb <- subset(pas, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  pas.pcb <- subset(pas.pcb, select = -c(A1016:DistanceToSource))
  # Log10 individual PCBs 
  pas.pcb <- log10(pas.pcb)
  # Replace -inf to NA
  pas.pcb <- do.call(data.frame,
                     lapply(pas.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  pas.pcb.1 <- pas.pcb[,
                       -which(colSums(is.na(pas.pcb))/nrow(pas.pcb) > 0.7)]
  # Change date format
  SampleDate <- as.Date(pas$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(SampleDate),
                                  min(as.Date(SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(pas$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add distance to source
  DistanceToSource <- pas$DistanceToSource
  # Add date and time to pas.pcb.1
  pas.pcb.1 <- cbind(pas.pcb.1, as.factor(pas$SiteID), SampleDate,
                      data.frame(time.day), season.s, DistanceToSource)
  # Remove site located in the ocean. Possible typo in original coordinates.
  pas.pcb.1 <- subset(pas.pcb.1, !(as.factor(pas$SiteID) %in% c("WCPCB-PAS022")))
  # Include flow data from USGS station Passaic River
  sitepasN1 <- "01381900" # No temp
  sitepasN3 <- "01389005" # No flow
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow.1 <- readNWISdv(sitepasN1, paramflow,
                       min(pas.pcb.1$SampleDate), max(pas.pcb.1$SampleDate))
  temp <- readNWISdv(sitepasN3, paramtemp,
                     min(pas.pcb.1$SampleDate), max(pas.pcb.1$SampleDate))
  # Add USGS data to pas.tpcb.1, matching dates, conversion to m3/s
  pas.pcb.1$flow.1 <- 0.03*flow.1$X_00060_00003[match(pas.pcb.1$SampleDate,
                                                      flow.1$Date)]
  pas.pcb.1$temp <- 273.15 + temp$X_.from.middle.intake_00010_00003[match(pas.pcb.1$SampleDate,
                                                                          temp$Date)]
  # Remove metadata not use in the random forest
  pas.pcb.1 <- pas.pcb.1[, !(names(pas.pcb.1) %in% c("SampleDate"))]
}

# Set the seed for reproducibility
set.seed(123)

# Find the numeric columns (columns starting with "PCB")
pcb_numeric_columns <- grep("^PCB", colnames(pas.pcb.1), value = TRUE)

# Find the corresponding character columns
char_columns <- setdiff(colnames(pas.pcb.1), pcb_numeric_columns)

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
num_trees_grid <- c(500, 750, 1000)
mtry_grid <- c(4, 5)
min_node_size_grid <- c(4, 5, 6)

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
        combined_data <- cbind(pas.pcb.1[, pcb_numeric_columns[i], drop = FALSE],
                               pas.pcb.1[, char_columns])
        
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
          Location = rep("Passaic River", nrow(test_data)),
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
rf_results <- cbind(Location = rep("Passaic River", nrow(rf_results)),
                    rf_results)

# Remove rows in all_results where R_squared < 0
all_results <- all_results %>%
  filter(R_squared >= 0)

# Remove the "R_squared" column from all_results
all_results <- all_results %>% select(-R_squared)

# Export results
write.csv(rf_results,
          file = "Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFPCB.csv",
          row.names = FALSE)

# Export combined results
write.csv(all_results,
          file = "Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFObsPredPCB.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(all_results, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 1, fill = "white") +
  scale_y_log10(limits = c(0.01, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.01, 10^6),
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
ggsave("Output/Plots/Sites/ObsPred/PassaicRiver/PassaicRiverRFPCB.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)


