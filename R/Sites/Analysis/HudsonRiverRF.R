## Water PCB concentrations data analysis per site
## Hudson River

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

# Select Hudson River data ---------------------------------------------------
hud <- wdc[str_detect(wdc$LocationName, 'Hudson River'),]
# PCBs were discharged to the river from the General Electric
# (GE) manufacturing plants in Hudson Falls and Fort Edward, NY
# Dredging from 2009 to 2015
# https://www.epa.gov/system/files/documents/2021-08/hudson_summer2021_floodplainrifs_factsheet_final.pdf

# Add distance to source --------------------------------------------------
{
  # Define source coordinates
  source1 <- c(Latitude = 43.295369, Longitude = -73.590631)  # GE Hudson Falls Plant
  source2 <- c(Latitude = 43.28639, Longitude = -73.588380)  # GE Fort Edward Plant
  # Create an sf point for the source
  source1_sf <- st_sfc(st_point(c(source1["Longitude"], source1["Latitude"])))
  source2_sf <- st_sfc(st_point(c(source2["Longitude"], source1["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(source1_sf) <- 4326
  st_crs(source2_sf) <- 4326
  # Transform source1_sf to UTM Zone 18N (EPSG:32618)
  source1_sf_utm <- st_transform(source1_sf, 32618)
  # Transform source2_sf to UTM Zone 18N (EPSG:32618)
  source2_sf_utm <- st_transform(source2_sf, 32618)
  # Convert the data frame to an sf object
  sf_hud <- st_as_sf(hud, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_hud <- st_set_crs(sf_hud, 4326)
  # Transform to UTM Zone 18N
  sf_hud_utm <- st_transform(sf_hud, 32618)
  # Calculate distances in meters from each location to source1
  distances_meters1 <- st_distance(sf_hud_utm, source1_sf_utm)
  # Convert distances to kilometers
  distances_km1 <- units::set_units(distances_meters1, "km")
  # Extract numeric values and assign to the DistanceToSource column
  hud$DistanceToSource1 <- as.numeric(distances_km1[, 1])
  # Calculate distances in meters from each location to source2
  distances_meters2 <- st_distance(sf_hud_utm, source2_sf_utm)
  # Convert distances to kilometers
  distances_km2 <- units::set_units(distances_meters2, "km")
  # Extract numeric values and assign to the DistanceToSource column
  hud$DistanceToSource2 <- as.numeric(distances_km2[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  hud$SampleDate <- as.Date(hud$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(hud$SampleDate),
                                  min(as.Date(hud$SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hud$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  hud.tpcb <- cbind(factor(hud$SiteID), hud$SampleDate, as.matrix(hud$tPCB),
                    data.frame(time.day), season.s, hud$DistanceToSource1,
                    hud$DistanceToSource2)
  # Add column names
  colnames(hud.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                          "DistanceSource1", "DistanceSource2")
}

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

# Random Forest Model tPCB ------------------------------------------------
# Remove columns not used here
# Use flow.2 and DistanceSource1
hud.tpcb.2 <- select(hud.tpcb.2, -c(date, flow.1, flow.3, flow.4,
                                    DistanceSource2))

# Set seed for reproducibility
set.seed(123)

# Train-Test Split
train_indices <- sample(1:nrow(hud.tpcb.2), 0.8 * nrow(hud.tpcb.2))
train_data <- hud.tpcb.2[train_indices, ]
test_data <- hud.tpcb.2[-train_indices, ]

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
  log10(tPCB) ~ time + SiteID + season + flow.2 + temp + DistanceSource1,
  data = train_data,
  method = "ranger",
  importance = 'permutation',
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best mtry
best_mtry <- rf_model$bestTune$mtry

final_rf_model <- ranger(
  formula = log10(tPCB) ~ time + SiteID + season + flow.2 + temp + DistanceSource1,
  data = train_data,
  num.trees = 5000, # need to manualy modify this parameter
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
          file = "Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFtPCB.csv",
          row.names = FALSE)

# Create a data frame for plotting
plot_data <- data.frame(
  Location = rep("Hudson River", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions
)

# Export results
write.csv(plot_data,
          file = "Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFObsPredtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(plot_data, aes(x = 10^(Actual), y = 10^(Predicted))) +
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
ggsave("Output/Plots/Sites/ObsPred/HudsonRiver/HudsonRiverRFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Random Forest Model individual PCBs -----------------------------------
{
  hud.pcb <- subset(hud, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  hud.pcb <- subset(hud.pcb, select = -c(A1016:DistanceToSource2))
  # Log10 individual PCBs 
  hud.pcb <- log10(hud.pcb)
  # Replace -inf to NA
  hud.pcb <- do.call(data.frame,
                     lapply(hud.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  hud.pcb.1 <- hud.pcb[,
                       -which(colSums(is.na(hud.pcb))/nrow(hud.pcb) > 0.7)]
  # Change date format
  SampleDate <- as.Date(hud$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(SampleDate),
                                  min(as.Date(SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hud$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add distance to source
  DistanceToSource1 <- hud$DistanceToSource1
  # Add date and time to hud.pcb.1
  hud.pcb.1 <- cbind(hud.pcb.1, as.factor(hud$SiteID), SampleDate,
                     data.frame(time.day), season.s, DistanceToSource1)
  # Remove site Bakers Falls. Upstream source
  # North Bakers Falls = WCPCB-HUD006 and
  # South Bakers Falls = WCPCB-HUD010.
  # Remove rows with SiteID equal to "WCPCB-HUD006" or "WCPCB-HUD010"
  hud.pcb.1 <- hud.pcb.1[!(hud.pcb.1$`as.factor(hud$SiteID)` %in% c("WCPCB-HUD006",
                                                              "WCPCB-HUD010")), ]
  # Include flow data from USGS station Hudson River
  # sitehudN3 for flow and sitehudN5 for water temperature
  sitehudN3 <- "01328770" # HUDSON RIVER AT THOMSON NY, no temp!
  sitehudN5 <- "01359139" # HUDSON RIVER AT ALBANY NY No flow!
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  # Retrieve USGS data
  # Flow (ft3/s)
  flow.2 <- readNWISdv(sitehudN2, paramflow,
                       min(hud.pcb.1$SampleDate), max(hud.pcb.1$SampleDate))
  # Water temperature in Celsius
  temp <- readNWISdv(sitehudN5, paramtemp,
                     min(hud.pcb.1$SampleDate), max(hud.pcb.1$SampleDate))
  
  # Add USGS data to hud.tpcb.1 matching dates
  hud.pcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(hud.pcb.1$SampleDate,
                                                      flow.2$Date)]
  hud.pcb.1$temp <- 273.15 + temp$X_00010_00003[match(hud.pcb.1$SampleDate,
                                                      temp$Date)]
  # Remove samples with flow.2 = NA
  hud.pcb.2 <- hud.pcb.1[!is.na(hud.pcb.1$flow.2), ]
  # Remove metadata not use in the random forest
  hud.pcb.2 <- hud.pcb.2[, !(names(hud.pcb.2) %in% c("SampleDate"))]
}

# Set seed for reproducibility
set.seed(123)

# Find the numeric columns (columns starting with "PCB")
pcb_numeric_columns <- grep("^PCB", colnames(hud.pcb.2), value = TRUE)

# Find the corresponding character columns
char_columns <- setdiff(colnames(hud.pcb.2), pcb_numeric_columns)

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
num_trees_grid <- c(1000, 4000, 5000)
mtry_grid <- c(4, 5, 6)
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
        combined_data <- cbind(hud.pcb.2[, pcb_numeric_columns[i], drop = FALSE],
                               hud.pcb.2[, char_columns])
        
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
          Location = rep("USA", nrow(test_data)),
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
rf_results <- cbind(Location = rep("Hudson River", nrow(rf_results)),
                    rf_results)

# Remove rows in all_results where R_squared < 0
all_results <- all_results %>%
  filter(R_squared >= 0)

# Remove the "R_squared" column from all_results
all_results <- all_results %>% select(-R_squared)

# Export results
write.csv(rf_results,
          file = "Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFPCB.csv",
          row.names = FALSE)

# Export combined results
write.csv(all_results,
          file = "Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFObsPredPCB.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(all_results, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 1, fill = "white") +
  scale_y_log10(limits = c(0.1, 10^5),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.1, 10^5),
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
ggsave("Output/Plots/Sites/ObsPred/HudsonRiver/HudsonRiverRFPCB.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)


