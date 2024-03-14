## Water PCB concentrations data analysis per site
## Chesapeake Bay & Delaware Canal

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
  library(sf) # Create file to be used in Google Earth
  library(units)
  library(ranger) # Random Forest functions
  library(caret) # For cross-validation
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# Select Chesapeake Bay & Delaware Canal data ---------------------------------------------------
che <- wdc[str_detect(wdc$LocationName, 'Chesapeake Bay'),]

# Create centroid and distance to locations -------------------------------
{
  # Calculate the mean latitude and longitude
  center_lat <- mean(che$Latitude)
  center_lon <- mean(che$Longitude)
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
  # Transform to UTM Zone 18 for the center
  sf_center_utm <- st_transform(sf_center, 32618)
  # Convert the data frame to an sf object
  sf_che <- st_as_sf(che, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_che <- st_set_crs(sf_che, 4326)
  # Transform to UTM Zone 18
  sf_che <- st_transform(sf_che, 32618)
  # Calculate distances in meters from each location to the center
  distances_meters <- st_distance(sf_che, sf_center_utm)
  # Convert distances to kilometers
  distances_km <- units::set_units(distances_meters, "km")
  # Extract numeric values and assign to the DistanceToCentroid column
  che$DistanceToCentroid <- as.numeric(distances_km[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  che$SampleDate <- as.Date(che$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(che$SampleDate),
                                  min(as.Date(che$SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(che$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  che.tpcb <- cbind(factor(che$SiteID), che$SampleDate, as.matrix(che$tPCB),
                    data.frame(time.day), season.s, che$DistanceToCentroid)
  # Add column names
  colnames(che.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                          "DistanceToCentroid")
}

# Add water temperature data ----------------------------------------------
# See code: R/ExtractingData/ChesapeakeBay/WaterTemp.R
{
  # Read water temperature
  wtp <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/Watertemp/ChesapeakeBayWT.csv")
  # Convert date columns to Date format
  wtp$Date <- as.Date(wtp$Date)
  # Add water temperature to grl.tpcb
  che.tpcb$temp <- wtp$WTMP_K[match(che.tpcb$date, wtp$Date)]
  # Remove samples with temp = NA
  che.tpcb <- na.omit(che.tpcb)
}

# Remove date column
che.tpcb <- select(che.tpcb, -c(date))

# Random Forest Model tPCB ------------------------------------------------
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(che.tpcb), 0.8 * nrow(che.tpcb))
train_data <- che.tpcb[train_indices, ]
test_data <- che.tpcb[-train_indices, ]

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
  log10(tPCB) ~ time + SiteID + season + temp + DistanceToCentroid,
  data = train_data,
  method = "ranger",
  importance = 'permutation',
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best mtry
best_mtry <- rf_model$bestTune$mtry

final_rf_model <- ranger(
  formula = log10(tPCB) ~ time + SiteID + season + temp +
    DistanceToCentroid, data = train_data,
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
compare_df <- data.frame(observed = test_data$tPCB,
                           predicted = 10^predictions)

# Estimate a factor of 2 between observations and predictions
compare_df$factor2 <- compare_df$observed/compare_df$predicted

# Calculate the percentage of observations within the factor of 2
factor2_percentage <- nrow(compare_df[compare_df$factor2 > 0.5 & compare_df$factor2 < 2, ])/nrow(compare_df)*100

# Create the data frame directly
performance_df <- data.frame(Heading = c("RMSE", "R2", "Factor2"),
                             Value = c(rmse, r_squared,
                                       factor2_percentage))

# Remove unnecessary columns
performance_df <- performance_df[, !(names(performance_df) %in% c("V1", "V2", "V3"))]

# Print the modified data frame
print(performance_df)

# Export results
write.csv(performance_df,
          file = "Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFtPCB.csv",
          row.names = FALSE)

# Create a data frame for plotting
plot_data <- data.frame(
  Location = rep("Chesapeake Bay", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions
)

# Export results
write.csv(plot_data,
          file = "Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFObsPredtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(plot_data, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(100, 10^5),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(100, 10^5),
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
ggsave("Output/Plots/Sites/ObsPred/ChesapeakeBay/ChesapeakeBayRFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Random Forest Model individual PCBs -------------------------------------
# Prepare data.frame
{
  che.pcb <- subset(che, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  che.pcb <- subset(che.pcb, select = -c(A1016:DistanceToCentroid))
  # Log10 individual PCBs 
  che.pcb <- log10(che.pcb)
  # Replace -inf to NA
  che.pcb <- do.call(data.frame,
                     lapply(che.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  che.pcb.1 <- che.pcb[,
                       -which(colSums(is.na(che.pcb))/nrow(che.pcb) > 0.7)]
  # Change date format
  SampleDate <- as.Date(che$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(SampleDate),
                                  min(as.Date(SampleDate)), units = "days"))
  # Add distance to the centroid
  centroid <- che$DistanceToCentroid
  # Include season
  yq.s <- as.yearqtr(as.yearmon(che$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to grl.pcb.1
  che.pcb.1 <- cbind(che.pcb.1, as.factor(che$SiteID), SampleDate,
                     data.frame(time.day), season.s, centroid)
  # Add water temperature to grl.pcb.1
  che.pcb.1$temp <- wtp$WTMP_K[match(che.pcb.1$SampleDate,
                                               wtp$Date)]
  # Remove samples with temperature = NA
  che.pcb.2 <- che.pcb.1[!is.na(che.pcb.1$temp), ]
  # Remove metadata not use in the random forest
  che.pcb.2 <- che.pcb.2[, !(names(che.pcb.2) %in% c("SampleDate"))]
}

# Set seed for reproducibility
set.seed(123)

# Find the numeric columns (columns starting with "PCB")
pcb_numeric_columns <- grep("^PCB", colnames(che.pcb.2), value = TRUE)

# Find the corresponding character columns
char_columns <- setdiff(colnames(che.pcb.2), pcb_numeric_columns)

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
mtry_grid <- c(3, 4)
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
        combined_data <- cbind(che.pcb.2[, pcb_numeric_columns[i], drop = FALSE],
                               che.pcb.2[, char_columns])
        
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
          Location = rep("Chesapeake Bay", nrow(test_data)),
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
rf_results <- cbind(Location = rep("Chesapeake Bay", nrow(rf_results)),
                    rf_results)

# Remove rows in all_results where R_squared < 0
all_results <- all_results %>%
  filter(R_squared >= 0)

# Remove the "R_squared" column from all_results
all_results <- all_results %>% select(-R_squared)

# Export results
write.csv(rf_results,
          file = "Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFPCB.csv",
          row.names = FALSE)

# Export combined results
write.csv(all_results,
          file = "Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFObsPredPCB.csv",
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
ggsave("Output/Plots/Sites/ObsPred/ChesapeakeBay/ChesapeakeBayRFPCB.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)
