## Water PCB concentrations data analysis per site
## Blue River
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
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# Select Bannister Federal Complex data ---------------------------------------------------
bfc <- wdc[str_detect(wdc$LocationName, 'Bannister Fed Complex'),]

# Calculate central location ---------------------------------------------
{
  # Calculate the mean latitude and longitude
  center_lat <- mean(bfc$Latitude)
  center_lon <- mean(bfc$Longitude)
  # Create a data frame for the center
  center_df <- data.frame( SiteID = "Center", Latitude = center_lat,
                           Longitude = center_lon)
  # Convert the center data frame to an sf object
  sf_center <- st_as_sf(center_df, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326) for the center
  sf_center <- st_set_crs(sf_center, 4326)
  # Transform to UTM Zone 15 for the center
  sf_center_utm <- st_transform(sf_center, 32615)
  # Convert the data frame to an sf object
  sf_bfc <- st_as_sf(bfc, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_bfc <- st_set_crs(sf_bfc, 4326)
  # Transform to UTM Zone 15
  sf_bfc <- st_transform(sf_bfc, 32615)
  # Calculate distances in meters from each location to the center
  distances_meters <- st_distance(sf_bfc, sf_center_utm)
  # Convert distances to kilometers
  distances_km <- units::set_units(distances_meters, "km")
  # Extract numeric values and assign to the DistanceToCentroid column
  bfc$DistanceToCentroid <- as.numeric(distances_km[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  bfc$SampleDate <- as.Date(bfc$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(bfc$SampleDate),
                                  min(as.Date(bfc$SampleDate)), units = "days"))
  # Create individual code for each site sampled
  site.numb <- bfc$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(bfc$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  bfc.tpcb <- cbind(factor(bfc$SiteID), bfc$SampleDate, as.matrix(bfc$tPCB),
                    data.frame(time.day), site.numb, season.s,
                    bfc$DistanceToCentroid)
  # Add column names
  colnames(bfc.tpcb) <- c("SiteID", "date", "tPCB", "time", "site.code",
                          "season", "DistanceToCentroid")
}

# Random Forest Model tPCB ------------------------------------------------
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(bfc.tpcb), 0.8 * nrow(bfc.tpcb))
train_data <- bfc.tpcb[train_indices, ]
test_data <- bfc.tpcb[-train_indices, ]

# Fit the GBM Model
gbm_model <- gbm(log10(tPCB) ~ time + SiteID + season + DistanceToCentroid,
                 data = train_data,
                 distribution = "gaussian",  # For regression tasks
                 n.trees = 5000,             # Number of trees
                 interaction.depth = 5,      # Interaction depth
                 shrinkage = 0.003)           # Learning rate

# Print summary of the model
summary(gbm_model)

# Make predictions
predictions <- predict(
  object = gbm_model,
  newdata = test_data,
  n.trees = 5000)

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

# Create the data frame directly
performance_RF <- data.frame(Heading = c("RMSE", "R2", "Factor2"),
                             Value = c(rmse, r_squared,
                                       factor2_percentage))
# Print the modified data frame
print(performance_RF)

# Export results
write.csv(performance_RF,
          file = "Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexRFtPCBV02.csv",
          row.names = FALSE)

# Create a data frame for plotting Observations vs Predictions
plot_data <- data.frame(
  Location = rep("Bannister Fed Complex", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions
)

# Export results
write.csv(plot_data,
          file = "Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexRFObsPredtPCBV02.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(plot_data, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10, 10^6),
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

# Print the plot
print(plotRF)

# Save plot in folder
ggsave("Output/Plots/Sites/ObsPred/BannisterFedComplex/BannisterFedComplexRFtPCBV02.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Random Forest Model individual PCBs -------------------------------------
# Prepare data.frame
{
  bfc.pcb <- subset(bfc, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  bfc.pcb <- subset(bfc.pcb, select = -c(A1016:DistanceToCentroid))
  # Log10 individual PCBs 
  bfc.pcb <- log10(bfc.pcb)
  # Replace -inf to NA
  bfc.pcb <- do.call(data.frame,
                     lapply(bfc.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  bfc.pcb.1 <- bfc.pcb[,
                       -which(colSums(is.na(bfc.pcb))/nrow(bfc.pcb) > 0.7)]
  
  # Create individual code for each site sampled
  site.numb <- factor(bfc$SiteID %>% as.factor() %>% as.numeric)
  # Change date format
  SampleDate <- as.Date(bfc$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(SampleDate, min(SampleDate), units = "days"))
  # Add distance to the centroid
  centroid <- bfc$DistanceToCentroid
  # Include season
  yq.s <- as.yearqtr(as.yearmon(bfc$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to bfc.pcb.1
  bfc.pcb.1 <- cbind(bfc.pcb.1, data.frame(time.day), site.numb,
                     season.s, centroid)
}

# Perform imputation of missing values with the adjusted lowest observed value
bfc.pcb.1_imputed <- bfc.pcb.1

# Iterate over each column (try one of the options above)
for (col in colnames(bfc.pcb.1_imputed)) {
  # Check if the column is numeric
  if (is.numeric(bfc.pcb.1_imputed[[col]])) {
    # Find the lowest observed value in the column
    lowest_value <- min(bfc.pcb.1_imputed[[col]], na.rm = TRUE)
    # Calculate the adjusted
    # Option 1: Lowest value
    # adjusted_value <- lowest_value
    # Option 2: Lowest value/sqrt(2)
    # adjusted_value <- lowest_value / sqrt(2)
    # Option 3: Mean Imputation
    adjusted_value <- mean(bfc.pcb.1_imputed[[col]], na.rm = TRUE)
    # Option 4: Median Imputation
    # adjusted_value <- median(bfc.pcb.1_imputed[[col]], na.rm = TRUE)
    # Option 5: Robust Scaling (using median and interquartile range)
    # Calculate the adjusted values using robust scaling
    # non_na_values <- bfc.pcb.1_imputed[[col]][!is.na(bfc.pcb.1_imputed[[col]])]
    # median_value <- median(non_na_values)
    # iqr_value <- IQR(non_na_values)
    # Calculate adjusted values for non-NA values
    # adjusted_values <- (non_na_values - median_value) / iqr_value
    # Replace missing values with the adjusted lowest value
    bfc.pcb.1_imputed[[col]][is.na(bfc.pcb.1_imputed[[col]])] <- adjusted_value
  }
}

# Set the seed for reproducibility
set.seed(123)

# Find the numeric columns (columns starting with "PCB")
pcb_numeric_columns <- grep("^PCB", colnames(bfc.pcb.1_imputed), value = TRUE)

# Find the corresponding character columns
char_columns <- setdiff(colnames(bfc.pcb.1_imputed), pcb_numeric_columns)

# Initialize the results matrix
rf_results <- data.frame(
  Congener = pcb_numeric_columns,
  RMSE = rep(NA, length(pcb_numeric_columns)),
  R_squared = rep(NA, length(pcb_numeric_columns)),
  Factor2_Percentage = rep(NA, length(pcb_numeric_columns))
)

# Create an empty data frame to store all predicted and actual data
all_results <- data.frame()

# Iterate over each numeric column
for (i in seq_along(pcb_numeric_columns)) {
  # Combine numeric and character data
  combined_data <- cbind(bfc.pcb.1_imputed[, pcb_numeric_columns[i], drop = FALSE],
                         bfc.pcb.1_imputed[, char_columns])
  
  # Convert combined_data to data frame
  combined_data <- as.data.frame(combined_data)
  
  # Sample indices for training
  train_indices <- sample(1:nrow(combined_data), 0.8 * nrow(combined_data))
  
  # Create separate training and testing sets
  train_data <- combined_data[train_indices, ]
  test_data <- combined_data[-train_indices, ]
  
  # Modeling code using gbm
  fit <- gbm(train_data[, 1] ~ ., data = train_data[, -1], 
             distribution = "gaussian",  # For regression tasks
             n.trees = 5000,             # Number of trees
             interaction.depth = 5,      # Interaction depth
             shrinkage = 0.003)          # Learning rate
  
  # Example: Make predictions on the test set
  predictions <- predict(
    object = fit,
    newdata = test_data,
    n.trees = 5000)
  
  # Calculate mean squared error (mse) for illustration
  mse <- mean((predictions - test_data[, 1])^2)
  
  # Calculate R-squared
  r_squared <- 1 - (sum((test_data[, 1] - predictions)^2) / sum((test_data[, 1] - mean(test_data[, 1]))^2))
  
  # Calculate factor2_percentage within the loop
  compare_df <- data.frame(
    observed = test_data[, 1],
    predicted = predictions
  )
  compare_df$factor2 <- compare_df$observed / compare_df$predicted
  factor2_percentage <- sum(compare_df$factor2 > 0.5 & compare_df$factor2 < 2) / nrow(compare_df) * 100
  
  # Store the results in the matrix
  rf_results[i, 2:4] <- c(mse, r_squared, factor2_percentage)
  
  # Create a data frame for each column's results
  col_results <- data.frame(
    Location = rep("Bannister Fed Complex", length(test_data[, 1])),
    Congener = rep(pcb_numeric_columns[i], length(test_data[, 1])),
    Actual = test_data[, 1],
    Predicted = predictions,
    R_squared = r_squared  # Add R_squared column
  )
  
  # Bind the data frame to the overall results
  all_results <- rbind(all_results, col_results)
}

# Remove congeners w/R2 < 0
rf_results <- rf_results %>%
  filter(R_squared >= 0)

# Add location name
rf_results <- cbind(Location = rep("Bannister Fed Complex", nrow(rf_results)),
                    rf_results)

# Remove rows in all_results where R_squared < 0
all_results <- all_results %>%
  filter(R_squared >= 0)

# Remove the "R_squared" column from all_results
all_results <- all_results %>% select(-R_squared)

# Export results
write.csv(rf_results,
          file = "Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexRFPCBV02.csv",
          row.names = FALSE)

# Export combined results
write.csv(all_results,
          file = "Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexRFObsPredPCBV02.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(all_results, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(1, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration PCBi (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) + # 1:1 line
  geom_abline(intercept = -1, slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = 1, slope = 1, col = "blue", linewidth = 0.7) +   # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Print the plot
print(plotRFPCBi)

# Save plot in folder
ggsave("Output/Plots/Sites/ObsPred/BannisterFedComplex/BannisterFedComplexRFPCBV02.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)

