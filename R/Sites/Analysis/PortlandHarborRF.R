## Water PCB concentrations data analysis per site
## Portland Harbor
## Random Forest model

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
  time.day <- data.frame(as.Date(por$SampleDate) - min(as.Date(por$SampleDate)))
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
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(por.tpcb.1), 0.8 * nrow(por.tpcb.1))
train_data <- por.tpcb.1[train_indices, ]
test_data <- por.tpcb.1[-train_indices, ]

# Fit the Model. Flow.2
rf_model.1 <- randomForest(log10(tPCB) ~ time + SiteID + season +
                             flow.2 + temp + DistanceToNorthernLocation,
                           data = train_data)

# Make Predictions
predictions.1 <- predict(rf_model.1, newdata = test_data)

# Evaluate Model Performance
mse.1 <- mean((predictions.1 - log10(test_data$tPCB))^2)
rmse.1 <- sqrt(mse.1)
r_squared.1 <- 1 - (sum((log10(test_data$tPCB) - predictions.1)^2)/sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2))

# Estimate a factor of 2 between observations and predictions
# Create a data frame with observed and predicted values
compare_df.1 <- data.frame(observed = test_data$tPCB,
                           predicted = 10^predictions.1)

# Estimate a factor of 2 between observations and predictions
compare_df.1$factor2 <- compare_df.1$observed/compare_df.1$predicted

# Calculate the percentage of observations within the factor of 2
factor2_percentage.1 <- nrow(compare_df.1[compare_df.1$factor2 > 0.5 & compare_df.1$factor2 < 2, ])/nrow(compare_df.1)*100

# Create the data frame directly
performance_df <- data.frame(Heading = c("RMSE", "R2", "Factor2"),
                             Value = c(rmse.1, r_squared.1,
                                       factor2_percentage.1))

# Remove unnecessary columns
performance_df <- performance_df[, !(names(performance_df) %in% c("V1", "V2", "V3"))]

# Print the modified data frame
print(performance_df)

# Export results
write.csv(performance_df,
          file = "Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFPerformancetPCBV02.csv")

# Feature Importance
importance.1 <- importance(rf_model.1)
barplot(importance.1[, 1], names.arg = rownames(importance.1),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting and exporting
plot_data.1 <- data.frame(
  Location = rep("Portland Harbor", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions.1
)

# Export results
write.csv(plot_data.1,
          file = "Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFObsPredtPCBV02.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(plot_data.1, aes(x = 10^(Actual), y = 10^(Predicted))) +
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
ggsave("Output/Plots/Sites/ObsPred/PortlandHarbor/PortlandHarborRFtPCBV02.png",
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
  # Add site ID
  SiteID <- factor(por$SiteID)
  # Change date format
  SampleDate <- as.Date(por$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(por$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add distance to northern location sample
  DistanceToNorthernLocation <- por$DistanceToNorthernLocation
  # Add date and time to por.pcb.1
  por.pcb.1 <- cbind(por.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     season.s, DistanceToNorthernLocation)
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

# Iterate over each numeric column
for (i in seq_along(pcb_numeric_columns)) {
  # Combine numeric and character data
  combined_data <- cbind(por.pcb.2[, pcb_numeric_columns[i],
                                   drop = FALSE], por.pcb.2[, char_columns])
  
  # Exclude rows with missing values
  combined_data <- na.omit(combined_data)
  
  # Sample indices for training
  train_indices <- sample(1:nrow(combined_data), 0.8 * nrow(combined_data))
  
  # Create separate training and testing sets
  train_data <- combined_data[train_indices, ]
  test_data <- combined_data[-train_indices, ]
  
  # Modeling code using randomForest
  fit <- randomForest(train_data[, 1] ~ ., data = train_data)
  
  # Example: Make predictions on the test set
  predictions <- predict(fit, newdata = test_data)
  
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
    Location = rep("Portland Harbor", length(test_data[, 1])),
    Congener = rep(pcb_numeric_columns[i], length(test_data[, 1])),
    Actual = test_data[, 1],
    Predicted = predictions
  )
  
  # Bind the data frame to the overall results
  all_results <- rbind(all_results, col_results)
}

# Export results
write.csv(rf_results,
          file = "Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFPerformancePCBV02.csv",
          row.names = FALSE)

# Export combined results
write.csv(all_results,
          file = "Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFObsPredPCBV02.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(all_results, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(0.0001, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.0001, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration PCBi (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = 0.30103, slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = -0.30103, slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Print the plot
print(plotRFPCBi)

# Save plot in folder
ggsave("Output/Plots/Sites/ObsPred/PortlandHarbor/PortlandHarborRFPCBV02.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)

