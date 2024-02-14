## Water PCB concentrations data analysis per site
## Hudson River

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
  time.day <- data.frame(as.Date(hud$SampleDate) - min(as.Date(hud$SampleDate)))
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

# Random Forest Model -----------------------------------------------------
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(hud.tpcb.2), 0.8 * nrow(hud.tpcb.2))
train_data <- hud.tpcb.2[train_indices, ]
test_data <- hud.tpcb.2[-train_indices, ]

# Fit the Model (4)
rf_model.1 <- randomForest(log10(tPCB) ~ time + SiteID + season + flow.3 +
                             temp + DistanceSource1, data = train_data)

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
          file = "Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFPerformancetPCB.csv",
          row.names = FALSE)

# Feature Importance
importance.1 <- importance(rf_model.1)
barplot(importance.1[, 1], names.arg = rownames(importance.1),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting
plot_data.1 <- data.frame(
  Location = rep("Hudson River", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions.1
)

# Export results
write.csv(plot_data.1,
          file = "Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFObsPredtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(plot_data.1, aes(x = 10^(Actual), y = 10^(Predicted))) +
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

# Print the plot
print(plotRF)

# Save plot in folder
ggsave("Output/Plots/Sites/ObsPred/HudsonRiver/HudsonRiverRFtPCBV01.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Individual PCB Analysis -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
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
  # Add site ID
  SiteID <- factor(hud$SiteID)
  # Change date format
  SampleDate <- as.Date(hud$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hud$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add distance to source
  DistanceToSource1 <- hud$DistanceToSource1
  # Add date and time to hud.pcb.1
  hud.pcb.1 <- cbind(hud.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     season.s, DistanceToSource1)
  # Remove site Bakers Falls. Upstream source
  # North Bakers Falls = WCPCB-HUD006 and
  # South Bakers Falls = WCPCB-HUD006.
  hud.pcb.1 <- subset(hud.pcb.1, SiteID != c("WCPCB-HUD006"))
  hud.pcb.1 <- subset(hud.pcb.1, SiteID != c("WCPCB-HUD010"))
  # Include flow data from USGS station Hudson River
  # sitehudN3 for flow and sitehudN5 for water temperature
  sitehudN3 <- "01328770" # HUDSON RIVER AT THOMSON NY, no temp!
  sitehudN5 <- "01359139" # HUDSON RIVER AT ALBANY NY No flow!
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  # Retrieve USGS data
  # Flow (ft3/s)
  flow.3 <- readNWISdv(sitehudN3, paramflow,
                       min(hud.pcb.1$SampleDate), max(hud.pcb.1$SampleDate))
  # Water temperature in Celsius
  temp <- readNWISdv(sitehudN5, paramtemp,
                     min(hud.pcb.1$SampleDate), max(hud.pcb.1$SampleDate))
  
  # Add USGS data to hud.tpcb.1 matching dates
  hud.pcb.1$flow.3 <- 0.03*flow.3$X_00060_00003[match(hud.pcb.1$SampleDate,
                                                      flow.3$Date)]
  hud.pcb.1$temp <- 273.15 + temp$X_00010_00003[match(hud.pcb.1$SampleDate,
                                                      temp$Date)]
  # Remove samples with flow.3 = NA
  hud.pcb.2 <- hud.pcb.1[!is.na(hud.pcb.1$flow.3), ]
  # Remove metadata not use in the random forest
  hud.pcb.2 <- hud.pcb.2[, !(names(hud.pcb.2) %in% c("SampleDate"))]
}

# Set the seed for reproducibility
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

# Iterate over each numeric column
for (i in seq_along(pcb_numeric_columns)) {
  # Combine numeric and character data
  combined_data <- cbind(hud.pcb.2[, pcb_numeric_columns[i],
                                   drop = FALSE], hud.pcb.2[, char_columns])
  
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
    Location = rep("Hudson River", length(test_data[, 1])),
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
rf_results <- cbind(Location = rep("Hudson River", nrow(rf_results)),
                    rf_results)

# Remove rows in all_results where R_squared < 0
all_results <- all_results %>%
  filter(R_squared >= 0)

# Remove the "R_squared" column from all_results
all_results <- all_results %>% select(-R_squared)

# Export results
write.csv(rf_results,
          file = "Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFPerformancePCB.csv",
          row.names = FALSE)

# Export combined results
write.csv(all_results,
          file = "Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFObsPredPCB.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(all_results, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(0.1, 10^5),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.1, 10^5),
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
ggsave("Output/Plots/Sites/ObsPred/HudsonRiver/HudsonRiverRFPCBV01.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)


