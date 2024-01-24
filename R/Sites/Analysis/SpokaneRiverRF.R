## Water PCB concentrations data analysis per site
## Spokane River
## Random forest model

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

# Select Spokane River data ---------------------------------------------------
spo <- wdc[str_detect(wdc$LocationName, 'Spokane River'),]

# Located eastern location & calculate distance to other locations ---------
{
  # Identify the eastern sample based on the maximum Latitude
  index_of_eastern_sample <- which.max(spo$Latitude)
  # Extract coordinates for the eastern sample
  eastern_sample_latitude <- spo$Latitude[index_of_eastern_sample]
  eastern_sample_longitude <- spo$Longitude[index_of_eastern_sample]
  # Define source coordinates for the eastern sample
  eastern_source <- c(Latitude = eastern_sample_latitude,
                      Longitude = eastern_sample_longitude)
  # Create an sf point for the eastern source
  eastern_source_sf <- st_sfc(st_point(c(eastern_source["Longitude"],
                                         eastern_source["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(eastern_source_sf) <- 4326
  # Transform eastern_source_sf to UTM Zone 10N (EPSG:32610)
  eastern_source_sf_utm <- st_transform(eastern_source_sf, 32610)
  # Convert the data frame to an sf object for the eastern sample (spo)
  sf_spo <- st_as_sf(spo, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_spo <- st_set_crs(sf_spo, 4326)
  # Transform to UTM Zone 10N
  sf_spo_utm <- st_transform(sf_spo, 32610)
  # Calculate distances in meters from each location to eastern source
  distances_meters_spo <- st_distance(sf_spo_utm, eastern_source_sf_utm)
  # Convert distances to kilometers
  distances_km_spo <- units::set_units(distances_meters_spo, "km")
  # Extract numeric values and assign to the DistanceToEasternSource column
  spo$DistanceToEasternLocation <- as.numeric(distances_km_spo[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  spo$SampleDate <- as.Date(spo$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(spo$SampleDate) - min(as.Date(spo$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(spo$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  spo.tpcb <- cbind(factor(spo$SiteID), spo$SampleDate,
                    as.matrix(spo$tPCB), data.frame(time.day), season.s,
                    spo$DistanceToEasternLocation)
  # Add column names
  colnames(spo.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                          "DistanceToEasternLocation")
  # Include flow data from USGS station Spokane River
  siteSpoN1 <- "12417650" # SPOKANE RIVER BLW BLACKWELL NR COEUR D ALENE ID
  siteSpoN2 <- "12419000" # Spokane River near Post Falls, ID
  siteSpoN3 <- "12422500" # Spokane River at Spokane, WA
  siteSpoN4 <- "12424000" # Hangman Creek at Spokane, WA
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  #paramtemp <- "00010" # water temperature, C No data
  # Retrieve USGS data
  flow.1 <- readNWISdv(siteSpoN1, paramflow,
                       min(spo.tpcb$date), max(spo.tpcb$date))
  flow.2 <- readNWISdv(siteSpoN2, paramflow,
                       min(spo.tpcb$date), max(spo.tpcb$date))
  flow.3 <- readNWISdv(siteSpoN3, paramflow,
                       min(spo.tpcb$date), max(spo.tpcb$date))
  flow.4 <- readNWISdv(siteSpoN4, paramflow,
                       min(spo.tpcb$date), max(spo.tpcb$date))
  
  # Add USGS data to spo.tpcb, matching dates
  spo.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(spo.tpcb$date,
                                                     flow.1$Date)]
  spo.tpcb$flow.2 <- 0.03*flow.2$X_00060_00003[match(spo.tpcb$date,
                                                     flow.2$Date)]
  spo.tpcb$flow.3 <- 0.03*flow.3$X_00060_00003[match(spo.tpcb$date,
                                                     flow.3$Date)]
  spo.tpcb$flow.4 <- 0.03*flow.4$X_00060_00003[match(spo.tpcb$date,
                                                     flow.4$Date)]
}

# Remove site -------------------------------------------------------------
## Sample sites not located at the Spokane River
{
  spo.tpcb.1 <- subset(spo.tpcb, SiteID != c("WCPCB-SPR002")) # City of Spokane WRF
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR005")) # Regional WRF
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR006")) # Inland Empire paper
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR008")) # Kaiser Aluminum
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR010")) # Liberty Lake sewer
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR011")) # Post Falls WWTP
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR013")) # Coeur d'Alene WWTP
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR015")) # Hangman Creek
}

# Random Forest Model -----------------------------------------------------
# Using all the data, spo.tpcb
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(spo.tpcb.1), 0.8 * nrow(spo.tpcb.1))
train_data <- spo.tpcb.1[train_indices, ]
test_data <- spo.tpcb.1[-train_indices, ]

# Fit the Model. Flow.3
rf_model.1 <- randomForest(log10(tPCB) ~ time + SiteID + season + flow.3 +
                             DistanceToEasternLocation, data = train_data)

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
          file = "Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFPerformancetPCBV02.csv")

# Feature Importance
importance.1 <- randomForest::importance(rf_model.1)

# Plot features
barplot(importance.1[, 1], names.arg = rownames(importance.1),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting and exporting
plot_data.1 <- data.frame(
  Location = rep("Spokane River", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions.1
)

# Export results
write.csv(plot_data.1,
          file = "Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFObsPredtPCBV02.csv",
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
ggsave("Output/Plots/Sites/ObsPred/SpokaneRiver/SpokaneRiverRFtPCBV02.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Individual PCB Analysis -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  spo.pcb <- subset(spo, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  spo.pcb <- subset(spo.pcb, select = -c(A1016:DistanceToEasternLocation))
  # Log10 individual PCBs 
  spo.pcb <- log10(spo.pcb)
  # Replace -inf to NA
  spo.pcb <- do.call(data.frame,
                     lapply(spo.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  spo.pcb.1 <- spo.pcb[, colSums(is.na(spo.pcb))/nrow(spo.pcb) <= 0.7]
  # Add site ID
  SiteID <- factor(spo$SiteID)
  # Change date format
  SampleDate <- as.Date(spo$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(spo$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add distance to eastern location
  eastern <- spo$DistanceToEasternLocation
  # Add date and time to spo.pcb.1
  spo.pcb.1 <- cbind(spo.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     season.s, eastern)
  # Include flow data from USGS station Spokane River
  siteSpoN3 <- "12422500" # Spokane River at Spokane, WA
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  # Retrieve USGS data
  flow.3 <- readNWISdv(siteSpoN3, paramflow,
                       min(spo.pcb.1$SampleDate), max(spo.pcb.1$SampleDate))
  # Add USGS data to spo.tpcb, matching dates
  spo.pcb.1$flow.3 <- 0.03*flow.3$X_00060_00003[match(spo.pcb.1$SampleDate,
                                                     flow.3$Date)]
  # Sample sites not located at the Spokane River
  spo.pcb.2 <- subset(spo.pcb.1, SiteID != c("WCPCB-SPR002")) # City of Spokane WRF
  spo.pcb.2 <- subset(spo.pcb.2, SiteID != c("WCPCB-SPR005")) # Regional WRF
  spo.pcb.2 <- subset(spo.pcb.2, SiteID != c("WCPCB-SPR006")) # Inland Empire paper
  spo.pcb.2 <- subset(spo.pcb.2, SiteID != c("WCPCB-SPR008")) # Kaiser Aluminum
  spo.pcb.2 <- subset(spo.pcb.2, SiteID != c("WCPCB-SPR010")) # Liberty Lake sewer
  spo.pcb.2 <- subset(spo.pcb.2, SiteID != c("WCPCB-SPR011")) # Post Falls WWTP
  spo.pcb.2 <- subset(spo.pcb.2, SiteID != c("WCPCB-SPR013")) # Coeur d'Alene WWTP
  spo.pcb.2 <- subset(spo.pcb.2, SiteID != c("WCPCB-SPR015")) # Hangman Creek
  
  # Remove metadata not use in the random forest
  spo.pcb.2 <- spo.pcb.2[, !(names(spo.pcb.2) %in% c("SampleDate"))]
}

# Set the seed for reproducibility
set.seed(123)

# Find the numeric columns (columns starting with "PCB")
pcb_numeric_columns <- grep("^PCB", colnames(spo.pcb.2), value = TRUE)

# Find the corresponding character columns
char_columns <- setdiff(colnames(spo.pcb.2), pcb_numeric_columns)

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
  combined_data <- cbind(spo.pcb.2[, pcb_numeric_columns[i],
                                   drop = FALSE], spo.pcb.2[, char_columns])
  
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
  
  # Create a data frame for each column's results, including R_squared
  col_results <- data.frame(
    Location = rep("Spokane River", length(test_data[, 1])),
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

# Remove rows in all_results where R_squared < 0
all_results <- all_results %>%
  filter(R_squared >= 0)

# Remove the "R_squared" column from all_results
all_results <- all_results %>% select(-R_squared)

# Export results
write.csv(rf_results,
          file = "Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFPerformancePCBV02.csv",
          row.names = FALSE)

# Export combined results
write.csv(all_results,
          file = "Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFObsPredPCBV02.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(all_results, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(0.01, 10^3),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.01, 10^3),
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
ggsave("Output/Plots/Sites/ObsPred/SpokaneRiver/SpokaneRiverRFPCBV02.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)

