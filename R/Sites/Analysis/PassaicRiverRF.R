## Water PCB concentrations data analysis per site
## Passaic River
## Source: Diamond Alkali. 80-120 Lister Avenue, Newark, NJ.

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

# Select Passaic River data ---------------------------------------------------
pass <- wdc[str_detect(wdc$LocationName, 'Passaic River'),]

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
  sf_pass <- st_as_sf(pass, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_pass <- st_set_crs(sf_pass, 4326)
  # Transform to UTM Zone 18N
  sf_pass_utm <- st_transform(sf_pass, 32618)
  # Calculate distances in meters from each location to source
  distances_meters <- st_distance(sf_pass_utm, source_sf_utm)
  # Convert distances to kilometers
  distances_km <- units::set_units(distances_meters, "km")
  # Extract numeric values and assign to the DistanceToSource column
  pass$DistanceToSource <- as.numeric(distances_km[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  pass$SampleDate <- as.Date(pass$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(pass$SampleDate) - min(as.Date(pass$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- pass$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(pass$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  pass.tpcb <- cbind(factor(pass$SiteID), pass$SampleDate, as.matrix(pass$tPCB),
                    data.frame(time.day), site.numb, season.s, pass$DistanceToSource)
  # Add column names
  colnames(pass.tpcb) <- c("SiteID", "date", "tPCB", "time", "site.code",
                           "season", "DistanceToSource")
}

# Include USGS flow and temperature data --------------------------------------------------
{
  # Include flow data from USGS station Passaic River
  sitepassN1 <- "01381900" # No temp
  sitepassN2 <- "01379500" # No temp
  sitepassN3 <- "01389005" # No flow
  sitepassN4 <- "01389010" # No temp
  sitepassN5 <- "01389500" # No temp
  sitepassN6 <- "01389890" # No temp

  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow.1 <- readNWISdv(sitepassN1, paramflow,
                     min(pass.tpcb$date), max(pass.tpcb$date))
  flow.2 <- readNWISdv(sitepassN2, paramflow,
                     min(pass.tpcb$date), max(pass.tpcb$date))
  flow.3 <- readNWISdv(sitepassN4, paramflow,
                     min(pass.tpcb$date), max(pass.tpcb$date))
  flow.4 <- readNWISdv(sitepassN5, paramflow,
                     min(pass.tpcb$date), max(pass.tpcb$date))
  flow.5 <- readNWISdv(sitepassN6, paramflow,
                     min(pass.tpcb$date), max(pass.tpcb$date))
  temp <- readNWISdv(sitepassN3, paramtemp,
                     min(pass.tpcb$date), max(pass.tpcb$date))
  
  # Add USGS data to pass.tpcb.2, matching dates, conversion to m3/s
  pass.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(pass.tpcb$date,
                                                      flow.1$Date)]
  pass.tpcb$flow.2 <- 0.03*flow.1$X_00060_00003[match(pass.tpcb$date,
                                                      flow.2$Date)]
  pass.tpcb$flow.3 <- 0.03*flow.1$X_00060_00003[match(pass.tpcb$date,
                                                      flow.3$Date)]
  pass.tpcb$flow.4 <- 0.03*flow.1$X_00060_00003[match(pass.tpcb$date,
                                                      flow.4$Date)]
  pass.tpcb$flow.5 <- 0.03*flow.1$X_00060_00003[match(pass.tpcb$date,
                                                      flow.5$Date)]
  pass.tpcb$temp <- 273.15 + temp$X_.from.middle.intake_00010_00003[match(pass.tpcb$date,
                                                       temp$Date)]
}

# Remove site -------------------------------------------------------------
# Remove site located in the ocean. Possible typo in original coordinates.
pass.tpcb.1 <- subset(pass.tpcb, SiteID != c("WCPCB-PASS022"))

# Random Forest Model -----------------------------------------------------
# Using flow.1
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(pass.tpcb.1), 0.8 * nrow(pass.tpcb.1))
train_data <- pass.tpcb.1[train_indices, ]
test_data <- pass.tpcb.1[-train_indices, ]

# Fit the Model
rf_model.1 <- randomForest(log10(tPCB) ~ time + site.code + season +
                             flow.1 + temp + DistanceToSource,
                           data = train_data)

# Make Predictions
predictions.1 <- predict(rf_model.1, newdata = test_data)

# Evaluate Model Performance
mse <- mean((predictions.1 - log10(test_data$tPCB))^2)
rmse <- sqrt(mse)
r_squared <- 1 - (sum((log10(test_data$tPCB) - predictions.1)^2)/sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2))

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
                             Value = c(rmse, r_squared,
                                       factor2_percentage.1))

# Remove unnecessary columns
performance_df <- performance_df[, !(names(performance_df) %in% c("V1", "V2", "V3"))]

# Print the modified data frame
print(performance_df)

# Export results
write.csv(performance_df,
          file = "Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFPerformancetPCB.csv",
          row.names = FALSE)

# Feature Importance
importance.1 <- importance(rf_model.1)
barplot(importance.1[, 1], names.arg = rownames(importance.1),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting and exporting
plot_data.1 <- data.frame(
  Location = rep("Passaic River", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions.1
)

# Export results
write.csv(plot_data.1,
          file = "Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFObsPredtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(plot_data.1, aes(x = 10^(Actual), y = 10^(Predicted))) +
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
ggsave("Output/Plots/Sites/ObsPred/PassaicRiver/PassaicRiverRFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Individual PCB Analysis -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  pass.pcb <- subset(pass, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  pass.pcb <- subset(pass.pcb, select = -c(A1016:DistanceToSource))
  # Log10 individual PCBs 
  pass.pcb <- log10(pass.pcb)
  # Replace -inf to NA
  pass.pcb <- do.call(data.frame,
                     lapply(pass.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  pass.pcb.1 <- pass.pcb[,
                       -which(colSums(is.na(pass.pcb))/nrow(pass.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(pass$SiteID)
  # Create individual code for each site sampled
  site.numb <- pass$SiteID %>% as.factor() %>% as.numeric
  # Change date format
  SampleDate <- as.Date(pass$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(pass$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add distance to source
  DistanceToSource <- pass$DistanceToSource
  # Add date and time to pass.pcb.1
  pass.pcb.1 <- cbind(pass.pcb.1, SiteID, site.numb, SampleDate,
                      data.frame(time.day), season.s, DistanceToSource)
  # Remove site located in the ocean. Possible typo in original coordinates.
  pass.pcb.1 <- subset(pass.pcb.1, SiteID != c("WCPCB-PASS022"))
  # Include flow data from USGS station Passaic River
  sitepassN1 <- "01381900" # No temp
  sitepassN3 <- "01389005" # No flow
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow.1 <- readNWISdv(sitepassN1, paramflow,
                       min(pass.pcb.1$SampleDate), max(pass.pcb.1$SampleDate))
  temp <- readNWISdv(sitepassN3, paramtemp,
                     min(pass.pcb.1$SampleDate), max(pass.pcb.1$SampleDate))
  # Add USGS data to pass.tpcb.1, matching dates, conversion to m3/s
  pass.pcb.1$flow.1 <- 0.03*flow.1$X_00060_00003[match(pass.pcb.1$SampleDate,
                                                      flow.1$Date)]
  pass.pcb.1$temp <- 273.15 + temp$X_.from.middle.intake_00010_00003[match(pass.pcb.1$SampleDate,
                                                                          temp$Date)]
  # Remove metadata not use in the random forest
  pass.pcb.1 <- pass.pcb.1[, !(names(pass.pcb.1) %in% c("SiteID", "SampleDate"))]
}

# Set the seed for reproducibility
set.seed(123)

# Find the numeric columns (columns starting with "PCB")
pcb_numeric_columns <- grep("^PCB", colnames(pass.pcb.1), value = TRUE)

# Find the corresponding character columns
char_columns <- setdiff(colnames(pass.pcb.1), pcb_numeric_columns)

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
  combined_data <- cbind(pass.pcb.1[, pcb_numeric_columns[i],
                                   drop = FALSE], pass.pcb.1[, char_columns])
  
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
    Location = rep("Passaic River", length(test_data[, 1])),
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
rf_results <- cbind(Location = rep("Passaic River", nrow(rf_results)),
                    rf_results)

# Remove rows in all_results where R_squared < 0
all_results <- all_results %>%
  filter(R_squared >= 0)

# Remove the "R_squared" column from all_results
all_results <- all_results %>% select(-R_squared)

# Export results
write.csv(rf_results,
          file = "Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFPerformancePCB.csv",
          row.names = FALSE)

# Export combined results
write.csv(all_results,
          file = "Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFObsPredPCB.csv",
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
ggsave("Output/Plots/Sites/ObsPred/PassaicRiver/PassaicRiverRFPCB.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)


