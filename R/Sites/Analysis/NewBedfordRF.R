## Water PCB concentrations data analysis per site
## New Bedford Harbor
## Background information
## When the cleanup began, the areas with the highest levels
## of PCBs were addressed first. A 5-acre northern portion
## of the Acushnet River estuary was identified as the
## "hot spot" area (about 14,000 yd3 of sediment exceeding
## a concentration of about 4,000 mg/kg total PCB) and was
## addressed prior to the start of the full scale dredging
## in the upper and lower harbor. This cleanup took place
## from 1994 to 1995 and the off-site disposal of the
## resulting highly contaminated material was completed in 2000.
## More info:
## https://19january2021snapshot.epa.gov/new-bedford-harbor/general-information-about-new-bedford-harbor-cleanup_.html
## https://semspub.epa.gov/work/01/100013466.pdf

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

# Select nbh River data ---------------------------------------------------
nbh <- wdc[str_detect(wdc$LocationName, 'New Bedford'),]

# Add distance to source --------------------------------------------------
{
  # Define source coordinates
  source <- c(Latitude = 41.675042, Longitude = -70.916295)  # Former Aerovox Site
  # Create an sf point for the source
  source_sf <- st_sfc(st_point(c(source["Longitude"], source["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(source_sf) <- 4326
  # Transform source1_sf to UTM Zone 19N (EPSG:32618)
  source_sf_utm <- st_transform(source_sf, 32619)
  # Convert the data frame to an sf object
  sf_nbh <- st_as_sf(nbh, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_nbh <- st_set_crs(sf_nbh, 4326)
  # Transform to UTM Zone 19N
  sf_nbh_utm <- st_transform(sf_nbh, 32619)
  # Calculate distances in meters from each location to source
  distances_meters <- st_distance(sf_nbh_utm, source_sf_utm)
  # Convert distances to kilometers
  distances_km <- units::set_units(distances_meters, "km")
  # Extract numeric values and assign to the DistanceToSource column
  nbh$DistanceToSource <- as.numeric(distances_km[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  nbh$SampleDate <- as.Date(nbh$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(nbh$SampleDate) - min(as.Date(nbh$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(nbh$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  nbh.tpcb <- cbind(factor(nbh$SiteID), nbh$SampleDate, as.matrix(nbh$tPCB),
                    data.frame(time.day), season.s, nbh$DistanceToSource)
  # Add column names
  colnames(nbh.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                          "DistanceToSource")
}

# Random Forest Model -----------------------------------------------------
# Using all the data, nbh.tpcb.1
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(nbh.tpcb), 0.8 * nrow(nbh.tpcb))
train_data <- nbh.tpcb[train_indices, ]
test_data <- nbh.tpcb[-train_indices, ]

# Fit the Model (1)
rf_model.1 <- randomForest(log10(tPCB) ~ time + SiteID + season +
                             DistanceToSource, data = train_data)

# Make Predictions
predictions.1 <- predict(rf_model.1, newdata = test_data)

# Evaluate Model Performance
mse <- mean((predictions.1 - log10(test_data$tPCB))^2)
rmse <- sqrt(mse)
r_squared <- 1 - (sum((log10(test_data$tPCB) - predictions.1)^2)/sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2))

# Estimate a factor of 2 between observations and predictions
# Create a data frame with observed and predicted values
compare_df <- data.frame(observed = test_data$tPCB,
                         predicted = 10^predictions.1)

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
          file = "Output/Data/Sites/csv/NewBedfordHarbor/NBHRFPerformancetPCB.csv",
          row.names = FALSE)

# Feature Importance
importance.1 <- importance(rf_model.1)
# Plot features
barplot(importance.1[, 1], names.arg = rownames(importance.1),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting and exporting
plot_data.1 <- data.frame(
  Location = rep("New Bedford Harbor", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions.1
)

# Export results
write.csv(plot_data.1,
          file = "Output/Data/Sites/csv/NewBedfordHarbor/NBHRFObsPredtPCB.csv",
          row.names = FALSE)

# Create the scatter plot using ggplot2
plotRF <- ggplot(plot_data.1, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10^2, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^2, 10^7),
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
ggsave("Output/Plots/Sites/ObsPred/NewBedfordHarbor/NBHRFtPCBV01.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Individual PCB Analysis -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  nbh.pcb <- subset(nbh, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  nbh.pcb <- subset(nbh.pcb, select = -c(A1016:DistanceToSource))
  # Log10 individual PCBs 
  nbh.pcb <- log10(nbh.pcb)
  # Replace -inf to NA
  nbh.pcb <- do.call(data.frame,
                     lapply(nbh.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  nbh.pcb.1 <- nbh.pcb[,
                       -which(colSums(is.na(nbh.pcb))/nrow(nbh.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(nbh$SiteID)
  # Change date format
  SampleDate <- as.Date(nbh$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(nbh$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add distance to source
  DistanceToSource <- nbh$DistanceToSource
  # Add date and time to nbh.pcb.1
  nbh.pcb.1 <- cbind(nbh.pcb.1, SiteID, SampleDate,
                     data.frame(time.day), season.s, DistanceToSource)
  # Remove metadata not use in the random forest
  nbh.pcb.1 <- nbh.pcb.1[, !(names(nbh.pcb.1) %in% c("SampleDate"))]
}

# Set the seed for reproducibility
set.seed(123)

# Find the numeric columns (columns starting with "PCB")
pcb_numeric_columns <- grep("^PCB", colnames(nbh.pcb.1), value = TRUE)

# Find the corresponding character columns
char_columns <- setdiff(colnames(nbh.pcb.1), pcb_numeric_columns)

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
  combined_data <- cbind(nbh.pcb.1[, pcb_numeric_columns[i],
                                   drop = FALSE], nbh.pcb.1[, char_columns])
  
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
    Location = rep("New Bedford Harbor", length(test_data[, 1])),
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
          file = "Output/Data/Sites/csv/NewBedfordHarbor/NBHRFPerformancePCB.csv",
          row.names = FALSE)

# Export combined results
write.csv(all_results,
          file = "Output/Data/Sites/csv/NewBedfordHarbor/NBHRFObsPredPCB.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(all_results, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(1, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^7),
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
ggsave("Output/Plots/Sites/ObsPred/NewBedfordHarbor/NBHRFPCBV01.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)

