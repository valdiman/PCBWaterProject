## Water PCB concentrations data analysis per site
## Chesapeake Bay & Delaware Canal

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
  time.day <- data.frame(as.Date(che$SampleDate) - min(as.Date(che$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- che$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(che$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  che.tpcb <- cbind(factor(che$SiteID), che$SampleDate, as.matrix(che$tPCB),
                    data.frame(time.day), site.numb, season.s,
                    che$DistanceToCentroid)
  # Add column names
  colnames(che.tpcb) <- c("SiteID", "date", "tPCB", "time", "site.code",
                          "season", "DistanceToCentroid")
}

# Add water temperature data ----------------------------------------------
{
  # Read water temperature
  wtp <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayWT.csv")
  # Convert date columns to Date format
  wtp$Date <- as.Date(wtp$Date)
  # Add water temperature to grl.tpcb
  che.tpcb$temp <- wtp$WTMP_K[match(che.tpcb$date, wtp$Date)]
  # Remove samples with temp = NA
  che.tpcb <- na.omit(che.tpcb)
}

# Random Forest Model -----------------------------------------------------
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(che.tpcb), 0.8 * nrow(che.tpcb))
train_data <- che.tpcb[train_indices, ]
test_data <- che.tpcb[-train_indices, ]

# Fit the Model (1)
rf_model.1 <- randomForest(log10(tPCB) ~ time + site.code + season + temp +
                             DistanceToCentroid,
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
          file = "Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFPerformancetPCB.csv")

# Feature Importance
importance.1 <- randomForest::importance(rf_model.1)
# Plot features
barplot(importance.1[, 1], names.arg = rownames(importance.1),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting
plot_data.1 <- data.frame(
  Location = rep("Chesapeake Bay", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions.1
)

# Export results
write.csv(plot_data.1,
          file = "Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFObsPredtPCB.csv")

# Create the scatter plot
plotRF <- ggplot(plot_data.1, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10, 10^5),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10, 10^5),
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
ggsave("Output/Plots/Sites/ObsPred/ChesapeakeBay/ChesapeakeBayRFtPCBV01.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Individual PCB Analysis -------------------------------------------------
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
  
  # Create individual code for each site sampled
  site.numb <- che$SiteID %>% as.factor() %>% as.numeric
  # Change date format
  SampleDate <- as.Date(che$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Add distance to the centroid
  centroid <- che$DistanceToCentroid
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(che$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to grl.pcb.1
  che.pcb.1 <- cbind(che.pcb.1, SampleDate, data.frame(time.day),
                     site.numb, season.s, centroid)
  # Add water temperature to grl.pcb.1
  che.pcb.1$temp <- wtp$WTMP_K[match(che.pcb.1$SampleDate,
                                               wtp$Date)]
  # Remove samples with temperature = NA
  che.pcb.2 <- che.pcb.1[!is.na(che.pcb.1$temp), ]
}

# Set the seed for reproducibility
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

# Iterate over each numeric column
for (i in seq_along(pcb_numeric_columns)) {
  # Combine numeric and character data
  combined_data <- cbind(che.pcb.2[, pcb_numeric_columns[i],
                                   drop = FALSE], che.pcb.2[, char_columns])
  
  # Exclude rows with missing values
  combined_data <- na.omit(combined_data)
  
  # Sample indices for training
  train_indices <- sample(1:nrow(combined_data), 0.8 * nrow(combined_data))
  
  # Create separate training and testing sets
  train_data <- combined_data[train_indices, ]
  test_data <- combined_data[-train_indices, ]
  
  # Remove the SampleDate column from the training data
  train_data <- train_data[, -which(names(train_data) == "SampleDate")]
  
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
    Location = rep("Chesapeake Bay", length(test_data[, 1])),
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
          file = "Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFPerformancePCBV01.csv",
          row.names = FALSE)

# Export combined results
write.csv(all_results,
          file = "Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFObsPredPCBV01.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(all_results, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(1, 10^5),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^5),
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
ggsave("Output/Plots/Sites/ObsPred/ChesapeakeBay/ChesapeakeBayRFPCBV01.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)
