## Water PCB concentrations data analysis
## Fox River 2005 - 2018
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

# Select Fox River data ---------------------------------------------------
fox <- wdc[str_detect(wdc$LocationName, 'Fox River'),]
# Lake Winnebago is a background site.
# Data preparation --------------------------------------------------------
{
  # Change date format
  fox$SampleDate <- as.Date(fox$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(fox$SampleDate) - min(as.Date(fox$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(fox$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  fox.tpcb <- cbind(factor(fox$SiteID), fox$SampleDate,
                    fox$Latitude, fox$Longitude, as.matrix(fox$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(fox.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "season")
}

# Remove site -------------------------------------------------------------
# Remove site Lake Winnebago (background site)
fox.tpcb.1 <- subset(fox.tpcb, SiteID != c("WCPCB-FOX001"))

# Include USGS flow and temperature data --------------------------------------------------
{
  # Include flow data from USGS station Fox River
  sitefoxN1 <- "04084445" # flow @ OX RIVER AT APPLETON, WI
  sitefoxN2 <- "040851385" # water temperature @ FOX RIVER AT OIL TANK DEPOT AT GREEN BAY, WI
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow <- readNWISdv(sitefoxN1, paramflow,
                     min(fox.tpcb.1$date), max(fox.tpcb.1$date))
  temp <- readNWISdv(sitefoxN2, paramtemp,
                     min(fox.tpcb.1$date), max(fox.tpcb.1$date))
  # Add USGS data to fox.tpcb.2, matching dates, conversion to m3/s
  fox.tpcb.1$flow <- 0.03*flow$X_.Primary.Stream.Flow._00060_00003[match(fox.tpcb.1$date,
                                                                         flow$Date)]
  fox.tpcb.1$temp <- 273.15 + temp$X_00010_00003[match(fox.tpcb.1$date,
                                                       temp$Date)]
  # Remove samples with temp = NA
  fox.tpcb.1 <- na.omit(fox.tpcb.1)
}

# Random Forest Model -----------------------------------------------------
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(fox.tpcb.1), 0.8 * nrow(fox.tpcb.1))
train_data <- fox.tpcb.1[train_indices, ]
test_data <- fox.tpcb.1[-train_indices, ]

# Fit the Model (1)
rf_model.1 <- randomForest(log10(tPCB) ~ time + SiteID + season + flow
                           + temp, data = train_data)

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
          file = "Output/Data/Sites/csv/FoxRiver/FoxRiverRFPerformancetPCB.csv")

# Feature Importance
importance.1 <- importance(rf_model.1)
# Plot features
barplot(importance.1[, 1], names.arg = rownames(importance.1),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting
plot_data.1 <- data.frame(
  Location = rep("Fox River", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions.1
)

# Export results
write.csv(plot_data.1,
          file = "Output/Data/Sites/csv/FoxRiver/FoxRiverRFObsPredtPCB.csv")

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
ggsave("Output/Plots/Sites/ObsPred/FoxRiver/FoxRiverRFtPCBV01.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Individual PCB Analysis -------------------------------------------------
# Prepare data.frame
{
  fox.pcb <- subset(fox, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  fox.pcb <- subset(fox.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  fox.pcb <- log10(fox.pcb)
  # Replace -inf to NA
  fox.pcb <- do.call(data.frame,
                     lapply(fox.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  fox.pcb.1 <- fox.pcb[,
                       -which(colSums(is.na(fox.pcb))/nrow(fox.pcb) > 0.7)]
  
  # Add site ID
  SiteID <- factor(fox$SiteID)
  # Change date format
  SampleDate <- as.Date(fox$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(fox$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to fox.pcb.1
  fox.pcb.1 <- cbind(fox.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     season.s)
  # Remove site Lake Winnebago (background site)
  fox.pcb.1 <- subset(fox.pcb.1, SiteID != c("WCPCB-FOX001"))
  # Include flow data from USGS station Fox River
  sitefoxN1 <- "04084445" # flow @ OX RIVER AT APPLETON, WI
  sitefoxN2 <- "040851385" # water temperature @ FOX RIVER AT OIL TANK DEPOT AT GREEN BAY, WI
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow <- readNWISdv(sitefoxN1, paramflow,
                     min(fox.pcb.1$SampleDate), max(fox.pcb.1$SampleDate))
  temp <- readNWISdv(sitefoxN2, paramtemp,
                     min(fox.pcb.1$SampleDate), max(fox.pcb.1$SampleDate))
  # Add USGS data to fox.pcb.1, matching dates, conversion to m3/s
  fox.pcb.1$flow <- 0.03*flow$X_.Primary.Stream.Flow._00060_00003[match(fox.pcb.1$SampleDate,
                                                                        flow$Date)]
  fox.pcb.1$temp <- 273.15 + temp$X_00010_00003[match(fox.pcb.1$SampleDate,
                                                      temp$Date)]
  # Remove samples with temperature = NA
  fox.pcb.2 <- fox.pcb.1[!is.na(fox.pcb.1$temp), ]
}

# Set the seed for reproducibility
set.seed(123)

# Find the numeric columns (columns starting with "PCB")
pcb_numeric_columns <- grep("^PCB", colnames(fox.pcb.2), value = TRUE)

# Find the corresponding character columns
char_columns <- setdiff(colnames(fox.pcb.2), pcb_numeric_columns)

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
  combined_data <- cbind(fox.pcb.2[, pcb_numeric_columns[i],
                                   drop = FALSE], fox.pcb.2[, char_columns])
  
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
  
  # Create a data frame for each column's results
  col_results <- data.frame(
    Location = rep("Fox River", length(test_data[, 1])),
    Congener = rep(pcb_numeric_columns[i], length(test_data[, 1])),
    Actual = test_data[, 1],
    Predicted = predictions
  )
  
  # Bind the data frame to the overall results
  all_results <- rbind(all_results, col_results)
}

# Export results
write.csv(rf_results,
          file = "Output/Data/Sites/csv/FoxRiver/FoxRiverRFPerformancePCB.csv",
          row.names = FALSE)

# Export combined results
write.csv(all_results,
          file = "Output/Data/Sites/csv/FoxRiver/FoxRiverRFObsPredPCB.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(all_results, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(0.01, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.01, 10^4),
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
ggsave("Output/Plots/Sites/ObsPred/FoxRiver/FoxRiverRFPCBV01.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)
