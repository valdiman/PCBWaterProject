## Water PCB concentrations analysis.
## Data were obtained from EPA and contractors from PCB Superfund
## sites in USA. Using log10 of the sum of PCB.
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

# Total Concentration Analysis --------------------------------------------
{
  # Change date format
  SampleDate <- as.Date(wdc$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Create individual code for each site sampled
  site.numb <- wdc$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(wdc$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  tpcb <- cbind(factor(wdc$SiteID), SampleDate, wdc$tPCB,
                data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(tpcb) <- c("SiteID", "date", "tPCB", "time", "site.code",
                      "season")
}

# Random Forest Model tPCB ------------------------------------------------
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(tpcb), 0.8 * nrow(tpcb))
train_data <- tpcb[train_indices, ]
test_data <- tpcb[-train_indices, ]

# Fit the Model (1)
rf_model <- randomForest(log10(tPCB) ~ time + site.code + season,
                           data = train_data)

# Make Predictions
predictions <- predict(rf_model, newdata = test_data)

# Evaluate Model Performance
mse <- mean((predictions - log10(test_data$tPCB))^2)
rmse <- sqrt(mse)
r_squared <- 1 - (sum((log10(test_data$tPCB) - predictions)^2)/sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2))

# Estimate a factor of 2 between observations and predictions for the normal values
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
          file = "Output/Data/Global/csv/RFtPCB.csv")

# Feature Importance
importance <- importance(rf_model)
# Plot features
barplot(importance[, 1], names.arg = rownames(importance),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting
plot_data.1 <- data.frame(
  Location = rep("USA", nrow(test_data)),
  Actual = log10(test_data$tPCB),
  Predicted = predictions
)

# Export results
write.csv(plot_data.1,
          file = "Output/Data/Global/csv/RFObsPredtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(plot_data.1, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(0.1, 10^8),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.1, 10^8),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue",
              linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Print the plot
print(plotRF)

# Save plot in folder
ggsave("Output/Plots/Global/RFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Random Forest Model individual PCBs -------------------------------------
{
  # Only consider congener data
  wdc.pcb <- subset(wdc, AroclorCongener == "Congener")
  # Create individual code for each site sampled
  site.numb <- wdc.pcb$SiteID %>% as.factor() %>% as.numeric
  # Change date format
  SampleDate <- as.Date(wdc.pcb$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(wdc.pcb$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Remove metadata
  wdc.pcb <- subset(wdc.pcb, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor and tPCB data
  wdc.pcb <- subset(wdc.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  wdc.pcb <- log10(wdc.pcb)
  # Replace -inf to NA
  wdc.pcb <- do.call(data.frame,
                     lapply(wdc.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  wdc.pcb.1 <- wdc.pcb[,
                       -which(colSums(is.na(wdc.pcb))/nrow(wdc.pcb) > 0.7)]
  # Add date and time to wdc.pcb.1
  wdc.pcb.1 <- cbind(wdc.pcb.1, data.frame(time.day), site.numb,
                     season.s)
}

# Set the seed for reproducibility
set.seed(123)

# Find the numeric columns (columns starting with "PCB")
pcb_numeric_columns <- grep("^PCB", colnames(wdc.pcb.1), value = TRUE)

# Find the corresponding character columns
char_columns <- setdiff(colnames(wdc.pcb.1), pcb_numeric_columns)

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
  combined_data <- cbind(wdc.pcb.1[, pcb_numeric_columns[i],
                                   drop = FALSE], wdc.pcb.1[, char_columns])
  
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
  R_squared <- 1 - (sum((test_data[, 1] - predictions)^2) / sum((test_data[, 1] - mean(test_data[, 1]))^2))
  
  # Calculate factor2_percentage within the loop with normal scale
  compare_df <- data.frame(
    observed = 10^(test_data[, 1]),
    predicted = 10^(predictions)
  )
  compare_df$factor2 <- compare_df$observed / compare_df$predicted
  factor2_percentage <- sum(compare_df$factor2 > 0.5 & compare_df$factor2 < 2) / nrow(compare_df) * 100
  
  # Store the results in the matrix
  rf_results[i, 2:4] <- c(mse, R_squared, factor2_percentage)
  
  # Create a data frame for each column's results
  col_results <- data.frame(
    Location = rep("USA", length(test_data[, 1])),
    Congener = rep(pcb_numeric_columns[i], length(test_data[, 1])),
    Actual = test_data[, 1],
    Predicted = predictions,
    R_squared = R_squared  # Add R_squared
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
          file = "Output/Data/Global/csv/RFPCB.csv",
          row.names = FALSE)

# Export combined results
write.csv(all_results,
          file = "Output/Data/Global/csv/RFObsPredPCB.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(all_results, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(0.001, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.001, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
  ylab(expression(bold("Predicted concentration PCBi (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# To check values inside the 2:1 and 2:1 lines
# Build the plot
plot_data <- ggplot_build(plotRFPCBi)$data[[1]]
# Count points inside the 2 lines
points_inside_lines <- sum(plot_data$y >= plot_data$x - log10(2) & plot_data$y <= plot_data$x + log10(2))
# Print the count
print(points_inside_lines)

# Count points outside the 2 lines
points_outside_lines <- sum(plot_data$y < plot_data$x - log10(2) | plot_data$y > plot_data$x + log10(2))
# Print the count
print(points_outside_lines)

# Print the plot
print(plotRFPCBi)

# Save plot in folder
ggsave("Output/Plots/Global/RFPCB.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)

