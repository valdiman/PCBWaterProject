## Water PCB concentrations data analysis per site
## Housatonic River
## Aroclors 1254 and 1260, no congener analysis
## GE facility map @ https://semspub.epa.gov/work/01/574882.pdf
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

# Select Housatonic River data ---------------------------------------------------
hou <- wdc[str_detect(wdc$LocationName, 'Housatonic River'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  hou$SampleDate <- as.Date(hou$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(hou$SampleDate) - min(as.Date(hou$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- hou$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hou$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  hou.tpcb <- cbind(factor(hou$SiteID), hou$SampleDate,
                    hou$Latitude, hou$Longitude, as.matrix(hou$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(hou.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Source locations + distance to samples
{
  source1 <- c(Latitude =  42.456479, Longitude = -73.217587)  # GE Pittsfield
  # Create sf points
  source1_sf <- st_point(c(source1["Latitude"], source1["Longitude"])) #Euclidean (straight-line) distance
  # Initialize a list to store point geometries
  point_list <- list()
  # Loop to create point geometries
  for (i in 1:nrow(hou.tpcb)) {
    point <- st_point(c(hou.tpcb[i, "Latitude"], hou.tpcb[i, "Longitude"])) # check this
    point_list[[i]] <- point
  }
  # Create an sf object from the list of point geometries
  hou.tpcb_sf <- st_sf(geometry = st_sfc(point_list), crs = 4326)
  # Ensure that both objects have the same CRS
  st_crs(hou.tpcb_sf) <- st_crs(source1_sf)
  # Calculate distances
  distances1 <- st_distance(hou.tpcb_sf, source1_sf) * 100
  # Add distances to hud.tpcb data.frame
  hou.tpcb$Distance_to_source1 <- distances1
}


# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Housatonic River
{
  siteHouN1 <- "01197000" # EAST BRANCH HOUSATONIC RIVER AT COLTSVILLE, MA
  siteHouN2 <- "01197500" # HOUSATONIC RIVER NEAR GREAT BARRINGTON, MA
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  # Retrieve USGS data
  flow.1 <- readNWISdv(siteHouN1, paramflow,
                       min(hou.tpcb$date), max(hou.tpcb$date))
  flow.2 <- readNWISdv(siteHouN2, paramflow,
                       min(hou.tpcb$date), max(hou.tpcb$date))
  # Add USGS data to hou.tpcb, matching dates (m3/s, 0.03 conversion factor)
  hou.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(hou.tpcb$date,
                                                     flow.1$Date)]
  hou.tpcb$flow.2 <- 0.03*flow.2$X_00060_00003[match(hou.tpcb$date,
                                                     flow.2$Date)]
}

# Random Forest Model -----------------------------------------------------
# Using all the data, hou.tpcb
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(hou.tpcb), 0.8 * nrow(hou.tpcb))
train_data <- hou.tpcb[train_indices, ]
test_data <- hou.tpcb[-train_indices, ]

# Fit the Model (1)
rf_model.1 <- randomForest(log10(tPCB) ~ time + site.code + season + flow.1 + flow.2 +
                             Distance_to_source1, data = train_data)

# Make Predictions
predictions.1 <- predict(rf_model.1, newdata = test_data)

# Evaluate Model Performance
mse.1 <- mean((predictions.1 - log10(test_data$tPCB))^2)
rmse.1 <- sqrt(mse.1) # Report
r_squared.1 <- 1 - (sum((log10(test_data$tPCB) - predictions.1)^2)/sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2)) # Report

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
          file = "Output/Data/Sites/csv/HousatonicRiver/HousatonicRiverRFPerformancetPCB.csv")

# Feature Importance
importance.1 <- importance(rf_model.1)
# Plot features
barplot(importance.1[, 1], names.arg = rownames(importance.1),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting
plot_data.1 <- data.frame(
  Actual = log10(test_data$tPCB),
  Predicted = predictions.1
)

# Export results
write.csv(plot_data.1,
          file = "Output/Data/Sites/csv/HousatonicRiver/HousatonicRiverRFtPCB.csv")

# Create the scatter plot using ggplot2
plotRF <- ggplot(plot_data.1, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10^3, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^3, 10^7),
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
ggsave("Output/Plots/Sites/ObsPred/HousatonicRiver/HousatonicRiverRFV01.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Fit the Model (2)
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(hou.tpcb), 0.8 * nrow(hou.tpcb))
train_data <- hou.tpcb[train_indices, ]
test_data <- hou.tpcb[-train_indices, ]

# Fit the Model
rf_model.2 <- randomForest(log10(tPCB) ~ time + site.code + season + flow.1 +
                             Distance_to_source1, data = train_data)

# Make Predictions
predictions.2 <- predict(rf_model.2, newdata = test_data)

# Evaluate Model Performance
mse.2 <- mean((predictions.2 - log10(test_data$tPCB))^2)
rmse.2 <- sqrt(mse.2)
mae.2 <- mean(abs(predictions.2 - log10(test_data$tPCB)))
r_squared.2 <- 1 - (sum((log10(test_data$tPCB) - predictions.2)^2)/sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2))

# Feature Importance
importance.2 <- importance(rf_model.2)
# Plot features
barplot(importance.2[, 1], names.arg = rownames(importance.2),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting
plot_data.2 <- data.frame(
  Actual = log10(test_data$tPCB),
  Predicted = predictions.2
)

# Create the scatter plot using ggplot2
ggplot(plot_data.2, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10^3, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^3, 10^7),
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

# Estimate a factor of 2 between observations and predictions
# Create a data frame with observed and predicted values
compare_df.2 <- data.frame(observed = test_data$tPCB,
                           predicted = 10^predictions.2)

# Estimate a factor of 2 between observations and predictions
compare_df.2$factor2 <- compare_df.2$observed/compare_df.2$predicted

# Calculate the percentage of observations within the factor of 2
factor2_percentage.2 <- nrow(compare_df.2[compare_df.2$factor2 > 0.5 & compare_df.2$factor2 < 2, ])/nrow(compare_df.2)*100

# Remove site -------------------------------------------------------------
## Remove sites upstream source
hou.tpcb.1 <- subset(hou.tpcb, SiteID != c("WCPCB-HOU004")) # Crane Paper Company
hud.tpcb.1 <- subset(hou.tpcb.1, SiteID != c("WCPCB-HOU014")) # Hubbard Ave. Bridge

# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(hou.tpcb.1), 0.8 * nrow(hou.tpcb.1))
train_data <- hou.tpcb.1[train_indices, ]
test_data <- hou.tpcb.1[-train_indices, ]

# Fit the Model (3)
rf_model.3 <- randomForest(log10(tPCB) ~ time + site.code + season + flow.1 +
                             Distance_to_source1, data = train_data)

# Make Predictions
predictions.3 <- predict(rf_model.3, newdata = test_data)

# Evaluate Model Performance
mse.3 <- mean((predictions.3 - log10(test_data$tPCB))^2)
rmse.3 <- sqrt(mse.3)
mae.3 <- mean(abs(predictions.3 - log10(test_data$tPCB)))
r_squared.3 <- 1 - (sum((log10(test_data$tPCB) - predictions.3)^2)/sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2))

# Feature Importance
importance.3 <- importance(rf_model.3)
# Plot features
barplot(importance.3[, 1], names.arg = rownames(importance.3),
        main = "Feature Importance", las = 2, cex.names = 0.7)

# Create a data frame for plotting
plot_data.3 <- data.frame(
  Actual = log10(test_data$tPCB),
  Predicted = predictions.3
)

# Create the scatter plot using ggplot2
ggplot(plot_data.3, aes(x = 10^(Actual), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10^3, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^3, 10^7),
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

# Estimate a factor of 2 between observations and predictions
# Create a data frame with observed and predicted values
compare_df.3 <- data.frame(observed = test_data$tPCB,
                           predicted = 10^predictions.3)

# Estimate a factor of 2 between observations and predictions
compare_df.3$factor2 <- compare_df.3$observed/compare_df.3$predicted

# Calculate the percentage of observations within the factor of 2
factor2_percentage.3 <- nrow(compare_df.3[compare_df.3$factor2 > 0.5 & compare_df.3$factor2 < 2, ])/nrow(compare_df.3)*100

# To predict future data --------------------------------------------------
# Assuming new_data is your new dataset without the tPCB variable
predictions <- predict(rf_model.3, newdata = new_data)


