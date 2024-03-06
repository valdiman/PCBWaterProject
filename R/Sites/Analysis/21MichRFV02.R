## Water PCB concentrations data analysis per site
## Data from DEQ Michigan
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

# Select 21 Michigan data ---------------------------------------------------
mic <- wdc[str_detect(wdc$LocationName, '21Mich'),]

# Calculate central location. ---------------------------------------------
{
  # Calculate the mean latitude and longitude
  center_lat <- mean(mic$Latitude)
  center_lon <- mean(mic$Longitude)
  # Create a data frame for the center
  center_df <- data.frame(
    SiteID = "Center",
    Latitude = center_lat,
    Longitude = center_lon
  )
  # Convert the center data frame to an sf object
  sf_center <- st_as_sf(center_df, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326) for the center
  sf_center <- st_set_crs(sf_center, 4326)
  # Transform to UTM Zone 16 for the center
  sf_center_utm <- st_transform(sf_center, 32616)
  # Convert the data frame to an sf object
  sf_mic <- st_as_sf(mic, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_mic <- st_set_crs(sf_mic, 4326)
  # Transform to UTM Zone 16
  sf_mic <- st_transform(sf_mic, 32616)
  # Calculate distances in meters from each location to the center
  distances_meters <- st_distance(sf_mic, sf_center_utm)
  # Convert distances to kilometers
  distances_km <- units::set_units(distances_meters, "km")
  # Extract numeric values and assign to the DistanceToCentroid column
  mic$DistanceToCentroid <- as.numeric(distances_km[, 1])
}

# tPCB data preparation ---------------------------------------------------
{
  # Change date format
  mic$SampleDate <- as.Date(mic$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(mic$SampleDate),
                                  min(as.Date(mic$SampleDate)), units = "days"))
  # Create individual code for each site sampled
  site.numb <- mic$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(mic$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  mic.tpcb <- cbind(factor(mic$SiteID), mic$SampleDate, as.matrix(mic$tPCB),
                    data.frame(time.day), site.numb, season.s,
                    mic$DistanceToCentroid)
  
  # Add column names
  colnames(mic.tpcb) <- c("SiteID", "date", "tPCB", "time", "site.code",
                          "season", "DistanceToCentroid")
}

# Random Forest Model tPCB ------------------------------------------------
# Set seed for reproducibility
set.seed(123)

# Train-test split
train_indices <- sample(1:nrow(mic.tpcb), 0.8 * nrow(mic.tpcb))
train_data <- mic.tpcb[train_indices, ]
test_data <- mic.tpcb[-train_indices, ]

# Fit the GBM model
gbm_model <- gbm(log10(tPCB) ~ time + SiteID + season + DistanceToCentroid,
                 data = train_data,
                 distribution = "gaussian",  # For regression tasks
                 n.trees = 4000,             # Number of trees
                 interaction.depth = 5,      # Interaction depth
                 shrinkage = 0.001)           # Learning rate

# Print summary of the model
summary(gbm_model)

# Make predictions
predictions <- predict(
  object = gbm_model,
  newdata = test_data,
  n.trees = 4000)

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
          file = "Output/Data/Sites/csv/21Mich/21MichRFtPCBV02.csv",
          row.names = FALSE)

# Create a data frame for plotting Observations vs Predictions
plot_data <- data.frame(
  Location = rep("21 Mich", nrow(test_data)),
  Observed = log10(test_data$tPCB),
  Predicted = predictions
)

# Export results
write.csv(plot_data,
          file = "Output/Data/Sites/csv/21Mich/21MichRFObsPredtPCBV02.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(plot_data, aes(x = 10^(Observed), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(1, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue",
              linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = -log10(2), slope = 1, col = "blue",
              linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Print the plot
print(plotRF)

# Save plot in folder
ggsave("Output/Plots/Sites/ObsPred/21Mich/21MichRFObsPredtPCBV02.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

