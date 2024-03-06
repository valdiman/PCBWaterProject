## Water PCB concentrations data analysis
## Fox River 2005 - 2018
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

# Select Fox River data ---------------------------------------------------
fox <- wdc[str_detect(wdc$LocationName, 'Fox River'),]

# Located eastern location & calculate distance to other locations ---------
{
  # Identify the eastern sample based on the maximum Latitude
  index_of_eastern_sample <- which.max(fox$Latitude)
  # Extract coordinates for the eastern sample
  eastern_sample_latitude <- fox$Latitude[index_of_eastern_sample]
  eastern_sample_longitude <- fox$Longitude[index_of_eastern_sample]
  # Define source coordinates for the eastern sample
  eastern_source <- c(Latitude = eastern_sample_latitude,
                      Longitude = eastern_sample_longitude)
  # Create an sf point for the eastern source
  eastern_source_sf <- st_sfc(st_point(c(eastern_source["Longitude"],
                                         eastern_source["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(eastern_source_sf) <- 4326
  # Transform eastern_source_sf to UTM Zone 16N (EPSG:32610)
  eastern_source_sf_utm <- st_transform(eastern_source_sf, 32616)
  # Convert the data frame to an sf object for the eastern sample (fox)
  sf_fox <- st_as_sf(fox, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_fox <- st_set_crs(sf_fox, 4326)
  # Transform to UTM Zone 16N
  sf_fox_utm <- st_transform(sf_fox, 32616)
  # Calculate distances in meters from each location to eastern source
  distances_meters_fox <- st_distance(sf_fox_utm, eastern_source_sf_utm)
  # Convert distances to kilometers
  distances_km_fox <- units::set_units(distances_meters_fox, "km")
  # Extract numeric values and assign to the DistanceToEasternSource column
  fox$DistanceToEasternLocation <- as.numeric(distances_km_fox[, 1])
}

# tPCB data preparation ---------------------------------------------------
{
  # Change date format
  fox$SampleDate <- as.Date(fox$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(fox$SampleDate),
                                  min(as.Date(fox$SampleDate)), units = "days"))
  # Create individual code for each site sampled
  site.numb <- fox$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(fox$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  fox.tpcb <- cbind(factor(fox$SiteID), fox$SampleDate, as.matrix(fox$tPCB),
                    data.frame(time.day), season.s, fox$DistanceToEasternLocation)
  # Add column names
  colnames(fox.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                          "DistanceToEasternLocation")
}

# Remove site -------------------------------------------------------------
# Remove site Lake Winnebago (background site)
fox.tpcb <- subset(fox.tpcb, SiteID != c("WCPCB-FOX001"))

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
                     min(fox.tpcb$date), max(fox.tpcb$date))
  temp <- readNWISdv(sitefoxN2, paramtemp,
                     min(fox.tpcb$date), max(fox.tpcb$date))
  # Add USGS data to fox.tpcb.2, matching dates, conversion to m3/s
  fox.tpcb$flow <- 0.03*flow$X_.Primary.Stream.Flow._00060_00003[match(fox.tpcb$date,
                                                                       flow$Date)]
  fox.tpcb$temp <- 273.15 + temp$X_00010_00003[match(fox.tpcb$date,
                                                     temp$Date)]
  # Remove samples with temp = NA
  fox.tpcb <- na.omit(fox.tpcb)
}

# Random Forest Model tPCB ------------------------------------------------
# Train-Test Split
set.seed(123)

# Train-test split
train_indices <- sample(1:nrow(fox.tpcb), 0.8 * nrow(fox.tpcb))
train_data <- fox.tpcb[train_indices, ]
test_data <- fox.tpcb[-train_indices, ]

# Fit the GBM model
gbm_model <- gbm(log10(tPCB) ~ time + SiteID + season + flow + temp +
                   DistanceToEasternLocation, data = train_data,
                 distribution = "gaussian",  # For regression tasks
                 n.trees = 5000,             # Number of trees
                 interaction.depth = 10,      # Interaction depth
                 shrinkage = 0.001)           # Learning rate

# Print summary of the model
summary(gbm_model)

# Make predictions
predictions <- predict(
  object = gbm_model,
  newdata = test_data,
  n.trees = 5000)

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
          file = "Output/Data/Sites/csv/FoxRiver/FoxRiverRFtPCBV02.csv",
          row.names = FALSE)

# Create a data frame for plotting Observations vs Predictions
plot_data <- data.frame(
  Location = rep("Fox River", nrow(test_data)),
  Observed = log10(test_data$tPCB),
  Predicted = predictions
)

# Export results
write.csv(plot_data,
          file = "Output/Data/Sites/csv/FoxRiver/FoxRiverRFObsPredtPCBV02.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(plot_data, aes(x = 10^(Observed), y = 10^(Predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(100, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(100, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue",
              linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue",
              linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Print the plot
print(plotRF)

# Save plot in folder
ggsave("Output/Plots/Sites/ObsPred/FoxRiver/FoxRiverRFtPCBV02.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)
