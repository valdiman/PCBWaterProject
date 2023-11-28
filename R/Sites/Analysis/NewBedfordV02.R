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

# Data preparation --------------------------------------------------------
{
  # Change date format
  nbh$SampleDate <- as.Date(nbh$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(nbh$SampleDate) - min(as.Date(nbh$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- nbh$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(nbh$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  nbh.tpcb <- cbind(factor(nbh$SiteID), nbh$SampleDate,
                    nbh$Latitude, nbh$Longitude, as.matrix(nbh$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(nbh.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Source locations
source1 <- c(Latitude = 41.675042, Longitude = -70.916295)  # Former Aerovox Site

# Create sf points
source1_sf <- st_point(c(source1["Latitude"], source1["Longitude"]))

# Initialize a list to store point geometries
point_list <- list()

# Loop to create point geometries
for (i in 1:nrow(nbh.tpcb)) {
  point <- st_point(c(nbh.tpcb[i, "Latitude"], nbh.tpcb[i, "Longitude"]))
  point_list[[i]] <- point
}

# Create an sf object from the list of point geometries
nbh.tpcb_sf <- st_sf(geometry = st_sfc(point_list), crs = 4326)

# Ensure that both objects have the same CRS
st_crs(nbh.tpcb_sf) <- st_crs(source1_sf)

# Calculate distances
distances <- st_distance(nbh.tpcb_sf, source1_sf) * 100

# Add distances to hud.tpcb data.frame
nbh.tpcb$Distance_to_source1 <- distances

# Random Forest Model -----------------------------------------------------
# Using all the data, hud.tpcb.1
# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(nbh.tpcb), 0.8 * nrow(nbh.tpcb))
train_data <- nbh.tpcb[train_indices, ]
test_data <- nbh.tpcb[-train_indices, ]

# Fit the Model (1)
rf_model.1 <- randomForest(log10(tPCB) ~ time + site.code + season +
                             Distance_to_source1, data = train_data)

# Make Predictions
predictions.1 <- predict(rf_model.1, newdata = test_data)

# Evaluate Model Performance
mse.1 <- mean((predictions.1 - log10(test_data$tPCB))^2)
rmse.1 <- sqrt(mse.1)
mae.1 <- mean(abs(predictions.1 - log10(test_data$tPCB)))
r_squared.1 <- 1 - (sum((log10(test_data$tPCB) - predictions.1)^2)/sum((log10(test_data$tPCB) - mean(log10(test_data$tPCB)))^2))

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

# Create the scatter plot using ggplot2
ggplot(plot_data.1, aes(x = 10^(Actual), y = 10^(Predicted))) +
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

ggplot(plot_data.1, aes(x = Actual, y = Predicted)) +
  geom_point(shape = 21, size = 3, fill = "white") +
  xlab("Observed concentration PCB (pg/L)") +
  ylab("Predicted concentration PCB (pg/L)") +
  geom_abline(intercept = 0.30103, slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = -0.30103, slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15)

# Create a data frame with observed and predicted values
compare_df <- data.frame(observed = test_data$tPCB,
                         predicted = 10^predictions.1)

# Estimate a factor of 2 between observations and predictions
compare_df$factor2 <- compare_df$observed/compare_df$predicted

# Calculate the percentage of observations within the factor of 2
factor2_percentage <- nrow(compare_df[compare_df$factor2 > 0.5 & compare_df$factor2 < 2, ])/nrow(compare_df) * 100



