# Code to combine all Passaic River data

# Install packages
install.packages("dplyr")

# Load libraries
library(dplyr)

# Read data ---------------------------------------------------------------
# Data in pg/L
{
  pass <- read.csv("Data/PassaicRiver/pass.csv")
  pass.2 <- read.csv("Data/PassaicRiver/pass02.csv")
  pass.3 <- read.csv("Data/PassaicRiver/pass03.csv")
  pass.4 <- read.csv("Data/PassaicRiver/pass04.csv")
  pass.5 <- read.csv("Data/PassaicRiver/pass05.csv")
  pass.6 <- read.csv("Data/PassaicRiver/pass06.csv")
  pass.7 <- read.csv("Data/PassaicRiver/pass07.csv")
  pass.8 <- read.csv("Data/PassaicRiver/pass08.csv")
  pass.9 <- read.csv("Data/PassaicRiver/pass09.csv")
  pass.10 <- read.csv("Data/PassaicRiver/pass10.csv")
}

# Combine all the dataset in one.
merged_pass <- rbind(pass, pass.2, pass.3, pass.4, pass.5, pass.6, pass.7,
                   pass.8, pass.9, pass.10)

# Delete the first column
merged_pass <- merged_pass[, -1]

# Names and values for the new columns
new_col_names <- c("SampleID", "EPARegion", "StateSampled", "LocationName")
new_col_values <- c("SampleIDValue", "R2", "NJ", "Passaic River")

# Add new columns at the beginning (from column 1)
merged_pass <- cbind(
  setNames(data.frame(matrix(NA, nrow = nrow(merged_pass),
                             ncol = length(new_col_names))), new_col_names),
  merged_pass
)

# Fill the new columns with values
for (i in 1:length(new_col_names)) {
  col_name <- new_col_names[i]
  col_value <- new_col_values[i]
  merged_pass[, col_name] <- col_value
}

# Name and value for the "SiteID" column
site_id_col_name <- "SiteID"
site_id_col_value <- "SiteIDValue"

# Add the "SiteID" column at position 6
merged_pass <- merged_pass %>%
  dplyr::mutate(!!site_id_col_name := site_id_col_value, .before = 6)

# Create a data frame with unique Latitude and Longitude combinations
unique_combinations <- merged_pass %>%
  distinct(Latitude, Longitude)

# Function to generate SiteID based on LATITUDE and LONGITUDE
generate_SiteID <- function(lat, lon) {
  index <- which(unique_combinations$Latitude == lat & unique_combinations$Longitude == lon)
  SiteID <- paste("WCPCB-PASS", sprintf("%03d", index), sep = "")
  return(SiteID)
}

# Add a SiteID column based on LATITUDE and LONGITUDE
merged_pass <- merged_pass %>%
  mutate(SiteID = mapply(generate_SiteID, Latitude, Longitude))

# Define the values for the new columns
PhaseMeasuredValue <- "SurfaceWater"
EPAMethodValue <- "M1668"
AroclorCongenerValue <- "Congener"

# Insert the new columns at position 11
merged_pass <- merged_pass %>%
  mutate(PhaseMeasured = PhaseMeasuredValue, EPAMethod = EPAMethodValue, AroclorCongener = AroclorCongenerValue) %>%
  select(1:10, PhaseMeasured, EPAMethod, AroclorCongener, everything())

# Define the column names
column_names <- c("A1016", "A1221", "A1232", "A1242", "A1248", "A1254", "A1260")

# Initialize these columns with "NA"
merged_pass[, column_names] <- NA

# Insert the columns at position 117
merged_pass <- merged_pass %>%
  select(1:116, all_of(column_names), everything())

# Export results
write.csv(merged_pass, file = "Data/PassaicRiver/PassaicRiverData.csv")


