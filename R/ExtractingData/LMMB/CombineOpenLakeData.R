# Code to combine all opebutary River data from LMMB

# Install packages
install.packages("dplyr")

# Load libraries
{
  library(dplyr)
}

# Read data ---------------------------------------------------------------
# Data in pg/L
{
  ope.1 <- read.csv("Data/LMMB/OpenLake/1994SOL.csv")
  ope.2 <- read.csv("Data/LMMB/OpenLake/1994SuOL.csv")
  ope.3 <- read.csv("Data/LMMB/OpenLake/1995SOL.csv")
  ope.4 <- read.csv("Data/LMMB/OpenLake/1995SuOL.csv")
}

# Combine all the dataset in one.
merged_ope <- rbind(ope.1, ope.2, ope.3, ope.4)

# Delete the first column
merged_ope <- merged_ope[, -1]

# Rename SAMPLE_ID to SiteName and change the name of the SAMPLE_ID
names(merged_ope)[names(merged_ope) == "SAMPLE_ID"] <- "SiteName"
merged_ope$SiteName <- sub("^[^-]*-", "", merged_ope$SiteName)

# Names and values for the new columns
new_col_names <- c("SampleID", "EPARegion", "StateSampled", "LocationName", "SiteID")
new_col_values <- c("SampleIDValue", "R5", "NA", "LMMB", "SiteIDValue")

# Add new columns at the beginning (from column 1)
merged_ope <- cbind(setNames(data.frame(matrix(NA, nrow = nrow(merged_ope),
                                               ncol = length(new_col_names))),
                             new_col_names), merged_ope)

# Fill the new columns with values
for (i in 1:length(new_col_names)) {
  col_name <- new_col_names[i]
  col_value <- new_col_values[i]
  merged_ope[, col_name] <- col_value
}

# Organize column order
desired_column_order <- c(
  "SampleID",
  "EPARegion",
  "StateSampled",
  "LocationName",
  "SiteName",
  "SiteID",
  "SampleDate",
  "Latitude",
  "Longitude",
  "Units"
)

merged_ope <- merged_ope %>%
  select(desired_column_order, everything())

## Until here !
# Add State to StateSampled column
merged_ope$StateSampled[merged_ope$SiteName %in% c(1, "1-7m", "1-8m", 3,
                                                   "3-10m", "3-8m")] <- "IN"
merged_ope$StateSampled[merged_ope$SiteName %in% c(5, "MB9")] <- "IL"
merged_ope$StateSampled[merged_ope$SiteName %in% c(17, "MB21", "MB21-5m", "MB21-24m",
                                                   280, 240,
                                                   "MB25", 31, 110, 140,
                                                   180, "40M", "MB38",
                                                   45, "GB100M", "GB24M",
                                                   "GB17")] <- "WI"





# Create a data frame with unique Latitude and Longitude combinations
unique_combinations <- merged_ope %>%
  distinct(Latitude, Longitude)

# Function to generate SiteID based on LATITUDE and LONGITUDE
generate_SiteID <- function(lat, lon) {
  index <- which(unique_combinations$Latitude == lat & unique_combinations$Longitude == lon)
  SiteID <- paste("WCPCB-LMM", sprintf("%03d", index), sep = "")
  return(SiteID)
}

# Add a SiteID column based on LATITUDE and LONGITUDE
merged_ope <- merged_ope %>%
  mutate(SiteID = mapply(generate_SiteID, Latitude, Longitude))

# Create a new column "SampleDate_yyyymmdd" with the date in yyyymmdd format
merged_ope$SampleDate_yyyymmdd <- format(as.Date(merged_ope$SampleDate,
                                                  format = "%m/%d/%y"), format = "%Y%m%d")

# Create SampleID by concatenating SiteID and SampleDate_yyyymmdd
merged_ope$SampleID <- paste(merged_ope$SiteID,
                              merged_ope$SampleDate_yyyymmdd, sep = "-")

# Create a new column "SampleCount" to count samples with the same SampleID
merged_ope <- merged_ope %>%
  group_by(SampleID) %>%
  mutate(SampleCount = row_number())

# Add a dot and SampleCount to SampleID if SampleCount is greater than 1
for (i in unique(merged_ope$SampleID)) {
  subset_df <- merged_ope[merged_ope$SampleID == i, ]
  num_samples <- nrow(subset_df)
  if (num_samples > 1) {
    for (j in 1:num_samples) {
      subset_df$SampleID[j] <- paste(subset_df$SampleID[j], ".", j, sep = "")
    }
    merged_ope[merged_ope$SampleID == i, ] <- subset_df
  } else {
    merged_ope[merged_ope$SampleID == i, "SampleID"] <- paste(i, ".1", sep = "")
  }
}

# Remove the SampleCount and SampleDate_yyyymmdd columns if no longer needed
merged_ope$SampleCount <- NULL
merged_ope$SampleDate_yyyymmdd <- NULL

# Define the values for the new columns
PhaseMeasuredValue <- "SurfaceWater"
EPAMethodValue <- "M1668"
AroclorCongenerValue <- "Congener"

# Insert the new columns at position 11
merged_ope <- merged_ope %>%
  mutate(PhaseMeasured = PhaseMeasuredValue, EPAMethod = EPAMethodValue,
         AroclorCongener = AroclorCongenerValue) %>%
  select(1:10, PhaseMeasured, EPAMethod, AroclorCongener, everything())

# Define the column names
column_names <- c("A1016", "A1221", "A1232", "A1242", "A1248", "A1254", "A1260")

# Initialize these columns with "NA"
merged_ope[, column_names] <- NA

# Insert the columns at position 118
merged_ope <- merged_ope %>%
  select(1:117, all_of(column_names), everything())

# Export results
write.csv(merged_pass, file = "Data/PassaicRiver/PassaicRiverData.csv")


