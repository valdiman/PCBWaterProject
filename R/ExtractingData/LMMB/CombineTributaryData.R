# Code to combine all Tributary River data from LMMB

# Install packages
install.packages("dplyr")

# Load libraries
{
  library(dplyr)
}

# Read data ---------------------------------------------------------------
# Data in pg/L
{
  tri.1 <- read.csv("Data/LMMB/Tributaries/1994T.csv")
  tri.2 <- read.csv("Data/LMMB/Tributaries/1995T.csv")
}

# Combine all the dataset in one.
merged_tri <- rbind(tri.1, tri.2)

# Delete the first column
merged_tri <- merged_tri[, -1]

# Define a mapping between codes and names
code_to_name <- c(
  "TFOXRB01" = "Tributary Fox River", #WI
  "TGRANH01" = "Tributary Grand River", # MI
  "TIHCAE01" = "Tributary Indiana Harbor Canal", # IN
  "TKALAG01" = "Tributary Kalamazoo River", # MI
  "TMANIK06" = "Tributary Manistique River", #WI?
  "TMENOA01" = "Tributary Menominee River", #WI
  "TMILWD04" = "Tributary Milwaukee River", #WI
  "TMUSKI01" = "Tributary Muskegon River", #MI
  "TPEREJ05" = "Tributary Pere Marquette River", #MI
  "TSHEBC01" = "Tributary Sheboygan River", #WI
  "TSTJOF05" = "Tributary St. Joseph River" #MI
)


# Names and values for the new columns
new_col_names <- c("SampleID", "EPARegion", "StateSampled", "LocationName",
                   "SiteName")
new_col_values <- c("SampleIDValue", "R5", "NA", "LMMB", "NA")

# Add new columns at the beginning (from column 1)
merged_tri <- cbind(setNames(data.frame(matrix(NA, nrow = nrow(merged_tri),
                                               ncol = length(new_col_names))),
                             new_col_names), merged_tri)

# Fill the new columns with values
for (i in 1:length(new_col_names)) {
  col_name <- new_col_names[i]
  col_value <- new_col_values[i]
  merged_tri[, col_name] <- col_value
}

# Name and value for the "SiteID" column
site_id_col_name <- "SiteID"
site_id_col_value <- "SiteIDValue"

# Add the "SiteID" column at position 6
merged_tri <- merged_tri %>%
  dplyr::mutate(!!site_id_col_name := site_id_col_value, .before = 6)

# Create a data frame with unique Latitude and Longitude combinations
unique_combinations <- merged_tri %>%
  distinct(Latitude, Longitude)

# Function to generate SiteID based on LATITUDE and LONGITUDE
generate_SiteID <- function(lat, lon) {
  index <- which(unique_combinations$Latitude == lat & unique_combinations$Longitude == lon)
  SiteID <- paste("WCPCB-LMM", sprintf("%03d", index), sep = "")
  return(SiteID)
}

# Add a SiteID column based on LATITUDE and LONGITUDE
merged_tri <- merged_tri %>%
  mutate(SiteID = mapply(generate_SiteID, Latitude, Longitude))

# Create a new column "SampleDate_yyyymmdd" with the date in yyyymmdd format
merged_tri$SampleDate_yyyymmdd <- format(as.Date(merged_tri$SampleDate,
                                                  format = "%m/%d/%y"), format = "%Y%m%d")

# Create SampleID by concatenating SiteID and SampleDate_yyyymmdd
merged_tri$SampleID <- paste(merged_tri$SiteID,
                              merged_tri$SampleDate_yyyymmdd, sep = "-")

# Create a new column "SampleCount" to count samples with the same SampleID
merged_tri <- merged_tri %>%
  group_by(SampleID) %>%
  mutate(SampleCount = row_number())

# Add a dot and SampleCount to SampleID if SampleCount is greater than 1
for (i in unique(merged_tri$SampleID)) {
  subset_df <- merged_tri[merged_tri$SampleID == i, ]
  num_samples <- nrow(subset_df)
  if (num_samples > 1) {
    for (j in 1:num_samples) {
      subset_df$SampleID[j] <- paste(subset_df$SampleID[j], ".", j, sep = "")
    }
    merged_tri[merged_tri$SampleID == i, ] <- subset_df
  } else {
    merged_tri[merged_tri$SampleID == i, "SampleID"] <- paste(i, ".1", sep = "")
  }
}

# Remove the SampleCount and SampleDate_yyyymmdd columns if no longer needed
merged_tri$SampleCount <- NULL
merged_tri$SampleDate_yyyymmdd <- NULL

# Define the values for the new columns
PhaseMeasuredValue <- "SurfaceWater"
EPAMethodValue <- "M1668"
AroclorCongenerValue <- "Congener"

# Insert the new columns at position 11
merged_tri <- merged_tri %>%
  mutate(PhaseMeasured = PhaseMeasuredValue, EPAMethod = EPAMethodValue,
         AroclorCongener = AroclorCongenerValue) %>%
  select(1:10, PhaseMeasured, EPAMethod, AroclorCongener, everything())

# Define the column names
column_names <- c("A1016", "A1221", "A1232", "A1242", "A1248", "A1254", "A1260")

# Initialize these columns with "NA"
merged_tri[, column_names] <- NA

# Insert the columns at position 118
merged_tri <- merged_tri %>%
  select(1:117, all_of(column_names), everything())

# Export results
write.csv(merged_pass, file = "Data/PassaicRiver/PassaicRiverData.csv")


