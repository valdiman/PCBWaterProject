 # https://cdxapps.epa.gov/cdx-glenda/action/querytool/querySystem

# Install packages
install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")
install.packages("jsonlite")
install.packages("lubridate")

# Load libraries
{
  library(dplyr)
  library(tidyr)
  library(readr)
  library(jsonlite)
  library(lubridate)
}

# Read the CSV file into a data frame
LMMB_data <- read_csv("Data/LMMB/results.csv")

# Remove the STAT_TYPE column
LMMB_data <- LMMB_data %>%
  select(-Row, -STAT_TYPE)

# Filter rows based on the SAMPLE_TYPE condition, removing QC samples
LMMB_data <- LMMB_data %>%
  filter(SAMPLE_TYPE == "Individual")

# Convert the datetime column to a DateTime object
LMMB_data$SAMPLING_DATE <- mdy_hm(LMMB_data$SAMPLING_DATE)

# Format the datetime column as "mm/dd/yyyy"
LMMB_data$SAMPLING_DATE <- format(LMMB_data$SAMPLING_DATE, format = "%m/%d/%Y")

# Select the metadata columns (from 'PROJECT' to 'SAMPLE_ID')
metadata <- LMMB_data %>%
  select(PROJECT:SAMPLE_ID)



# Combine all chemical columns into a single dataset
chemical_data <- LMMB_data %>%
  select(ANL_CODE_1:`...782`)  # Adjust column names as needed


# Assuming you have 'metadata' and 'chemical_data' data frames as defined earlier

# Combine metadata and chemical data by adding a unique ID to both data frames
metadata$ID <- seq_len(nrow(metadata))
chemical_data$ID <- seq_len(nrow(chemical_data))

# Split the chemical data into groups of 6 columns each (excluding ID)
num_chemicals <- 6
chemical_data_list <- lapply(seq(2, ncol(chemical_data), by = num_chemicals), function(start_col) {
  end_col <- min(start_col + num_chemicals - 1, ncol(chemical_data))
  select(chemical_data, starts_with("ANL_CODE_"))[, start_col:end_col, drop = FALSE]
})



# Combine the split chemical data into a single data frame
chemical_data_combined <- do.call(cbind, chemical_data_list)

# Merge metadata and combined chemical data using the 'ID' column
final_data <- merge(metadata, chemical_data_combined, by = "ID")

# Remove the 'ID' column
final_data <- final_data[, -which(names(final_data) == "ID")]

# Now, 'final_data' contains the metadata and chemical data with 6 columns each.
