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

# Generate a vector of column names from "VALUE_1" to "VALUE_127"
value_columns <- paste0("VALUE_", 1:127)

# Columns to fix, including logic columns
columns_to_fix <- c("VALUE_123", "VALUE_124", "VALUE_125", "VALUE_126", "VALUE_127")

# Define a function to handle the transformation
transform_column <- function(x) {
  # Replace "*INVALID" and numbers with an asterisk (*) with NA
  x <- ifelse(grepl("^\\*[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", x), NA_real_, x)
  # Replace values starting with "*" with NA
  x <- ifelse(grepl("^\\*", x), NA_real_, x)
  # Remove asterisks and convert to numeric (including scientific notation)
  as.numeric(gsub("[*]", "", x, perl = TRUE))
}

# Apply transformations using lapply for each column
LMMB_data[, value_columns] <- lapply(LMMB_data[, value_columns], transform_column)

# Apply transformations for logic columns
LMMB_data[, columns_to_fix] <- lapply(LMMB_data[, columns_to_fix], function(x) as.integer(x != "INVALID"))

# Remove asterisks from all cells in the data frame
LMMB_data <- data.frame(lapply(LMMB_data, function(x) gsub("\\*", "", x)))

# Work on reshaping the data
# Remove the "STAT_TYPE" column
LMMB_data <- LMMB_data %>%
  select(-STAT_TYPE)

# Identify the range of columns for chemical data (columns 20 to 781)
pivot_columns <- names(LMMB_data)[20:781]

# Pivot the data while keeping metadata columns and maintaining the "ANAL_CODE" column
LMMB_data_long <- LMMB_data %>%
  pivot_longer(
    cols = all_of(pivot_columns),
    names_to = c(".value", "Chemical"),
    names_sep = "_"
  ) %>%
  filter(!is.na(VALUE))  # Remove rows with NA values in VALUE column

# Remove unnecessary columns
LMMB_data_long <- LMMB_data_long %>%
  select(-Row, -PROJECT, -PROJ_CODE, -YEAR, -MONTH, -SEASON, -CRUISE_ID,
         -TIME_ZONE, -STN_DEPTH_M, -SAMPLE_DEPTH_M, -VISIT_ID, -Chemical, -ANL, -RESULT)

# Check for data problems
data_problems <- problems(LMMB_data_long)

# Display data problems, if any
if (nrow(data_problems) > 0) {
  print(data_problems)
} else {
  cat("No data problems found.\n")
}

# Data look good!
# Select correct samples
LMMB_data_long <- LMMB_data_long %>%
  filter(FRACTION == "Filtrate", SAMPLE_TYPE == "Individual", UNITS == "ng/l")

# Change date format
LMMB_data_long$SAMPLING_DATE <- sub(" \\d+:\\d+", "", LMMB_data_long$SAMPLING_DATE)




