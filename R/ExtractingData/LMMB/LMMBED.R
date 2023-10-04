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

# Check issues with the data before any modifications
data_problems <- problems(LMMB_data)

# Display data problems, if any
if (nrow(data_problems) > 0) {
  print(data_problems)
} else {
  cat("No initial data problems found.\n")
}

# Fix issues
# Remove asterisks from all cells in the data frame
LMMB_data <- data.frame(lapply(LMMB_data, function(x) gsub("\\*", "", x)))

# Make values to numeric values
# Select columns that start with "VALUE"
value_columns <- grep("^VALUE", names(LMMB_data), value = TRUE)

# Check if warnings exist and clear them
if (length(warnings()) > 0) {
  warnings()
  warnings(NULL)
}

# Check issues with the data before any modifications
data_problems <- problems(LMMB_data)

# Transform non-numeric columns to numeric values
LMMB_data <- LMMB_data %>%
  mutate_at(vars(all_of(value_columns)), ~{
    numeric_value <- as.numeric(as.character(.))
    ifelse(is.na(numeric_value) | is.infinite(numeric_value), NA, numeric_value)
  })

# Check for problems in the dataframe
data_problems <- problems(LMMB_data)

# Print the rows and columns with problems
print(data_problems)

# No problems, but still some columns that need to be numeric are not numeric
# Check the format of each selected column
format_info <- sapply(LMMB_data[, value_columns], class)
print(format_info)

# Select columns that need fixing
columns_to_fix <- c("VALUE_123", "VALUE_124", "VALUE_125", "VALUE_126", "VALUE_127")

# Loop through the selected columns and convert problematic values
for (col in columns_to_fix) {
  LMMB_data[[col]] <- as.numeric(as.character(LMMB_data[[col]]))
  # Identify non-numeric values and replace them with NA
  LMMB_data[[col]][!is.na(LMMB_data[[col]]) & !is.finite(LMMB_data[[col]])] <- NA
}

# Verify that the columns are now numeric
str(LMMB_data[, columns_to_fix])

# Work on the reshaping the data
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

# Remove a few columns not needed.
LMMB_data_long <- LMMB_data_long %>%
  select(-Row, -PROJECT, -PROJ_CODE, -YEAR, -MONTH, -SEASON, -CRUISE_ID,
         -TIME_ZONE, -STN_DEPTH_M, -SAMPLE_DEPTH_M, -VISIT_ID, -Chemical, -ANL, -RESULT)

LMMB_data_long <- LMMB_data_long %>%
  filter(FRACTION == "Filtrate", SAMPLE_TYPE == "Individual", UNITS == "ng/l")

# Change date format
LMMB_data_long$SAMPLING_DATE <- sub(" \\d+:\\d+", "", LMMB_data_long$SAMPLING_DATE)



