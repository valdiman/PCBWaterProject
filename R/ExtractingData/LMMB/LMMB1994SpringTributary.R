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
LMMB_data <- read_csv("Data/LMMB/1994SpringTributary.csv")

# Select only samples with PCB measurements
LMMB_data <- LMMB_data[LMMB_data$ANL_CODE_21 == "PCB-tot", ]

# Remove columns with other chemicals
LMMB_data <- LMMB_data[, -c(1:9, 13:15, 20:140)]

# Function to clean and convert VALUE_ columns to numeric
clean_and_convert_value <- function(x) {
  # Remove asterisks and any non-numeric characters except '.'
  cleaned_value <- gsub("[^0-9.*-]", "", x)
  # Replace asterisk with an empty string
  cleaned_value <- gsub("\\*", "", cleaned_value)
  # Handle empty cells by converting them to NA
  cleaned_value[cleaned_value == ""] <- NA
  # Convert to numeric
  as.numeric(cleaned_value)
}

# Identify the VALUE_ columns
value_columns <- grep("^VALUE_", names(LMMB_data), value = TRUE)

# Clean and convert each VALUE_ column
LMMB_data[value_columns] <- lapply(LMMB_data[value_columns],
                                   clean_and_convert_value)

# Define the pivot columns
pivot_columns <- c(paste0("ANL_CODE_", 21:568),
                   paste0("ANALYTE_", 21:568),
                   paste0("VALUE_", 21:568),
                   paste0("UNITS_", 21:568),
                   paste0("FRACTION_", 21:568),
                   paste0("RESULT_REMARK_", 21:568))

# Pivot the data
reshaped_data <- LMMB_data %>%
  pivot_longer(cols = all_of(pivot_columns), 
               names_to = c(".value", "set"),
               names_sep = "_"
  )

# Select samples with Fraction as Filtrate
filtered_data <- reshaped_data %>%
  filter(FRACTION == "Filtrate")

