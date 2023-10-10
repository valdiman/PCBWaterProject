 # https://cdxapps.epa.gov/cdx-glenda/action/querytool/querySystem

# Install packages
install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")
install.packages("jsonlite")
install.packages("lubridate")
install.packages("stringr")

# Load libraries
{
  library(dplyr)
  library(tidyr)
  library(readr)
  library(jsonlite)
  library(lubridate)
  library(stringr)
}


# Read the CSV file into a data frame
LMMB_data <- read_csv("Data/LMMB/1994SpringTributary.csv")

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

# Identify the range of columns for chemical data
start_column <- which(names(LMMB_data) == "ANL_CODE_1")
pivot_columns <- names(LMMB_data)[start_column:length(LMMB_data)]

# Pivot the data while keeping metadata columns and maintaining the "ANAL_CODE" column
LMMB_data_long <- LMMB_data %>%
  pivot_longer(
    cols = all_of(pivot_columns),
    names_to = c(".value", "Chemical"),
    names_pattern = "^(ANL_CODE|ANALYTE|VALUE|UNITS|FRACTION|RESULT_REMARK)_([0-9]+)$"
  ) %>%
  filter(!is.na(VALUE))

# Select samples
LMMB_data_long <- LMMB_data_long %>%
  filter(FRACTION == "Filtrate", UNITS == "ng/l",
         grepl("^PCB", ANL_CODE),
         !grepl("^field blank", QC_TYPE),
         !grepl("^Composite QC", SAMPLE_TYPE))

# Remove the "Chemical" column
LMMB_data_long <- LMMB_data_long %>%
  select(-Row, -PROJECT, -PROJ_CODE, -YEAR, -MONTH, -SEASON,
         -CRUISE_ID, -VISIT_ID, -STATION_ID, -TIME_ZONE,
         -STN_DEPTH_M, -SAMPLE_DEPTH_M, -SAMPLE_TYPE, -QC_TYPE,
         -Chemical, -FRACTION, -RESULT_REMARK)

# Change date format
LMMB_data_long$SAMPLING_DATE <- sub(" \\d+:\\d+", "",
                                    LMMB_data_long$SAMPLING_DATE)

# Step 1: Extract into separate columns
LMMB_data_long <- LMMB_data_long %>%
  separate(ANL_CODE, into = c("Prefix", "FirstNum", "SecondNum", "ThirdNum"), sep = "_|\\+")

# Step 2: Remove leading zeros
LMMB_data_long$FirstNum <- sub("^0+", "", LMMB_data_long$FirstNum)
LMMB_data_long$SecondNum <- sub("^0+", "", LMMB_data_long$SecondNum)
LMMB_data_long$ThirdNum <- sub("^0+", "", LMMB_data_long$ThirdNum)

# Step 3: Replace NAs with empty strings for all numeric columns
numeric_cols <- c("FirstNum", "SecondNum", "ThirdNum")
LMMB_data_long[numeric_cols] <- lapply(LMMB_data_long[numeric_cols], function(x) ifelse(is.na(x), "", x))

# Step 4: Concatenate with a dot, considering all numbers
LMMB_data_long$Modified_ANL_CODE <- paste0(LMMB_data_long$Prefix,
                                           ifelse(LMMB_data_long$FirstNum == "", "", LMMB_data_long$FirstNum),
                                           ifelse(LMMB_data_long$SecondNum == "", "", paste0(".", LMMB_data_long$SecondNum)),
                                           ifelse(LMMB_data_long$ThirdNum == "", "", paste0(".", LMMB_data_long$ThirdNum)))

# Replace "PCB-tot" with "tPCB" in Modified_ANL_CODE column
LMMB_data_long$Modified_ANL_CODE <- gsub("PCB-tot", "tPCB", LMMB_data_long$Modified_ANL_CODE)

# Remove temporary columns if needed
LMMB_data_long <- LMMB_data_long %>%
  select(-Prefix, -FirstNum, -SecondNum, -ThirdNum,
         -ANALYTE)




# Read the JSON file with new congener list from code NewPCBList.R
pcb_groups <- read_json("Data/pcb_groups.json")



# Reshape the dataset
transposed_data <- LMMB_data_long %>%
  pivot_wider(
    id_cols = c(LATITUDE, LONGITUDE, SAMPLING_DATE, SAMPLE_ID, UNITS),
    names_from = ANALYTE,
    values_from = AVG_VALUE
  )
