 # https://cdxapps.epa.gov/cdx-glenda/action/querytool/querySystem

# Install packages
install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")
install.packages("jsonlite")
install.packages("lubridate")
install.packages("strings")

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

# Remove the metadata column not necessaries
LMMB_data_long <- LMMB_data_long %>%
  select(-Row, -PROJECT, -PROJ_CODE, -YEAR, -MONTH, -SEASON,
         -CRUISE_ID, -VISIT_ID, -STATION_ID, -TIME_ZONE,
         -STN_DEPTH_M, -SAMPLE_DEPTH_M, -SAMPLE_TYPE, -QC_TYPE,
         -Chemical, -FRACTION, -RESULT_REMARK)

# Remove PCB congeners with 2 measurements (both measurements are the same)
LMMB_data_long <- LMMB_data_long %>%
  group_by(LATITUDE, LONGITUDE, SAMPLING_DATE, MEDIUM, SAMPLE_ID, ANL_CODE, ANALYTE) %>%
  mutate(VALUE = ifelse(row_number() == 1, VALUE, NA)) %>%
  ungroup() %>%
  na.omit()

# Change date format
LMMB_data_long$SAMPLING_DATE <- sub(" \\d+:\\d+", "",
                                    LMMB_data_long$SAMPLING_DATE)

# Function to change PCB names
transform_ANL_CODE <- function(anl_code) {
  # Replace "PCB-" with "PCB"
  anl_code <- gsub("^PCB-", "PCB", anl_code)
  
  # Replace "PCB" with "PCB"
  anl_code <- gsub("^PCB", "PCB", anl_code)
  
  # Extract numbers using regular expressions
  numbers <- str_extract_all(anl_code, "\\d+")
  
  if (length(numbers) > 0) {
    # Sort numbers in ascending order
    sorted_numbers <- sort(as.integer(unlist(numbers)))
    
    # Replace dots with plus signs and combine sorted numbers
    transformed_code <- paste("PCB", gsub("\\.", "+", paste(sorted_numbers, collapse = ".")), sep = "")
  } else {
    transformed_code <- anl_code
  }
  
  return(transformed_code)
}

# Apply function to change the PCB names
LMMB_data_long$Transformed_ANL_CODE <- sapply(LMMB_data_long$ANL_CODE,
                                              transform_ANL_CODE)

# Reshape the dataset
transposed_data <- LMMB_data_long %>%
  pivot_wider(
    id_cols = c(LATITUDE, LONGITUDE, SAMPLING_DATE, SAMPLE_ID, UNITS),
    names_from = Transformed_ANL_CODE,
    values_from = VALUE
  )

# Multiply PCB columns by 1000
pcb_columns <- names(transposed_data)[grep("^PCB\\d+", names(transposed_data))]
transposed_data[, pcb_columns] <- transposed_data[, pcb_columns] * 1000

# Multiply all values in the "PCB" column by 1000
transposed_data$PCB <- transposed_data$PCB * 1000

# Update RESULT_UNIT to "pg/l"
transposed_data$UNITS <- "pg/l"

# Rename the "PCB" column to "tPCB"
colnames(transposed_data)[colnames(transposed_data) == "PCB"] <- "tPCB"

# Reorder the columns so that "tPCB" is at the end
transposed_data <- transposed_data %>%
  select(-tPCB, everything())

# Export results
write.csv(transposed_data, file = "Data/LMMB/1994ST.csv")


## Not sure about this in this code
# Define a mapping between codes and names
code_to_name <- c(
  "TFOXRB01" = "Tributary Fox River",
  "TFOXRB02" = "Tributary Fox River",
  "TFOXRB03" = "Tributary Fox River",
  "TFOXRB04" = "Tributary Fox River",
  "TFOXRB05" = "Tributary Fox River",
  "TGRANH01" = "Tributary Grand River",
  "TGRANH02" = "Tributary Grand River",
  "TGRANH03" = "Tributary Grand River",
  "TGRANH04" = "Tributary Grand River",
  "TGRANH05" = "Tributary Grand River",
  "TGRANH06" = "Tributary Grand River",
  "TGRANH07" = "Tributary Grand River",
  "TGRANH08" = "Tributary Grand River",
  "TKALAG01" = "Tributary Kalamazoo River",
  "TKALAG02" = "Tributary Kalamazoo River",
  "TKALAG03" = "Tributary Kalamazoo River",
  "TKALAG05" = "Tributary Kalamazoo River",
  "TKALAG06" = "Tributary Kalamazoo River",
  "TMANIK01" = "Tributary Manistique River",
  "TMANIK02" = "Tributary Manistique River",
  "TMANIK03" = "Tributary Manistique River",
  "TMANIK04" = "Tributary Manistique River",
  "TMANIK05" = "Tributary Manistique River",
  "TMANIK06" = "Tributary Manistique River",
  "TMENOA01" = "Tributary Menominee River",
  "TMENOA02" = "Tributary Menominee River",
  "TMENOA03" = "Tributary Menominee River",
  "TMENOA04" = "Tributary Menominee River",
  "TMILWD01" = "Tributary Milwaukee River",
  "TMILWD02" = "Tributary Milwaukee River",
  "TMILWD03" = "Tributary Milwaukee River",
  "TMILWD04" = "Tributary Milwaukee River",
  "TMUSKI01" = "Tributary Muskegon River",
  "TMUSKI02" = "Tributary Muskegon River",
  "TMUSKI03" = "Tributary Muskegon River",
  "TMUSKI04" = "Tributary Muskegon River",
  "TMUSKI05" = "Tributary Muskegon River",
  "TPEREJ01" = "Tributary Pere Marquette River",
  "TPEREJ02" = "Tributary Pere Marquette River",
  "TPEREJ03" = "Tributary Pere Marquette River",
  "TPEREJ04" = "Tributary Pere Marquette River",
  "TPEREJ05" = "Tributary Pere Marquette River",
  "TSHEBC01" = "Tributary Sheboygan River",
  "TSHEBC02" = "Tributary Sheboygan River",
  "TSHEBC03" = "Tributary Sheboygan River",
  "TSHEBC04" = "Tributary Sheboygan River",
  "TSHEBC06" = "Tributary Sheboygan River",
  "TSTJOF01" = "Tributary St. Joseph River",
  "TSTJOF02" = "Tributary St. Joseph River",
  "TSTJOF03" = "Tributary St. Joseph River",
  "TSTJOF04" = "Tributary St. Joseph River",
  "TSTJOF05" = "Tributary St. Joseph River"
)

# Create SiteName column
transposed_data <- transposed_data %>%
  mutate(SiteName = code_to_name[SAMPLE_ID])

# Create SiteID Column
transposed_data <- transposed_data %>%
  mutate(SiteID = sprintf("WCPCB-LMM%03d",
                          as.integer(factor(SiteName, levels = unique(SiteName)))))
         



