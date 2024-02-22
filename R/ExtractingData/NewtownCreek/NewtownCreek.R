
# Install packages
install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")

# Load libraries
{
  library(dplyr)
  library(tidyr)
  library(readr)
  library(jsonlite)
  library(stringr)
}

# Read the CSV file into a data frame
NC_data <- read_csv("Data/NewtownCreek/Bi-B10-8_Analytical-Results.csv")

# Arrange the data to ensure it's ordered correctly (if needed)
NC_data <- NC_data %>% arrange(sys_sample_code, sample_date)

# Filter and keep only the columns you need
NC_data <- NC_data %>%
  select(sys_sample_code, sample_date, matrix_code, fraction, analytic_method,
         y_coord, x_coord, target_unit, chemical_name, result_value, dilution_factor)

# Filter rows based on the analytic_method and fraction condition
NC_data <- NC_data %>%
  filter(analytic_method == "E1668A",
         fraction == "D")

# Modify the CHEMICAL_NAME column to extract "PCB X" where X is the number
NC_data <- NC_data %>%
  mutate(chemical_name = sub(".*\\(PCB (\\d-)\\).*", "PCB \\1", chemical_name))

transform_chemical_name <- function(chem_code) {
  # Replace "PCB-" with "PCB"
  chem_code <- gsub("PCB-", "PCB", chem_code)
  
  # Extract numbers using regular expressions
  numbers <- str_extract_all(chem_code, "\\d+")
  
  if (length(numbers) > 0) {
    # Sort numbers in ascending order
    sorted_numbers <- sort(as.integer(unlist(numbers)))
    
    # Combine sorted numbers with "+"
    transformed_code <- paste("PCB", paste(sorted_numbers, collapse = "+"), sep = "")
  } else {
    transformed_code <- chem_code
  }
  
  return(transformed_code)
}

# Apply function to change the PCB names
NC_data$Transformed_chemical_name <- sapply(NC_data$chemical_name,
                                              transform_chemical_name)

# Create a new data frame with transposed values
transposed_data <- NC_data %>%
  pivot_wider(
    id_cols = c(sys_sample_code, sample_date, matrix_code, fraction,
                analytic_method, y_coord, x_coord, target_unit),
    names_from = Transformed_chemical_name,
    values_from = result_value
  )

# Multiply PCB columns by 1000
pcb_columns <- names(transposed_data)[grep("^PCB\\d+", names(transposed_data))]
transposed_data[, pcb_columns] <- transposed_data[, pcb_columns] * 1000

# Update target_unit to "pg/l"
transposed_data$target_unit <- "pg/l"


# Here!


# Read the JSON file with new congener list from code NewPCBList.R
pcb_groups <- read_json("Data/pcb_groups.json")

# Create an empty data frame to store the grouped data
grouped_data <- NC_data %>%
  distinct(sys_sample_code, sample_date, matrix_code, fraction,
           analytic_method, y_coord, x_coord, target_unit, dilution_factor)

# Remove spaces and special characters from column names in transposed_data
colnames(transposed_data) <- gsub("[^[:alnum:]]", "", colnames(transposed_data))

# Iterate through the pcb_groups list and sum columns
for (group_name in names(pcb_groups)) {
  group_columns <- pcb_groups[[group_name]]
  group_columns <- gsub("[^[:alnum:]]", "", group_columns)  # Clean group column names
  grouped_data[group_name] <- rowSums(transposed_data[group_columns], na.rm = TRUE)
}

# Remove SAMPLE_TYPE_CODE and MATRIX_CODE columns from grouped_data
grouped_data <- grouped_data %>%
  select(-SAMPLE_TYPE_CODE, -MATRIX_CODE, -ANALYTIC_METHOD)

# Multiply PCB columns by 1000
pcb_columns <- names(grouped_data)[grep("^PCB\\d+", names(grouped_data))]
grouped_data[, pcb_columns] <- grouped_data[, pcb_columns] * 1000

# Update RESULT_UNIT to "pg/l"
grouped_data$REPORT_RESULT_UNIT <- "pg/l"




# Create a new column named "tPCB" that sums columns 6 to 109
grouped_data <- grouped_data %>%
  mutate(tPCB = rowSums(select(., starts_with("PCB")), na.rm = TRUE))

# Replace 0s with NA for columns 6 to 109
grouped_data <- grouped_data %>%
  mutate_at(vars(starts_with("PCB")), ~ ifelse(. == 0, NA, .))

# Change SAMPLE_NAME
grouped_data <- grouped_data %>%
  mutate(SAMPLE_NAME = substr(SAMPLE_NAME, 1, nchar(SAMPLE_NAME) - 3))

# Remove SAMPLE_NAME with no name (NA)
grouped_data <- grouped_data %>%
  filter(!is.na(SAMPLE_NAME))

# Change the name of columns to be consistent
colnames(grouped_data)[colnames(grouped_data) == "SAMPLE_DATE"] <- "SampleDate"
colnames(grouped_data)[colnames(grouped_data) == "Y_COORD"] <- "Latitude"
colnames(grouped_data)[colnames(grouped_data) == "X_COORD"] <- "Longitude"
colnames(grouped_data)[colnames(grouped_data) == "REPORT_RESULT_UNIT"] <- "Units"

# Export results
write.csv(grouped_data, file = "Data/NewtownCreek/ntc.csv")

