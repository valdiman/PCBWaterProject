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
LMMB_data <- read_csv("Data/LMMB/Tributaries/1994Tributary.csv")

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
         QC_TYPE == "routine field sample",
         SAMPLE_TYPE == "Composite",
         grepl("^PCB", ANL_CODE))

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

# Working on ordering and few odd PCB combinations
# Need to change the PCB order
custom_order <- c("PCB3", "PCB4+10", "PCB5+8", "PCB6", "PCB7", "PCB7+9",
                  "PCB15+17", "PCB16+32", "PCB17", "PCB18", "PCB19", "PCB22",
                  "PCB24+27", "PCB25", "PCB26", "PCB28+31", "PCB33", "PCB37+42",
                  "PCB40", "PCB41+64+71", "PCB44", "PCB45", "PCB46", "PCB47+48",
                  "PCB49", "PCB51", "PCB52", "PCB53", "PCB56+60", "PCB63", "PCB66",
                  "PCB66+95", "PCB70+76", "PCB74", "PCB77+110", "PCB82", "PCB83",
                  "PCB84+92", "PCB85", "PCB87", "PCB89", "PCB91", "PCB95", "PCB97",
                  "PCB99", "PCB101", "PCB105+132+153", "PCB118", "PCB123+149",
                  "PCB128", "PCB132+153", "PCB135+144", "PCB136", "PCB137+176",
                  "PCB138+163", "PCB141", "PCB146", "PCB149", "PCB151", "PCB158",
                  "PCB167", "PCB170+190", "PCB171+202", "PCB172", "PCB172+197",
                  "PCB174", "PCB177", "PCB178", "PCB180", "PCB182+187", "PCB183",
                  "PCB185", "PCB193", "PCB194", "PCB195+208", "PCB196+203", "PCB198",
                  "PCB199", "PCB201", "PCB206", "PCB207", "tPCB")

# Reorder the data frame based on the custom order
transposed_data <- transposed_data %>%
  select(all_of(1:5), all_of(custom_order))

# Change NA to 0s
transposed_data[is.na(transposed_data)] <- 0

# Sum PCB7 and PCB7+9 and create a new column PCB7+9
transposed_data$PCB7_plus_9 <- transposed_data$PCB7 + transposed_data$`PCB7+9`

# Remove the original PCB7 and PCB7+9 columns
transposed_data <- transposed_data[, !colnames(transposed_data) %in% c("PCB7", "PCB7+9")]

# Rename the new column
colnames(transposed_data)[colnames(transposed_data) == "PCB7_plus_9"] <- "PCB7+9"

# Reorder columns to place "PCB7+9" at position 10
transposed_data <- transposed_data %>% 
  select(1:9, "PCB7+9", 10:ncol(transposed_data))

# Modify PCB15 and PCB17
# Create new columns newPCB15 and newPCB17
transposed_data$newPCB15 <- transposed_data$'PCB15+17' * 0.5
transposed_data$newPCB17 <- transposed_data$'PCB15+17' * 0.5

# Add the values of newPCB17 to PCB17
transposed_data$PCB17 <- transposed_data$PCB17 + transposed_data$newPCB17

# Remove the original PCB15+17 and newPCB17 columns
transposed_data <- transposed_data[, !colnames(transposed_data) %in% c("PCB15+17", "newPCB17")]

# Rename newPCB15 to PCB15
colnames(transposed_data)[colnames(transposed_data) == "newPCB15"] <- "PCB15"

# Reorder columns to place "PCB15" at position 11
transposed_data <- transposed_data %>% 
  select(1:10, "PCB15", 11:ncol(transposed_data))

# Reorder columns to place "PCB17" at position 13
transposed_data <- transposed_data %>% 
  select(1:12, "PCB17", 13:ncol(transposed_data))

# Change PCB18 to PCB18+30
colnames(transposed_data)[colnames(transposed_data) == "PCB18"] <- "PCB18+30"

# Create PCB20+21+28+31+33+50+53
transposed_data$`PCB20+21+28+31+33+50+53` <- transposed_data$`PCB28+31` + transposed_data$PCB33 + transposed_data$PCB53

# Remove the original PCB28+31, PCB33 and PCB53 columns
transposed_data <- transposed_data[, !colnames(transposed_data) %in% c("PCB28+31", "PCB33", "PCB53")]

# Reorder columns to place "PCB20+21+28+31+33+50+53" at position 16
transposed_data <- transposed_data %>% 
  select(1:15, "PCB20+21+28+31+33+50+53", 16:ncol(transposed_data))

# Change PCB26 to PCB26+29
colnames(transposed_data)[colnames(transposed_data) == "PCB26"] <- "PCB26+29"

# Create PCB40+41+64+71+72
transposed_data$`PCB40+41+64+71+72` <- transposed_data$`PCB41+64+71` + transposed_data$PCB40

# Remove the original PCB28+31, PCB33 and PCB53 columns
transposed_data <- transposed_data[, !colnames(transposed_data) %in% c("PCB41+64+71", "PCB40")]

# Reorder columns to place "PCB40+41+64+71+72" at position 22
transposed_data <- transposed_data %>% 
  select(1:21, "PCB40+41+64+71+72", 22:ncol(transposed_data))

# Create PCB43+49+52+69+73
transposed_data$`PCB43+49+52+69+73` <- transposed_data$PCB49 + transposed_data$PCB52

# Remove the original PCB28+31, PCB33 and PCB53 columns
transposed_data <- transposed_data[, !colnames(transposed_data) %in% c("PCB49", "PCB52")]

# Reorder columns to place "PCB43+49+52+69+73" at position 23
transposed_data <- transposed_data %>% 
  select(1:22, "PCB43+49+52+69+73", 23:ncol(transposed_data))

# Change PCB44 to PCB44+47+65
colnames(transposed_data)[colnames(transposed_data) == "PCB44"] <- "PCB44+47+65"

# Create PCB45+51
transposed_data$`PCB45+51` <- transposed_data$PCB45 + transposed_data$PCB51

# Remove the original PCB45 and PCB51 columns
transposed_data <- transposed_data[, !colnames(transposed_data) %in% c("PCB45", "PCB51")]

# Reorder columns to place "PCB45+51" at position 25
transposed_data <- transposed_data %>% 
  select(1:24, "PCB45+51", 25:ncol(transposed_data))

# Change PCB47+48 to PCB48+59+62+75
colnames(transposed_data)[colnames(transposed_data) == "PCB47+48"] <- "PCB48+59+62+75"

# Modify PCB66 and PCB66+95
# Add the values of newPCB66+95
transposed_data$`newPCB66+95` <- transposed_data$PCB66 + transposed_data$`PCB66+95` + transposed_data$PCB95

# Remove the original PCB66, PCB66+95 and PCB95 columns
transposed_data <- transposed_data[, !colnames(transposed_data) %in% c("PCB66", "PCB66+95", "PCB95")]

# Rename newPCB66+95 to PCB66+95
colnames(transposed_data)[colnames(transposed_data) == "newPCB66+95"] <- "PCB66+95"

# Reorder columns to place "PCB66+95" at position 35
transposed_data <- transposed_data %>% 
  select(1:34, "PCB66+95", 35:ncol(transposed_data))

# Create PCB61+66+70+74+76+93+95+98+100+102
transposed_data$`PCB61+66+70+74+76+93+95+98+100+102` <- transposed_data$`PCB66+95` + transposed_data$`PCB70+76` + transposed_data$PCB74

# Remove the original PCB66+95, PCB70+76,  and PCB74 columns
transposed_data <- transposed_data[, !colnames(transposed_data) %in% c("PCB66+95", "PCB70+76", "PCB74")]

# Reorder columns to place "PCB61+66+70+74+76+93+95+98+100+102" at position 29
transposed_data <- transposed_data %>% 
  select(1:28, "PCB61+66+70+74+76+93+95+98+100+102", 29:ncol(transposed_data))

# Create PCB77+85+110+111+115+116+117
transposed_data$`PCB77+85+110+111+115+116+117` <- transposed_data$`PCB77+110` + transposed_data$PCB85

# Remove the original PCB77+110, and PCB85 columns
transposed_data <- transposed_data[, !colnames(transposed_data) %in% c("PCB77+110", "PCB85")]

# Reorder columns to place "PCB77+85+110+111+115+116+117" at position 31
transposed_data <- transposed_data %>% 
  select(1:30, "PCB77+85+110+111+115+116+117", 31:ncol(transposed_data))









# Modify PCB105+132+153 and PCB132+153, and create a new PCB105
# Rename PCB105+132+153 to PCB105
colnames(transposed_data)[colnames(transposed_data) == "PCB105+132+153"] <- "PCB105"

# Sum PCB172 and PCB172+197 and create a new column PCB172+197
transposed_data$PCB172_plus_197 <- transposed_data$PCB172 + transposed_data$`PCB172+197`

# Remove the original PCB172 and PCB172+197 columns
transposed_data <- transposed_data[, !colnames(transposed_data) %in% c("PCB172", "PCB172+197")]

# Rename the new column
colnames(transposed_data)[colnames(transposed_data) == "PCB172_plus_197"] <- "PCB172+197"

# Reorder columns to place "PCB172+197" at position 65
transposed_data <- transposed_data %>% 
  select(1:64, "PCB172+197", 65:ncol(transposed_data))



# Export results
write.csv(transposed_data, file = "Data/LMMB/Tributaries/1994T.csv")
## need to continue working on excel. The congener list is too complicated
## to be fixed here. 



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
         



