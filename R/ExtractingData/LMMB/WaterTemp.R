# Code to read and change format of water temperatures for Lake Michigan 93-95
# https://www.ndbc.noaa.gov/station_history.php?station=45002

# Read the specific columns (YY, MM, DD, WTMP)
# Temperature in C
WT1993 <- read.table('Data/LMMB/OpenLake/WaterTemp/45002h1993.txt',
                     header = TRUE,
                     colClasses = c("numeric", "numeric",
                                    "numeric",
                                    "numeric"))[, c("YY", "MM", "DD", "WTMP")]

WT1994 <- read.table('Data/LMMB/OpenLake/WaterTemp/45002h1994.txt',
                     header = TRUE,
                     colClasses = c("numeric", "numeric",
                                    "numeric",
                                    "numeric"))[, c("YY", "MM", "DD", "WTMP")]

WT1995 <- read.table('Data/LMMB/OpenLake/WaterTemp/45002h1995.txt',
                     header = TRUE,
                     colClasses = c("numeric", "numeric",
                                    "numeric",
                                    "numeric"))[, c("YY", "MM", "DD", "WTMP")]

# Daily average, not including 999s
WT1993 <- aggregate(WTMP ~ YY + MM + DD, data = WT1993, mean, subset = WTMP != 999)
WT1994 <- aggregate(WTMP ~ YY + MM + DD, data = WT1994, mean, subset = WTMP != 999)
WT1995 <- aggregate(WTMP ~ YY + MM + DD, data = WT1995, mean, subset = WTMP != 999)

# Transforming temperature from C to K
WT1993$WTMP_K <- WT1993$WTMP + 273.15
WT1994$WTMP_K <- WT1994$WTMP + 273.15
WT1995$WTMP_K <- WT1995$WTMP + 273.15

# Create a new column with date
WT1993$Date <- as.Date(paste(WT1993$MM, WT1993$DD, WT1993$YY,
                             sep = "/"), format = "%m/%d/%y")
WT1994$Date <- as.Date(paste(WT1994$MM, WT1994$DD, WT1994$YY,
                             sep = "/"), format = "%m/%d/%y")
WT1995$Date <- as.Date(paste(WT1995$MM, WT1995$DD, WT1995$YY,
                             sep = "/"), format = "%m/%d/%y")

# Export results
write.csv(WT1993,
          file = "Output/Data/Sites/csv/GreatLakes/WT1993.csv")
write.csv(WT1994,
          file = "Output/Data/Sites/csv/GreatLakes/WT1994.csv")
write.csv(WT1995,
          file = "Output/Data/Sites/csv/GreatLakes/WT1995.csv")

