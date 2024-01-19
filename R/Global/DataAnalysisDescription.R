## Water PCB concentrations analysis.
## Data were obtained from EPA and contractors from PCB Superfund
## sites in USA. Using log10 of the sum of PCB.

# Install packages
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("robustbase")
install.packages("dplyr")
install.packages("tibble")
install.packages("Matrix")
install.packages("zoo")
install.packages("reshape")
install.packages("sf")
install.packages("sfheaders")

# Load libraries
{
  library(ggplot2)
  library(scales) # function trans_breaks
  library(stringr) # str_detect
  library(robustbase) # function colMedians
  library(dplyr) # performs %>%
  library(tibble) # adds a column
  library(zoo) # yields seasons
  library(reshape)
  library(sf)
  library(sfheaders) # Create file to be used in Google Earth
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# General information -----------------------------------------------------
# Number of locations and number of site per location
location_count <- wdc %>%
  group_by(LocationName) %>%
  summarise(count = n())

print(location_count)
# Median amount of samples per location
median(location_count$count)

# Number of locations per states
state_count <- wdc %>%
  group_by(StateSampled) %>%
  summarise(count = n())

print(state_count)

# Number of sites and number of replicates per site
site_count <- wdc %>%
  group_by(SiteID) %>%
  summarise(count = n())

print(site_count)

# Media of number of site available
median(site_count$count)

# Number of replicates per sites from the same day
site_repli_count <- wdc %>%
  group_by(SampleID) %>%
  summarise(count = n())

print(site_repli_count)

# Media of number of replicates per sites from the same day
median(site_repli_count$count)

# Find the SampleID with the highest count
max_count_sample <- site_repli_count %>%
  filter(count == max(count)) %>%
  pull(SampleID)

filtered_wdc <- wdc %>%
  filter(SampleID == max_count_sample)

# Extract Site Name
extracted_string <- filtered_wdc %>%
  select(SiteName) %>%
  pull()

# Display Site Name
print(extracted_string[1])

# Aroclor Congener Summary ------------------------------------------------
# Create a new data frame with 'Year' and 'AroclorCongener' columns
wdc.aroclorcongener <- wdc %>%
  mutate(Year = as.integer(format(as.Date(SampleDate, format = "%m/%d/%y"), "%Y"))) %>%
  filter(AroclorCongener %in% c("Congener", "Aroclor")) %>%
  select(Year, AroclorCongener)

# Calculate counts for 'Congener' and 'Aroclor' for each year
yearly_counts <- wdc.aroclorcongener %>%
  group_by(Year, AroclorCongener) %>%
  summarize(Count = n())

# Calculate percentages for 'Congener' and 'Aroclor' for each year
percentages <- yearly_counts %>%
  group_by(Year) %>%
  mutate(Percent = Count / sum(Count))

# Precompute the percentages for 'Congener' and 'Aroclor' for each year
precomputed_percentages <- percentages %>%
  group_by(Year, AroclorCongener) %>%
  summarize(Percent = sum(Percent))

# Calculate the number of samples per year for 'Congener' and 'Aroclor'
sample_counts <- yearly_counts %>%
  group_by(Year) %>%
  summarize(TotalSamples = sum(Count))

# Convert Year to character
precomputed_percentages$Year <- as.character(precomputed_percentages$Year)

# Create a stacked bar plot with percentages as percentages (0 to 100)
plot.aroclor.congener <- ggplot(precomputed_percentages,
                                aes(x = factor(Year,
                                               levels = unique(precomputed_percentages$Year)))) +
  geom_bar(aes(y = Percent * 100, fill = AroclorCongener), stat = "identity") +
  geom_text(data = sample_counts, aes(label = TotalSamples, y = 100), size = 2.5,
            fontface = "bold") +
  theme_classic() +
  labs(title = "Percentage of Congener and Aroclor Over the Years",
       x = "Year",
       y = "Percentage") +
  scale_fill_manual(values = c("Congener" = "grey50", "Aroclor" = "grey90")) +
  scale_y_continuous(limits = c(0, 100)) +  # Set the y-axis limits
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black")) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 11)) +
  coord_cartesian(ylim = c(0, 100))

# See the plot
print(plot.aroclor.congener)

# Save plot in folder
ggsave("Output/Plots/Global/AroclorCongenerV2.png", plot = plot.aroclor.congener,
       width = 10, height = 3, dpi = 300)

# Aroclor summary ---------------------------------------------------------
# Number of samples analyzed using Aroclor method
count_Aroclor <- sum(wdc$AroclorCongener == "Aroclor")
total_samples <- length(wdc[,1])
percent_aroclor <- count_Aroclor/total_samples*100

# Calculate sample % for each Aroclor mixtures
aroclors <- c('A1016', 'A1221', 'A1232', 'A1242', 'A1248',
              'A1254', 'A1260')

# Calculate the number of non-NA values (frequency of numbers) in each Aroclor
frequency_aroclors <- lapply(wdc[aroclors], function(column) {
  length(na.omit(column))
})

# Percentage of Aroclor mixtures in relation to all Aroclors
for (i in seq_along(aroclors)) {
  column_name <- aroclors[i]
  print(paste(column_name))
  print(frequency_aroclors[[i]]/count_Aroclor*100)
}

# Congener frequency ------------------------------------------------------
{
  # Remove metadata
  wdc.nmd <- subset(wdc, select = -c(SampleID:AroclorCongener)) #nmd = no meta data
  # Remove Aroclor data
  wdc.nmd <- subset(wdc.nmd, select = -c(A1016:A1260))
  # (2) Only consider congener data
  wdc.cong <- subset(wdc, AroclorCongener == "Congener")
  # Remove metadata
  wdc.cong.1 <- subset(wdc.cong, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data and tPCB
  wdc.cong.1 <- subset(wdc.cong.1, select = -c(A1016:tPCB))
}

# Create a frequency detection plot
{
  wdc.cong.freq <- colSums(! is.na(wdc.cong.1) & (wdc.cong.1 !=0))/nrow(wdc.cong.1)
  wdc.cong.freq <- data.frame(wdc.cong.freq)
  colnames(wdc.cong.freq) <- c("PCB.frequency")
  congener <- row.names(wdc.cong.freq)
  wdc.cong.freq <- cbind(congener, wdc.cong.freq$PCB.frequency)
  colnames(wdc.cong.freq) <- c("congener", "PCB.frequency")
  wdc.cong.freq <- data.frame(wdc.cong.freq)
  wdc.cong.freq$congener <- as.character(wdc.cong.freq$congener)
  wdc.cong.freq$congener <- gsub('\\.', '+', wdc.cong.freq$congener) # replace dot for +
  wdc.cong.freq$PCB.frequency <- as.numeric(as.character(wdc.cong.freq$PCB.frequency))
  wdc.cong.freq$congener <- factor(wdc.cong.freq$congener,
                                   levels = rev(wdc.cong.freq$congener)) # change the order to be plotted.
}

# Summary statistic of frequency of detection
summary(wdc.cong.freq$PCB.frequency)

# Frequency detection plot
plot.cong.freq <- ggplot(wdc.cong.freq, aes(x = 100*PCB.frequency, y = congener)) +
  geom_bar(stat = "identity", fill = "#66ccff", color = "black") +
  geom_vline(xintercept = 100*mean(wdc.cong.freq$PCB.frequency),
             color = "red") +
  ylab("") +
  theme_bw() +
  xlim(c(0,100)) +
  theme(aspect.ratio = 20/5) +
  xlab(expression(bold("Frequency detection (%)"))) +
  theme(axis.text.x = element_text(face = "bold", size = 8),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.text.y = element_text(face = "bold", size = 7))

# See plot
print(plot.cong.freq)

# Save map in folder
ggsave("Output/Plots/Global/FreqPCBV02.png", plot = plot.cong.freq,
       width = 5, height = 10, dpi = 300)

# Total PCB description ---------------------------------------------------
# Data preparation
{
  # Change date format
  SampleDate <- as.Date(wdc$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Create individual code for each site sampled
  site.numb <- wdc$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(wdc$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  tpcb <- cbind(factor(wdc$SiteID), SampleDate,
                wdc$Latitude, wdc$Longitude, wdc$tPCB,
                data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                      "tPCB", "time", "site.code", "season")
}

# Get coordinates per site to plot in Google Earth
location <- tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = location, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Global/PCBSampleLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Summary statistic of total PCB (congeners + Aroclor) in pg/L not including 0s
summary(wdc$tPCB)

# Highest sample location and time
max_sample <- wdc[which.max(wdc$tPCB), ]
max_sample <- max_sample[c("LocationName", "SampleDate", "SiteName")]
print(max_sample)

# Individual congeners description
summary(wdc.cong.1, na.rm = T, zero = T)
# Get the max value for each congener
cong.max <-as.numeric(sub('.*:', '',
                          summary(wdc.cong.1, na.rm = T,
                                  zero = T)[6,]))
# Add congener
cong.max <- cbind(congener, data.frame(cong.max))

# Obtain the median for each individual congener
cong.median <- as.numeric(sub('.*:',
                              '', summary(wdc.cong.1, na.rm = T,
                                          zero = T)[3,]))
# Add congener
cong.median <- cbind(congener, data.frame(cong.median))
# Min
print(min(cong.median$cong.median))
#Max
print(max(cong.median$cong.median))
# Mean
print(mean(cong.median$cong.median))

# Global plots ------------------------------------------------------------
# (1) Histogram
hist(tpcb$tPCB)
hist(log10(tpcb$tPCB))

## (2) Total PCBs in 1 box plot
## include 64 and 640 pg/L from EPA
plot.box.tPCB <- ggplot(tpcb, aes(x = "", y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "white") +
  geom_boxplot(lwd = 0.8, width = 0.7, outlier.shape = NA, alpha = 0) +
  theme_classic() +
  theme(aspect.ratio = 14/2) +
  xlab(expression(bold(Sigma*"PCB")))+
  ylab(expression(bold("Water Concentration 1979 - 2020 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 16),
        axis.title.y = element_text(face = "bold", size = 16)) +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, vjust = 5)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_hline(yintercept = 640, color = "#9999CC",
             linewidth = 0.8) +
  geom_hline(yintercept = 64, color = "#CC6666",
             linewidth = 0.8)

# See plot
print(plot.box.tPCB)

# Save map in folder
ggsave("Output/Plots/Global/tPCBBoxPlotV03.png", plot = plot.box.tPCB,
       width = 5, height = 10, dpi = 300)

# Calculate % samples above both EPA thresholds
EPA640 <- sum(tpcb$tPCB > 640)/nrow(tpcb)*100
EPA64 <- sum(tpcb$tPCB > 64)/nrow(tpcb)*100

# (3) Total PCBs for selected locations
sites_to_include <- c("Housatonic River", "New Bedford Harbor", "Passaic River",
                      "Hudson River", "Kalamazoo River", "Fox River",
                      "Portland Harbor", "Lake Michigan Mass Balance",
                      "Spokane River", "Chesapeake Bay")

# Filter the data to include only the specified sites
filtered_data <- wdc %>%
  filter(LocationName %in% sites_to_include)

tpcb.site <- ggplot(filtered_data, aes(x = factor(LocationName),
                                       y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 20/15) +
  ylab(expression(bold("Water Concentration " *Sigma*"PCB (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 12,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "white") +
  geom_boxplot(lwd = 0.5, width = 0.7, outlier.shape = NA, alpha = 0)

# See plot
print(tpcb.site)

# Save plot in folder
ggsave("Output/Plots/Global/tPCBSiteV02.png", plot = tpcb.site,
       width = 5, height = 10, dpi = 300)

# (4) Box plot for individual PCBs
PCBi_boxplot <- ggplot(stack(wdc.cong.1), aes(x = ind, y = values)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot(width = 0.6, shape = 21, outlier.fill = "white",
               fill = "white", outlier.shape = 21) +
  scale_x_discrete(labels = wdc.cong.freq$congener) + # use to change the "." to "+"
  theme_bw() +
  theme(aspect.ratio = 25/135) +
  xlab(expression("")) +
  ylab(expression(bold("PCB congener concentration (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8,
                                   color = "black"),
        axis.title.y = element_text(face = "bold", size = 8,
                                    color = "black")) +
  theme(axis.text.x = element_text(face = "bold", size = 6,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.6, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# See plot
print(PCBi_boxplot)

# Save map in folder
ggsave("Output/Plots/Global/PCBiBoxPlotV03.png", plot = PCBi_boxplot,
       width = 10, height = 5, dpi = 300)

# (5) Individual PCBs for selected locations
# Filter out rows with NA and 0 values in the 'PCBi' column
# (5.1) PCB11 plot
filtered_datai <- wdc %>%
  filter(LocationName %in% sites_to_include, !is.na(PCB11),
         !(PCB11 == 0))

pcbi.site <- ggplot(filtered_datai, aes(x = factor(LocationName),
                                        y = PCB11)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 20/15) +
  ylab(expression(bold("PCB 11 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 14)) +
  theme(axis.text.x = element_text(face = "bold", size = 12,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "white") +
  geom_boxplot(lwd = 0.5, width = 0.7, outlier.shape = NA, alpha = 0)

# See plot
print(pcbi.site)

# Save plot in folder
ggsave("Output/Plots/Global/PCB11Site.png", plot = pcbi.site,
       width = 5, height = 10, dpi = 300)

# (5.2) PCB20.21.28.31.33.50.53 plot
filtered_datai <- wdc %>%
  filter(LocationName %in% sites_to_include, !is.na(PCB20.21.28.31.33.50.53),
         !(PCB20.21.28.31.33.50.53 == 0))

pcbi.site <- ggplot(filtered_datai, aes(x = factor(LocationName),
                                        y = PCB20.21.28.31.33.50.53)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 20/15) +
  ylab(expression(bold("PCB 20+21+28+31+33+50+53 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 14,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "white") +
  geom_boxplot(lwd = 0.5, width = 0.7, outlier.shape = NA, alpha = 0)

# See plot
print(pcbi.site)

# Save plot in folder
ggsave("Output/Plots/Global/PCB20Site.png", plot = pcbi.site,
       width = 5, height = 10, dpi = 300)

# (5.3) PCBB44+47+65 plot
filtered_datai <- wdc %>%
  filter(LocationName %in% sites_to_include, !is.na(PCB44.47.65),
         !(PCB44.47.65 == 0))

pcbi.site <- ggplot(filtered_datai, aes(x = factor(LocationName),
                                        y = PCB44.47.65)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 20/15) +
  ylab(expression(bold("PCB 44+47+65 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 14,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "white") +
  geom_boxplot(lwd = 0.5, width = 0.7, outlier.shape = NA, alpha = 0)

# See plot
print(pcbi.site)

# Save plot in folder
ggsave("Output/Plots/Global/PCB44Site.png", plot = pcbi.site,
       width = 5, height = 10, dpi = 300)

# (5.4) PCBB67 plot
filtered_datai <- wdc %>%
  filter(LocationName %in% sites_to_include, !is.na(PCB67),
         !(PCB67 == 0))

pcbi.site <- ggplot(filtered_datai, aes(x = factor(LocationName),
                                        y = PCB67)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 20/15) +
  ylab(expression(bold("PCB 67 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 14,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "white") +
  geom_boxplot(lwd = 0.5, width = 0.7, outlier.shape = NA, alpha = 0)

# See plot
print(pcbi.site)

# Save plot in folder
ggsave("Output/Plots/Global/PCB67Site.png", plot = pcbi.site,
       width = 5, height = 10, dpi = 300)

# (6) Time trend plots
plot.time.tPCB <- ggplot(tpcb, aes(y = tPCB,
                                   x = format(date,'%Y'))) +
  geom_point(shape = 21, size = 2.5, fill = "white") +
  xlab("") +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB 1979 - 2020 (pg/L)"))))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  annotation_logticks(sides = "l") +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black")) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 11))

# See plot
print(plot.time.tPCB)

# Save plot in folder
ggsave("Output/Plots/Global/tPCBTimeV03.png", plot = plot.time.tPCB,
       width = 10, height = 3, dpi = 300)

# (7) Individual PCB trend plots
# (7.1) PCB5.8
pcb5.8 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB5.8,
                data.frame(time.day), site.numb, season.s)
# Add column names
colnames(pcb5.8) <- c("SiteID", "date", "PCB5.8", "time",
                      "site.code", "season")
# Remove 0s and NA values
pcb5.8 <- pcb5.8[complete.cases(pcb5.8$PCB5.8) & pcb5.8$PCB5.8 != 0, ]

plot.time.pcb5.8 <- ggplot(pcb5.8, aes(y = PCB5.8,
                                     x = format(date,'%Y'))) +
  geom_point(shape = 21, size = 2.5, fill = "white") +
  xlab("") +
  ylab(expression(bold("Water Concentration PCB 5+8 1994 - 2019 (pg/L)"))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  annotation_logticks(sides = "l") +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black")) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 11))

# See plot
print(plot.time.pcb5.8)

# (7.2) PCB11
pcb11 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB11,
               data.frame(time.day), site.numb, season.s)
# Add column names
colnames(pcb11) <- c("SiteID", "date", "PCB11", "time",
                     "site.code", "season")
# Remove 0s and NA values
pcb11 <- pcb11[complete.cases(pcb11$PCB11) & pcb11$PCB11 != 0, ]

plot.time.pcb11 <- ggplot(pcb11, aes(y = PCB11,
                                   x = format(date,'%Y'))) +
  geom_point(shape = 21, size = 2.5, fill = "white") +
  xlab("") +
  ylab(expression(bold("Water Concentration PCB11 2006 - 2019 (pg/L)"))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  annotation_logticks(sides = "l") +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black")) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 11))

# See plot
print(plot.time.pcb11)

# (7.3) PCB18.30
pcb18.30 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB18.30,
                  data.frame(time.day), site.numb, season.s)
# Add column names
colnames(pcb18.30) <- c("SiteID", "date", "PCB18.30", "time",
                        "site.code", "season")
# Remove 0s and NA values
pcb18.30 <- pcb18.30[complete.cases(pcb18.30$PCB18.30) & pcb18.30$PCB18.30 != 0, ]

plot.time.pcb18.30 <- ggplot(pcb18.30, aes(y = PCB18.30,
                                     x = format(date,'%Y'))) +
  geom_point(shape = 21, size = 2.5, fill = "white") +
  xlab("") +
  ylab(expression(bold("Water Concentration PCB 18+30 1994 - 2019 (pg/L)"))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  annotation_logticks(sides = "l") +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black")) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 11))

# See plot
print(plot.time.pcb18.30)

# (7.4) PCB20.21.28.31.33.50.53
pcb20 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB20.21.28.31.33.50.53,
               data.frame(time.day), site.numb, season.s)
# Add column names
colnames(pcb20) <- c("SiteID", "date", "PCB20", "time",
                     "site.code", "season")
# Remove 0s and NA values
pcb20 <- pcb20[complete.cases(pcb20$PCB20) & pcb20$PCB20 != 0, ]

plot.time.pcb20 <- ggplot(pcb20, aes(y = PCB20,
                                     x = format(date,'%Y'))) +
  geom_point(shape = 21, size = 2.5, fill = "white") +
  xlab("") +
  ylab(expression(bold("Water Concentration PCB 20.. 1994 - 2019 (pg/L)"))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  annotation_logticks(sides = "l") +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black")) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 11))

# See plot
print(plot.time.pcb20)

# (7.5) PCB44+47+65
pcb44 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB44.47.65,
               data.frame(time.day), site.numb, season.s)
# Add column names
colnames(pcb44) <- c("SiteID", "date", "PCB44", "time",
                     "site.code", "season")
# Remove 0s and NA values
pcb44 <- pcb44[complete.cases(pcb44$PCB44) & pcb44$PCB44 != 0, ]

plot.time.pcb44 <- ggplot(pcb44, aes(y = PCB44,
                                     x = format(date,'%Y'))) +
  geom_point(shape = 21, size = 2.5, fill = "white") +
  xlab("") +
  ylab(expression(bold("Water Concentration PCB 44.. 1994 - 2019 (pg/L)"))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  annotation_logticks(sides = "l") +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black")) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 11))

# See plot
print(plot.time.pcb44)

# (7.6) PCB 67
pcb67 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB67,
               data.frame(time.day), site.numb, season.s)
# Add column names
colnames(pcb67) <- c("SiteID", "date", "PCB67", "time",
                     "site.code", "season")
# Remove 0s and NA values
pcb67 <- pcb67[complete.cases(pcb67$PCB67) & pcb67$PCB67 != 0, ]

plot.time.pcb67 <- ggplot(pcb67, aes(y = PCB67,
                                     x = format(date,'%Y'))) +
  geom_point(shape = 21, size = 2.5, fill = "white") +
  xlab("") +
  ylab(expression(bold("Water Concentration PCB 67 2004 - 2019 (pg/L)"))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  annotation_logticks(sides = "l") +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black")) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 11))

# See plot
print(plot.time.pcb67)

# (7.7) PCB 106+118
pcb106.118 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB106.118,
                    data.frame(time.day), site.numb, season.s)
# Add column names
colnames(pcb106.118) <- c("SiteID", "date", "PCB106.118", "time",
                          "site.code", "season")
# Remove 0s and NA values
pcb106.118 <- pcb106.118[complete.cases(pcb106.118$PCB106.118) & pcb106.118$PCB106.118 != 0, ]

plot.time.pcb106.118 <- ggplot(pcb106.118, aes(y = PCB106.118,
                                     x = format(date,'%Y'))) +
  geom_point(shape = 21, size = 2.5, fill = "white") +
  xlab("") +
  ylab(expression(bold("Water Concentration PCB 106+118 2004 - 2019 (pg/L)"))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  annotation_logticks(sides = "l") +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black")) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 11))

# See plot
print(plot.time.pcb106.118)

# Extra plots -------------------------------------------------------------
# (7) Seasonality
ggplot(tpcb, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB 1979 - 2020 (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

