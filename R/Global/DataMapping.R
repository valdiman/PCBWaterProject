## Water PCB concentrations mapping.

# Install packages
install.packages("ggplot2")
install.packages("devtools")
install.packages("dplyr")
install.packages("stringr")
install.packages("maps")
install.packages("mapdata")
install.packages("ggmap")
install.packages("usethis")
install.packages("GISTools")
install.packages("rgeos")
install.packages("ggsn")
install.packages("ggrepel")
install.packages("ggpp")
install.packages("scales")
install.packages("viridis")
install.packages("tidyr")

# Load libraries
{
  library(dplyr)
  library(usethis)
  library(devtools)
  library(ggplot2)
  library(ggmap) # function map_data
  library(maps)
  library(leaflet)
  library(rgeos)
  library(ggsn)
  library(ggrepel)
  library(reshape2)
  library(ggpmisc)
  library(scales) # add commas in legend in maps
  library(cowplot)
  library(viridis) # customize color legend
  library(tidyr)
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")
# Extract sample site locations -------------------------------------------
# Average tPCB per sample site
tpcb.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = wdc, mean)

# USA/State maps -------------------------------------------------------------
us <- map_data("usa")
states <- map_data("state")
counties <- map_data("county")

# Find number of samples per StateSampled
wdc.2 <- wdc %>%
  group_by(StateSampled) %>%
  summarise(n = n())

# Convert wdc.2 to long format
wdc.3 <- gather(wdc.2, key = "Variable", value = "Value", -StateSampled)

# Find number of samples per LocationName
wdc.4 <- wdc %>%
  group_by(LocationName) %>%
  summarise(n = n())

# Convert wdc.4 to long format
wdc.5 <- gather(wdc.4, key = "Variable", value = "Value", -LocationName)

# (1) Map of US with locations
maploc <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = "transparent") +
  coord_fixed(1.3) +
  theme_void() +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_polygon(color = "black", fill = NA) +
  geom_point(data = tpcb.ave, aes(x = Longitude, y = Latitude),
             color = "black",
             size = 2.2, shape = 20) +
  geom_text(data = wdc.3, aes(x = -66.5, y = 50.2 - seq_along(StateSampled),
                              label = paste(StateSampled, Value), hjust = 0,
                              vjust = 1), size = 3) +
  geom_text(data = wdc.5, aes(x = -138, y = 50.2 - seq_along(LocationName),
                              label = paste(LocationName, Value), hjust = 0,
                              vjust = 1), size = 3) +
  geom_text(aes(x = -69.5, y = 50.5, label = "States/Samples", hjust = 0,
                vjust = 1), size = 3, fontface = "bold") +
  geom_text(aes(x = -138, y = 50.5, label = "Location/Samples", hjust = 0,
                vjust = 1), size = 3, fontface = "bold")

print(maploc)

# Save map in folder
ggsave("Output/Maps/Global/maplocV05.pdf", plot = maploc,
       width = 12, height = 6, dpi = 300)

# (2) Map + tPCB
maptPCB <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_polygon(data = counties, aes(x = long, y = lat, group = group),
               color = "grey", fill = NA) +  # County boundaries
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = tpcb.ave, aes(x = Longitude, y = Latitude,
                                  fill = tPCB), alpha = 1, color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(
    name = element_blank(),
    limits = c(1, 40000000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position = c(1.1, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18)  # Adjust the size of legend title
  )

print(maptPCB)  # Print the plot

# Save the plot as PDF
ggsave("Output/Maps/Global/maptPCBV06.pdf", plot = maptPCB,
       width = 14, height = 4)

# Individual PCB Maps -----------------------------------------------------
# PCB5+8
pcb5.8 <- wdc[wdc$PCB5.8 != 0, ]
# Average PCB5+8 per sample site
pcb5.8.ave <- aggregate(PCB5.8 ~ SiteID + Latitude + Longitude,
                      data = pcb5.8, mean)

# Plot
mapPCB5.8 <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = pcb5.8.ave, aes(x = Longitude, y = Latitude,
                                  fill = PCB5.8), alpha = 1, color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(
    name = expression(bold("PCBs 5+8")),
    limits = c(0.5, 1500000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position = c(1.05, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18)  # Adjust the size of legend title
  )

print(mapPCB5.8)  # Print the plot

# Save map in folder
ggsave("Output/Maps/Global/mapPCB5_8.pdf", plot = mapPCB5.8,
       width = 14, height = 4)

# PCB11
pcb11 <- wdc[wdc$PCB11 != 0, ]

# Average PCB11 per sample site
pcb.11.ave <- aggregate(PCB11 ~ SiteID + Latitude + Longitude,
                      data = pcb11, mean)

# Plot
mapPCB11 <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = pcb.11.ave, aes(x = Longitude, y = Latitude,
                                  fill = PCB11), alpha = 1, color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(
    name = expression(bold("PCB 11")),
    limits = c(0.1, 2000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position = c(1.1, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18)  # Adjust the size of legend title
  )
        
print(mapPCB11)  # Print the plot

# Save map in folder
ggsave("Output/Maps/Global/mapPCB11.pdf", plot = mapPCB11,
       width = 14, height = 4)

# Filter out rows with 0 values for PCB 20+21+28+31+33+50+53
pcb20 <- wdc[wdc$PCB20.21.28.31.33.50.53 != 0, ]

# Average PCB11 per sample site
pcb.20.ave <- aggregate(PCB20.21.28.31.33.50.53 ~ SiteID + Latitude + Longitude,
                      data = pcb20, mean)

# Plot
mapPCB20 <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = pcb.20.ave, aes(x = Longitude, y = Latitude,
                                    fill = PCB20.21.28.31.33.50.53),
             alpha = 1, color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(name = expression(bold(atop("PCBs 20+21+28",
                                                   paste("+31+33+50+53")))),
    limits = c(1, 1000000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position = c(1.15, 0.5),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18)  # Adjust the size of legend title
  )

print(mapPCB20)  # Print the plot

# Save map in folder
ggsave("Output/Maps/Global/mapPCB20.pdf", plot = mapPCB20,
       width = 14, height = 4)

# Filter out rows with 0 values for PCB44+47+65
pcb44 <- wdc[wdc$PCB44.47.65 != 0, ]

# Average PCB11 per sample site
pcb.44.ave <- aggregate(PCB44.47.65 ~ SiteID + Latitude + Longitude,
                      data = pcb44, mean)

# Plot
mapPCB44 <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = pcb.44.ave, aes(x = Longitude, y = Latitude,
                                    fill = PCB44.47.65), alpha = 1,
             color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(
    name = expression(bold("PCBs 44+47+65")),
    limits = c(0.5, 200000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position = c(1.12, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18)  # Adjust the size of legend title
  )

print(mapPCB44)  # Print the plot

# Save map in folder
ggsave("Output/Maps/Global/mapPCB44.pdf", plot = mapPCB44,
       width = 14, height = 4)

# Filter out rows with 0 values for PCB67
pcb67 <- wdc[wdc$PCB67 != 0, ]

# Average PCB67 per sample site
pcb.67.ave <- aggregate(PCB67 ~ SiteID + Latitude + Longitude,
                      data = pcb67, mean)

# Plot
mapPCB67 <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = pcb.67.ave, aes(x = Longitude, y = Latitude,
                                    fill = PCB67), alpha = 1,
             color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(
    name = expression(bold("PCB 67")),
    limits = c(1, 500),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position = c(1.05, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18)  # Adjust the size of legend title
  )

print(mapPCB67)  # Print the plot

# Save map in folder
ggsave("Output/Maps/Global/mapPCB67.pdf", plot = mapPCB67,
       width = 14, height = 4)

# Filter out rows with 0 values for PCB106 + 118
pcb106 <- wdc[wdc$PCB106.118 != 0, ]

# Average PCB67 per sample site
pcb.106.ave <- aggregate(PCB106.118 ~ SiteID + Latitude + Longitude,
                        data = pcb106, mean)

# Plot
mapPCB106 <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = pcb.106.ave, aes(x = Longitude, y = Latitude,
                                    fill = PCB106.118), alpha = 1,
             color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(
    name = expression(bold("PCBs 106+118")),
    limits = c(0.3, 50000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position = c(1.1, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18)  # Adjust the size of legend title
  )

print(mapPCB106)  # Print the plot

# Save map in folder
ggsave("Output/Maps/Global/mapPCB106.pdf", plot = mapPCB106,
       width = 14, height = 4)

# Filter out rows with 0 values for PCB182+187
pcb182 <- wdc[wdc$PCB182.187 != 0, ]

# Average PCB182+187 per sample site
pcb.182.ave <- aggregate(PCB182.187 ~ SiteID + Latitude + Longitude,
                        data = pcb182, mean)

# Plot
mapPCB182 <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = pcb.182.ave, aes(x = Longitude, y = Latitude,
                                    fill = PCB182.187), alpha = 1,
             color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(
    name = expression(bold("PCBs 182+187")),
    limits = c(0.2, 20000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position = c(1.1, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18)  # Adjust the size of legend title
  )

print(mapPCB182)  # Print the plot

# Save map in folder
ggsave("Output/Maps/Global/mapPCB182.pdf", plot = mapPCB182,
       width = 14, height = 4)

# Specific locations ------------------------------------------------------
# Version 1
# Midwest
mid.west <- subset(wdc, LocationName %in% c("Fox River", "Kalamazoo River",
                                            "Lake Michigan Mass Balance",
                                            "Indiana Harbor and Ship Canal"))

mid.west.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                          data = mid.west, mean)

maptPCB <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_polygon(data = counties, aes(x = long, y = lat, group = group),
               color = "grey", fill = NA) +  # County boundaries
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = mid.west.ave, aes(x = Longitude, y = Latitude,
                                      fill = tPCB), alpha = 1, color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(
    name = element_blank(),
    limits = c(30, 1000000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3, xlim = c(-89, -84.6), ylim = c(41.5, 45.9)) +  # Adjust these values accordingly
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position = c(1.25, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold",
                              margin = margin(b = 10))  # Adjust title position
  ) +
  ggtitle("(X) Fox & Kalamazzo Rivers, LMMB & IHSC")

print(maptPCB)

# Save the plot as PDF
ggsave("Output/Maps/Global/maptPCBLMV01.pdf", plot = maptPCB,
       width = 14, height = 4)

# Housatonic Hudson & Passaic Rivers
hhp <- subset(wdc, LocationName %in% c("Housatonic River", "Hudson River",
                                       "Passaic River"))

hhp.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                     data = hhp, mean)

maptPCB <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_polygon(data = counties, aes(x = long, y = lat, group = group),
               color = "grey", fill = NA) +  # County boundaries
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = hhp.ave, aes(x = Longitude, y = Latitude,
                                 fill = tPCB), alpha = 1, color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(
    name = element_blank(),
    limits = c(40, 600000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3, xlim = c(-74.9, -70.4), ylim = c(38.8, 43.3)) +  # Adjust these values accordingly
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position = c(1.3, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold",
                              margin = margin(b = 10))  # Adjust title position
  ) +
  ggtitle("(X) Housatonic, Hudson & Passaic Rivers")

print(maptPCB)

# Save the plot as PDF
ggsave("Output/Maps/Global/maptPCBHHPV02.pdf", plot = maptPCB,
       width = 14, height = 4)

# North West
now <- subset(wdc, LocationName %in% c("Portland Harbor", "Lower Duwamish"))

now.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                     data = now, mean)

maptPCB <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_polygon(data = counties, aes(x = long, y = lat, group = group),
               color = "grey", fill = NA) +  # County boundaries
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = now.ave, aes(x = Longitude, y = Latitude,
                                 fill = tPCB), alpha = 1, color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(
    name = element_blank(),
    limits = c(80, 2000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3, xlim = c(-124.2, -121.7), ylim = c(45.5, 47.8)) +  # Adjust these values accordingly
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position = c(1.2, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold",
                              margin = margin(b = 10))  # Adjust title position
  ) +
  ggtitle("(X) Lower Duwamish River & Portland Harbor")

print(maptPCB)

# Save the plot as PDF
ggsave("Output/Maps/Global/maptPCBNOWV02.pdf", plot = maptPCB,
       width = 14, height = 4)

# Version 2
# Portland Harbor ---------------------------------------------------------
{
  # Select only from Portland Harbor
  wdc.PO <- subset(wdc, LocationName == "Portland Harbor")
  # Increase the longitude range to make the map wider
  lon_range <- 0.01  # Modify this value to control the width
  
  # Create a new bounding box with the adjusted longitude range
  PO.box <- make_bbox(
    lon = c(min(wdc.PO$Longitude) - lon_range, max(wdc.PO$Longitude) + lon_range),
    lat = wdc.PO$Latitude,
    f = 0.5)
  # Fetch the map using the new bounding box
  PO.map <- get_stamenmap(bbox = PO.box, zoom = 10)
  
  # Plot map with sites
  # Prepare data
  # Get tPCB and coordinates
  tPCB.PO <- data.frame(cbind(wdc.PO$SiteID, wdc.PO$Latitude,
                              wdc.PO$Longitude, wdc.PO$tPCB))
  # Name the columns
  colnames(tPCB.PO) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.PO$Latitude <- as.numeric(tPCB.PO$Latitude)
  tPCB.PO$Longitude <- as.numeric(tPCB.PO$Longitude)
  tPCB.PO$tPCB <- as.numeric(tPCB.PO$tPCB)
  # Average tPCB per site
  tPCB.PO.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                           data = tPCB.PO, FUN = mean)
}

# (1) Plot map + locations
ggmap(PO.map) +
  geom_point(data = tPCB.PO.ave, aes(x = Longitude, y = Latitude),
             shape = 21, color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.PO.ave, family = 'Times New Roman', size = 3, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
maptPCBPO <- ggmap(PO.map) +
  geom_point(data = tPCB.PO.ave, aes(x = Longitude, y = Latitude,
                                     size = tPCB), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotate('text', x = -122.77, y = 45.7,
           label = 'Portland Harbor (OR)', colour = 'black', size = 3.5,
           fontface = 2) +
  scale_size_area(breaks = c(100, 250, 500, 1000, 1500), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean) 2007-2019 (pg/L)"))),
                  max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(1.62, 0.75),  # Adjust the legend.position values
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))

# Save map in folder
ggsave("Output/Maps/Sites/maptPCBPortlandAveV01.png", plot = maptPCBPO,
       width = 8, height = 4, dpi = 300)

# Fox River ---------------------------------------------------------------
{
  # Select only from Fox River
  wdc.Fox <- subset(wdc, LocationName == "Fox River")
  # Create general map
  Fox.box <- make_bbox(lon = wdc.Fox$Longitude, lat = wdc.Fox$Latitude,
                       f = 0.22)
  Fox.map <- get_stamenmap(bbox = Fox.box, zoom = 10)
  
  # Plot map with sites
  # Prepare data
  # Remove samples with total PCBs  = 0
  wdc.Fox.1 <- wdc.Fox[!(wdc.Fox$tPCB) == 0, ]
  # Get tPCB and coordinates
  tPCB.Fox <- data.frame(cbind(wdc.Fox.1$SiteID, wdc.Fox.1$Latitude,
                               wdc.Fox.1$Longitude, wdc.Fox.1$tPCB))
  # Name the columns
  colnames(tPCB.Fox) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.Fox$Latitude <- as.numeric(tPCB.Fox$Latitude)
  tPCB.Fox$Longitude <- as.numeric(tPCB.Fox$Longitude)
  tPCB.Fox$tPCB <- as.numeric(tPCB.Fox$tPCB)
  # Average tPCB per site
  tPCB.Fox.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                             data = tPCB.Fox, FUN = mean)
}

# (1) Plot map + locations
ggmap(Fox.map) +
  geom_point(data = tPCB.Fox.ave, aes(x = Longitude, y = Latitude),
             shape = 21, color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.Fox.ave, family = 'Times New Roman', size = 2.8, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
maptPCBFoxRiver <- ggmap(Fox.map) +
  geom_point(data = tPCB.Fox.ave, aes(x = Longitude, y = Latitude,
                                      size = tPCB), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotate('text', x = -88.35, y = 44.6,
           label = 'Fox River/Green Bay (WI)', colour = 'black', size = 3.5,
           fontface = 2) +
  scale_size_area(breaks = c(100, 500, 1500, 30000, 65000), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean) 2005-2018 (pg/L)"))),
                  max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(1.44, 0.74),  # Adjust the legend.position values
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))

# Save map in folder
ggsave("Output/Maps/Sites/maptPCBFoxRiverAveV01.png", plot = maptPCBFoxRiver,
       width = 9, height = 4, dpi = 300)

# Hudson River ------------------------------------------------------------
{
  # Select only from Hudson River
  wdc.Hud <- subset(wdc, LocationName == "Hudson River")
  # Increase the longitude range to make the map wider
  lon_range <- 0.5  # Modify this value to control the width
  
  # Create a new bounding box with the adjusted longitude range
  Hud.box <- make_bbox(
    lon = c(min(wdc.Hud$Longitude) - lon_range, max(wdc.Hud$Longitude) + lon_range),
    lat = wdc.Hud$Latitude,
    f = 0.1)
  # Fetch the map using the new bounding box
  Hud.map <- get_stamenmap(bbox = Hud.box, zoom = 8)
  
  # Plot map with sites
  # Prepare data
  # Get tPCB and coordinates
  tPCB.Hud <- data.frame(cbind(wdc.Hud$SiteID, wdc.Hud$Latitude,
                               wdc.Hud$Longitude, wdc.Hud$tPCB))
  # Name the columns
  colnames(tPCB.Hud) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.Hud$Latitude <- as.numeric(tPCB.Hud$Latitude)
  tPCB.Hud$Longitude <- as.numeric(tPCB.Hud$Longitude)
  tPCB.Hud$tPCB <- as.numeric(tPCB.Hud$tPCB)
  # Average tPCB per site
  tPCB.Hud.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                             data = tPCB.Hud, FUN = mean)
}

# (1) Plot map + locations
ggmap(Hud.map) +
  geom_point(data = tPCB.Hud.ave, aes(x = Longitude, y = Latitude),
             shape = 21, color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.Hud.ave, family = 'Times', size = 1.8, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
maptPCBHudsonRiver <- ggmap(Hud.map) +
  geom_point(data = tPCB.Hud.ave, aes(x = Longitude, y = Latitude, size = tPCB),
             alpha = 1, color = "black", shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotate('text', x = -74, y = 43.4,
           label = 'Hudson River (NY)', colour = 'black', size = 3.5,
           fontface = 2) +
  scale_size_area(
    breaks = c(2000, 5000, 10000, 20000, 30000, 40000),
    labels = comma,
    name = expression(bold(Sigma*"PCBs (mean) 2005-2017 (pg/L)")),
    max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = "right",  # Move the legend to the bottom
        legend.justification = c(0.0, 1.05),  # Center the legend horizontally
        legend.title = element_text(margin = margin(b = -1, unit = "pt")))

# Save map in folder
ggsave("Output/Maps/Sites/maptPCBHudsonRiverAveV01.png", plot = maptPCBHudsonRiver,
       width = 6, height = 4, dpi = 300)

# Housatonic River --------------------------------------------------------
{
  # Select only from Housotonic River
  wdc.Hou <- subset(wdc, LocationName == "Housatonic River")
  
  # Increase the longitude range to make the map wider
  lon_range <- 0.1  # Modify this value to control the width
  
  # Create a new bounding box with the adjusted longitude range
  Hou.box <- make_bbox(
    lon = c(min(wdc.Hou$Longitude) - lon_range, max(wdc.Hou$Longitude) + lon_range),
    lat = wdc.Hou$Latitude,
    f = 0.4)
  # Fetch the map using the new bounding box
  Hou.map <- get_stamenmap(bbox = Hou.box, zoom = 8)
  
  # Plot map with sites
  # Prepare data
  # Get tPCB and coordinates
  tPCB.Hou <- data.frame(cbind(wdc.Hou$SiteID, wdc.Hou$Latitude,
                               wdc.Hou$Longitude, wdc.Hou$tPCB))
  # Name the columns
  colnames(tPCB.Hou) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.Hou$Latitude <- as.numeric(tPCB.Hou$Latitude)
  tPCB.Hou$Longitude <- as.numeric(tPCB.Hou$Longitude)
  tPCB.Hou$tPCB <- as.numeric(tPCB.Hou$tPCB)
  # Average tPCB per site
  tPCB.Hou.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                            data = tPCB.Hou, FUN = mean)
}

# (1) Plot map + locations
ggmap(Hou.map) +
  geom_point(data = tPCB.Hou.ave, aes(x = Longitude, y = Latitude),
             shape = 21, color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.Hou.ave, family = 'Times New Roman', size = 3, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# Plot map with sites and tPCB
maptPCBHouRiver <- ggmap(Hou.map) +
  geom_point(data = tPCB.Hou.ave, aes(x = Longitude, y = Latitude, size = tPCB),
             alpha = 1, color = "black", shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotate('text', x = -73.34, y = 42.57,
           label = 'Housatonic River (CT+MA)', colour = 'black', size = 2.8,
           fontface = 2) +
  scale_size_area(
    breaks = c(50000, 100000, 200000, 300000, 400000),
    labels = comma,
    name = expression(bold(Sigma*"PCBs (mean) 1979-2020 (pg/L)")),
    max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(1.8, 0.78),  # Adjust the legend position (x, y)
        legend.box.just = "center",  # Center the legend inside the bounding box
        legend.title = element_text(margin = margin(b = -1, unit = "pt")))

# Save map in folder
ggsave("Output/Maps/Sites/maptPCBHousatonicRiverAveV02.png",
       plot = maptPCBHouRiver, width = 8, height = 4, dpi = 300)

# Kalamazoo River ---------------------------------------------------------
{
  # Select only from Kalamazoo  River
  wdc.Kal <- subset(wdc, LocationName == "Kalamazoo River")
  
  # Create general map
  Kal.box <- make_bbox(lon = wdc.Kal$Longitude, lat = wdc.Kal$Latitude, f = 0.1)
  Kal.map <- get_stamenmap(bbox = Kal.box, zoom = 10)
  
  # Plot map with sites
  # Prepare data
  # Get tPCB and coordinates
  tPCB.Kal <- data.frame(cbind(wdc.Kal$SiteID, wdc.Kal$Latitude,
                               wdc.Kal$Longitude, wdc.Kal$tPCB))
  # Name the columns
  colnames(tPCB.Kal) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.Kal$Latitude <- as.numeric(tPCB.Kal$Latitude)
  tPCB.Kal$Longitude <- as.numeric(tPCB.Kal$Longitude)
  tPCB.Kal$tPCB <- as.numeric(tPCB.Kal$tPCB)
  # Average tPCB per site
  tPCB.Kal.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                         data = tPCB.Kal, FUN = mean)
}

# (1) Plot map + locations
ggmap(Kal.map) +
  geom_point(data = tPCB.Kal.ave, aes(x = Longitude, y = Latitude),
             shape = 21, color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.Kal.ave, family = 'Times New Roman', size = 3, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
maptPCBKalRiver <- ggmap(Kal.map) +
  geom_point(data = tPCB.Kal.ave, aes(x = Longitude, y = Latitude,
                                      size = tPCB), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotate('text', x = -85.6, y = 42.67,
           label = 'Kalamazoo River (MI)', colour = 'black', size = 2.8,
           fontface = 2) +
  scale_size_area(breaks = c(100, 1000, 10000, 50000, 700000), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean) 1994-2010 (pg/L)"))),
                  max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(1.33, 0.75),  # Adjust the legend.position values
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))

# Save map in folder
ggsave("Output/Maps/Sites/maptPCBKalamazooAveV01.png",
       plot = maptPCBKalRiver, width = 12, height = 4, dpi = 300)

# New Bedford -------------------------------------------------------------
{
  # Select only from Hudson River
  wdc.NB <- subset(wdc, LocationName == "New Bedford Harbor")
  
  # Increase the longitude range to make the map wider
  lon_range <- 0.021  # Modify this value to control the width
  
  # Create a new bounding box with the adjusted longitude range
  NB.box <- make_bbox(
    lon = c(min(wdc.NB$Longitude) - lon_range, max(wdc.NB$Longitude) + lon_range),
    lat = wdc.NB$Latitude,
    f = 0.08)
  NB.map <- get_stamenmap(bbox = NB.box, zoom = 12)
  
  # Plot map with sites
  # Prepare data
  # Get tPCB and coordinates
  tPCB.NB <- data.frame(cbind(wdc.NB$SiteID, wdc.NB$Latitude, wdc.NB$Longitude,
                              wdc.NB$tPCB))
  # Name the columns
  colnames(tPCB.NB) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.NB$Latitude <- as.numeric(tPCB.NB$Latitude)
  tPCB.NB$Longitude <- as.numeric(tPCB.NB$Longitude)
  tPCB.NB$tPCB <- as.numeric(tPCB.NB$tPCB)
  # Average tPCB per site
  tPCB.NB.ave <- aggregate(tPCB ~ SiteID+ Latitude + Longitude,
                           data = tPCB.NB, FUN = mean)
}

# (1) Plot map + locations
ggmap(NB.map) +
  geom_point(data = tPCB.NB.ave, aes(x = Longitude, y = Latitude), shape = 21,
             color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                 data = tPCB.NB.ave, family = 'Times', size = 1.8, 
                 box.padding = 0.2, point.padding = 0.3,
                 segment.color = 'grey50')

# (2) Plot map + tPCB
maptPCBNB <- ggmap(NB.map) +
  geom_point(data = tPCB.NB.ave, aes(x = Longitude, y = Latitude,
                                     size = tPCB), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotate('text', x = -70.918, y = 41.699,
           label = 'New Bedford Harbor (MA)', colour = 'black', size = 3.4,
           fontface = 2) +
  scale_size_area(breaks = c(10000, 500000, 1000000, 3000000), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean) 2006-2016 (pg/L)"))),
                  max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(1.65, 0.78),  # Adjust the legend.position values
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))

# Save map in folder
ggsave("Output/Maps/Sites/maptPCBNBHAveV01.png",
       plot = maptPCBNB, width = 9, height = 4, dpi = 300)

# Passaic River -----------------------------------------------------------
{
  # Select only from Hudson River
  wdc.Pas <- subset(wdc, LocationName == "Passaic River")
  
  # Increase the longitude range to make the map wider
  lon_range <- 0.1  # Modify this value to control the width
  
  # Create a new bounding box with the adjusted longitude range
  Pas.box <- make_bbox(
    lat = wdc.Pas$Latitude,
    lon = c(min(wdc.Pas$Longitude) - lon_range, max(wdc.Pas$Longitude) + lon_range),
    f = 0.12)
  Pas.map <- get_stamenmap(bbox = Pas.box, zoom = 10)
  
  # Plot map with sites
  # Prepare data
  # Get tPCB and coordinates
  tPCB.Pas <- data.frame(cbind(wdc.Pas$SiteID, wdc.Pas$Latitude,
                               wdc.Pas$Longitude, wdc.Pas$tPCB))
  # Name the columns
  colnames(tPCB.Pas) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.Pas$Latitude <- as.numeric(tPCB.Pas$Latitude)
  tPCB.Pas$Longitude <- as.numeric(tPCB.Pas$Longitude)
  tPCB.Pas$tPCB <- as.numeric(tPCB.Pas$tPCB)
  # Average tPCB per site
  tPCB.Pas.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                            data = tPCB.Pas, FUN = mean)
}

# (1) Plot map + locations
ggmap(Pas.map) +
  geom_point(data = tPCB.Pas.ave, aes(x = Longitude, y = Latitude), shape = 21,
             color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.Pas.ave, family = 'Times', size = 1.8, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
maptPCBPas <- ggmap(Pas.map) +
  geom_point(data = tPCB.Pas.ave, aes(x = Longitude, y = Latitude,
                                      size = tPCB), alpha = 1,
             color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotate('text', x = -117.43, y = 47.83,
           label = 'Passaic River (NJ)', colour = 'black', size = 3.4,
           fontface = 2) +
  scale_size_area(breaks = c(10, 100, 1000, 10000, 100000), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean) 2014-2016 (pg/L)"))),
                  max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(1.23, 0.75),  # Adjust the legend.position values
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))

# Save map in folder
ggsave("Output/Maps/Sites/maptPCBSpoAveV01.png",
       plot = maptPCBSpo, width = 12, height = 4, dpi = 300)

# Spokane River -----------------------------------------------------------
{
  # Select only from Hudson River
  wdc.Spo <- subset(wdc, LocationName == "Spokane River")
  # Increase the longitude range to make the map wider
  lat_range <- 0.1  # Modify this value to control the width
  
  # Create a new bounding box with the adjusted longitude range
  Spo.box <- make_bbox(
    lon = wdc.Spo$Longitude,
    lat = c(min(wdc.Spo$Latitude) - lat_range, max(wdc.Spo$Latitude) + lat_range),
    f = 0.12)
  Spo.map <- get_stamenmap(bbox = Spo.box, zoom = 10)
  
  # Plot map with sites
  # Prepare data
  # Get tPCB and coordinates
  tPCB.Spo <- data.frame(cbind(wdc.Spo$SiteID, wdc.Spo$Latitude,
                               wdc.Spo$Longitude, wdc.Spo$tPCB))
  # Name the columns
  colnames(tPCB.Spo) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.Spo$Latitude <- as.numeric(tPCB.Spo$Latitude)
  tPCB.Spo$Longitude <- as.numeric(tPCB.Spo$Longitude)
  tPCB.Spo$tPCB <- as.numeric(tPCB.Spo$tPCB)
  # Average tPCB per site
  tPCB.Spo.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                            data = tPCB.Spo, FUN = mean)
}

# (1) Plot map + locations
ggmap(Spo.map) +
  geom_point(data = tPCB.Spo.ave, aes(x = Longitude, y = Latitude), shape = 21,
             color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.Spo.ave, family = 'Times', size = 1.8, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
maptPCBSpo <- ggmap(Spo.map) +
  geom_point(data = tPCB.Spo.ave, aes(x = Longitude, y = Latitude,
                                      size = tPCB), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotate('text', x = -117.43, y = 47.83,
           label = 'Spokane River (WA)', colour = 'black', size = 3.4,
           fontface = 2) +
  scale_size_area(breaks = c(200, 1000, 2000, 4000, 8000), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean) 2014-2016 (pg/L)"))),
                  max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(1.23, 0.75),  # Adjust the legend.position values
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))

# Save map in folder
ggsave("Output/Maps/Sites/maptPCBSpoAveV01.png",
       plot = maptPCBSpo, width = 12, height = 4, dpi = 300)

# Chesapeake Bay ----------------------------------------------------------
{
  # Select only from Chesapeake Bay
  wdc.Che <- subset(wdc, LocationName == "Chesapeake Bay")
  # Increase the longitude range to make the map wider
  lat_range <- 0.05  # Modify this value to control the width
  
  # Create a new bounding box with the adjusted longitude range
  Che.box <- make_bbox(
    lon = wdc.Che$Longitude,
    lat = c(min(wdc.Che$Latitude) - lat_range, max(wdc.Che$Latitude) + lat_range),
    f = 0.15)
  Che.map <- get_stamenmap(bbox = Che.box, zoom = 9)
  
  # Plot map with sites
  # Prepare data
  # Get tPCB and coordinates
  tPCB.Che <- data.frame(cbind(wdc.Che$SiteID, wdc.Che$Latitude,
                               wdc.Che$Longitude, wdc.Che$tPCB))
  # Name the columns
  colnames(tPCB.Che) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.Che$Latitude <- as.numeric(tPCB.Che$Latitude)
  tPCB.Che$Longitude <- as.numeric(tPCB.Che$Longitude)
  tPCB.Che$tPCB <- as.numeric(tPCB.Che$tPCB)
  # Average tPCB per site
  tPCB.Che.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                            data = tPCB.Che, FUN = mean)
}

# (1) Plot map + locations
ggmap(Che.map) +
  geom_point(data = tPCB.Che.ave, aes(x = Longitude, y = Latitude), shape = 21,
             color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.Spo.ave, family = 'Times', size = 1.8, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB 
maptPCBChe <- ggmap(Che.map) +
  geom_point(data = tPCB.Che.ave, aes(x = Longitude, y = Latitude,
                                      size = tPCB), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotate('text', x = -76.15, y = 40.00,
           label = 'Chesapeake Bay (MA)', colour = 'black', size = 2.9,
           fontface = 2) +
  scale_size_area(breaks = c(1000, 5000, 10000, 20000, 30000), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean) 2001-2015 (pg/L)"))),
                  max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(2.05, 0.78),  # Adjust the legend.position values
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))

# Save map in folder
ggsave("Output/Maps/Sites/maptPCBCheAveV01.png",
       plot = maptPCBChe, width = 7, height = 4, dpi = 300)

# Blue River --------------------------------------------------------------
{
  # Select only from Blue River
  wdc.Blu <- subset(wdc, LocationName == "Blue River")
  
  # Create general map
  Blu.box <- make_bbox(lon = wdc.Blu$Longitude, lat = wdc.Blu$Latitude, f = 0.12)
  Blu.map <- get_stamenmap(bbox = Blu.box, zoom = 17)
  
  # Plot map with sites
  # Prepare data
  # Remove samples with total PCBs  = 0
  wdc.Blu.1 <- wdc.Blu[!(wdc.Blu$tPCB) == 0, ]
  # Get tPCB and coordinates
  tPCB.Blu <- data.frame(cbind(wdc.Blu.1$SiteID, wdc.Blu.1$Latitude,
                               wdc.Blu.1$Longitude, wdc.Blu.1$tPCB))
  # Name the columns
  colnames(tPCB.Blu) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.Blu$Latitude <- as.numeric(tPCB.Blu$Latitude)
  tPCB.Blu$Longitude <- as.numeric(tPCB.Blu$Longitude)
  tPCB.Blu$tPCB <- as.numeric(tPCB.Blu$tPCB)
  # Average tPCB per site
  tPCB.Blue.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                             data = tPCB.Blu, FUN = mean)
}

# (1) Plot map + locations
ggmap(Blu.map) +
  geom_point(data = tPCB.Blue.ave, aes(x = Longitude, y = Latitude), shape = 21,
             color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.Blue.ave, family = 'Times', size = 1.8, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(Blu.map) +
  geom_point(data = tPCB.Blue.ave, aes(x = Longitude, y = Latitude,
                                       size = tPCB), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotate('text', x = -94.572, y = 38.9685,
           label = 'Blue River (MO)', colour = 'black', size = 3.4,
           fontface = 2) +
  scale_size_area(breaks = c(1000, 10000, 80000, 120000, 160000), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean) 2004-2019 (pg/L)"))),
                  max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(1.45, 0.75),  # Adjust the legend.position values
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))


# Extra -------------------------------------------------------------------
# Prepare congener data for plotting
# Get coordinates per site
LW <- subset(w.WI, w.WI$SiteSampled == 'LakeWinnebago')
LW <- data.frame(c(LW[1,6], LW[1,7]))
OU1 <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit1')
OU1 <- data.frame(c(OU1[1,6], OU1[1,7]))
OU2A <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit2A')
OU2A <- data.frame(c(OU2A[1,6], OU2A[1,7]))
OU2B <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit2B')
OU2B <- data.frame(c(OU2B[1,6], OU2B[1,7]))
OU2C <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit2C')
OU2C <- data.frame(c(OU2C[1,6], OU2C[1,7]))
OU3 <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit3')
OU3 <- data.frame(c(OU3[1,6], OU3[1,7]))
wi.coord <- rbind(LW, OU1, OU2A, OU2B, OU2C, OU3)

# Total PCBs
# # remove samples (rows) with total PCBs  = 0
w.WI.t <- w.WI[!(rowSums(w.WI[,
                           c(12:115)],
                         na.rm = TRUE)==0),] # sum of PCB1 to PCB209
site.sampled <- w.WI.t$SiteSampled
w.WI.t <- subset(w.WI.t, select = -c(ID:AroclorCongener))
w.WI.t <- subset(w.WI.t, select = -c(AroclorA1016:AroclorA1260))
# Get mean congener per site, excluding zeros
tPCB <- rowSums(w.WI.t, na.rm = TRUE)
tPCB <- data.frame(cbind(site.sampled, tPCB))
tPCB$tPCB <- as.numeric(as.character(tPCB$tPCB))
tPCB.mean <- aggregate(tPCB ~ site.sampled, data = tPCB, mean)
# add coordinates
tPCB.mean <- data.frame(c(tPCB.mean, wi.coord))

# (3) Plot map + tPCB
ggmap(wi.map) +
  geom_point(data = tPCB.mean, aes(x = Long, y = Lat,
                              size = tPCB), alpha = 0.5) +
  scale_size_area(breaks = c(250, 500, 750, 1000, 1500),
                  labels = c(250, 500, 750, 1000, 1500),
                  name = "PCBs ng/L") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(title = "Fox River PCBs water concentration (pg/L) 2010-2018")

# Congener maps
# Select congener and remove samples with = 0 and NA for selected congener
w.WI.2 <- subset(w.WI, w.WI$PCB1 != 0 & w.WI$PCB1 != "NA")
# Get mean congener per site, excluding zeros
PCB1 <- aggregate(PCB1 ~ SiteSampled, data = w.WI.2, mean)
PCB1 <- data.frame(c(PCB1, wi.coord))

# (4) Plot map + congener
ggmap(wi.map) +
  geom_point(data = PCB1, aes(x = Long, y = Lat,
                              size = PCB1), alpha = 0.5) +
  scale_size_area(breaks = c(0.1, 1, 2, 4, 6),
                 labels = c(0.1, 1, 2, 4, 6),
                 name = "PCB 1 ng/L") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(title = "Fox River PCB 1 water concentration (pg/L) 2010-2018")
  #geom_label_repel(data = PCB1, aes(x = Long, y = Lat, label = SiteSampled),
  #                 fill = "white", box.padding = unit(0.3, "lines"),
  #                 label.padding = unit(0.15, "lines"),
  #                 segment.color = "black", segment.size = 1)
                   
                   
