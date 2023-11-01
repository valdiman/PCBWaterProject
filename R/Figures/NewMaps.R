# Install packages
install.packages("leaflet")
install.packages("dplyr")

# Load libraries
library(leaflet)
library(dplyr)

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

wdc.PO <- subset(wdc, LocationName == "Portland Harbor")
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

# Create a color palette based on tPCB values
color_palette <- colorNumeric(
  palette = "YlOrRd",  # Choose a color palette (you can change this)
  domain = tPCB.PO.ave$tPCB
)

# Create a map
m <- leaflet(data = tPCB.PO.ave) %>%
  addTiles()  # Use OpenStreetMap tiles

# Add markers for each data point with colored circles based on tPCB values
m <- m %>%
  addCircleMarkers(
    lng = ~Longitude,
    lat = ~Latitude,
    radius = 8,
    color = "black",  # Border color
    fillColor = ~color_palette(tPCB),  # Fill color based on tPCB values
    fillOpacity = 0.8,
    weight = 1,
    popup = ~SiteID
  )

# Define the range of tPCB values for the legend
tPCB_range <- range(tPCB.PO.ave$tPCB)

m <- m %>%
  addLegend(
    "bottomright",
    title = "ΣPCB (pg/L)",
    values = tPCB_range,
    pal = color_palette,
    opacity = 0.8,
    labels = format(tPCB_range, scientific = FALSE)
  )

# Display the map
m
