# Code to get coordinates per site to plot in Google Earth

# Install packages
install.packages("tidyverse")
install.packages("dplyr")
install.packages("sf")
install.packages("sfheaders")

# Load libraries
{
  library(dplyr) # performs %>%
  library(stringr) # str_detect
  library(sf) # Create file to be used in Google Earth
  library(sfheaders) # Create file to be used in Google Earth
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# Total PCB description ---------------------------------------------------
# Data preparation
tpcb <- cbind(wdc$SiteID, wdc$Latitude, wdc$Longitude, wdc$tPCB)
tpcb <- as.data.frame(tpcb)
# Transform characters to numeric values
tpcb[, 2:4] <- lapply(tpcb[, 2:4], as.numeric)
# Add column names
colnames(tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Global/PCBSampleLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select 21 Michigan data ---------------------------------------------------
mic <- wdc[str_detect(wdc$LocationName, '21Mich'),]

# Data preparation --------------------------------------------------------
# Create data frame
mic.tpcb <- cbind(mic$SiteID, mic$Latitude, mic$Longitude, mic$tPCB)
mic.tpcb <- as.data.frame(mic.tpcb)
# Transform characters to numeric values
mic.tpcb[, 2:4] <- lapply(mic.tpcb[, 2:4], as.numeric)
# Add column names
colnames(mic.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = mic.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/21MichLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Here!

# Select Anacostia River data ---------------------------------------------------
anr <- wdc[str_detect(wdc$LocationName, 'Anacostia River'),]

# Data preparation --------------------------------------------------------
{
  # Create individual code for each site sampled
  site.numb <- anr$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  anr.tpcb <- cbind(factor(anr$SiteID), anr$Latitude, anr$Longitude,
                    as.matrix(anr$tPCB))
  # Add column names
  colnames(anr.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = anr.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/AnacostiaRiverLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select Bannister Federal Complex data ---------------------------------------------------
bfc <- wdc[str_detect(wdc$LocationName, 'Bannister Fed Complex'),]

# Data preparation --------------------------------------------------------
{
  site.numb <- bfc$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  bfc.tpcb <- cbind(factor(bfc$SiteID), bfc$Latitude, bfc$Longitude,
                    as.matrix(bfc$tPCB))
  # Add column names
  colnames(bfc.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = bfc.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/BannisterFedComplexLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select Chesapeake Bay & Delaware Canal data ---------------------------------------------------
che <- wdc[str_detect(wdc$LocationName, 'Chesapeake Bay'),]

# Data preparation --------------------------------------------------------
{
  # Create individual code for each site sampled
  site.numb <- che$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  che.tpcb <- cbind(factor(che$SiteID), che$Latitude, che$Longitude,
                    as.matrix(che$tPCB))
  # Add column names
  colnames(che.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = che.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/ChesapeakeBayLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select Fox River data ---------------------------------------------------
fox <- wdc[str_detect(wdc$LocationName, 'Fox River'),]
# Lake Winnebago is a background site.
# Data preparation --------------------------------------------------------
{
  # Create individual code for each site sampled
  site.numb <- fox$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  fox.tpcb <- cbind(factor(fox$SiteID), fox$Latitude, fox$Longitude,
                    as.matrix(fox$tPCB))
  # Add column names
  colnames(fox.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = fox.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/FoxRiverLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select LMMB and Great Lakes data ---------------------------------------------------
grl <- wdc[str_detect(wdc$LocationName, 'Lake Michigan Mass Balance|Great Lakes'), ]

# (1) Just get lake data, remove data from tributaries
grl <- grl[!grepl("^Tributary", grl$SiteName), ]

# Data preparation --------------------------------------------------------
{
  # Create individual code for each site sampled
  site.numb <- grl$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  grl.tpcb <- cbind(factor(grl$SiteID), grl$Latitude, grl$Longitude,
                    as.matrix(grl$tPCB))
  # Add column names
  colnames(grl.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = grl.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/GreatLakesLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select Housatonic River data ---------------------------------------------------
hou <- wdc[str_detect(wdc$LocationName, 'Housatonic River'),]

# Data preparation --------------------------------------------------------
{
  # Create individual code for each site sampled
  site.numb <- hou$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  hou.tpcb <- cbind(factor(hou$SiteID), hou$Latitude, hou$Longitude,
                    as.matrix(hou$tPCB))
  # Add column names
  colnames(hou.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = hou.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/HousatonicRiverLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select Hudson River data ---------------------------------------------------
hud <- wdc[str_detect(wdc$LocationName, 'Hudson River'),]
# PCBs were discharged to the river from the General Electric
# (GE) manufacturing plants in Hudson Falls and Fort Edward, NY
# Dredging from 2009 to 2015
# https://www.epa.gov/system/files/documents/2021-08/hudson_summer2021_floodplainrifs_factsheet_final.pdf

# Data preparation --------------------------------------------------------
{
  # Create individual code for each site sampled
  site.numb <- hud$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  hud.tpcb <- cbind(factor(hud$SiteID), hud$Latitude, hud$Longitude,
                    as.matrix(hud$tPCB))
  # Add column names
  colnames(hud.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = hud.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/HudsonRiverLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select Kalamazoo data ---------------------------------------------------
kal <- wdc[str_detect(wdc$LocationName, 'Kalamazoo River'),]
# Superfund site from Morrow Dam (Kalamazoo River) to Lake Michigan
# and 30 miles of Portage Creek (south), Cork St and Portage Creek Cork St sites
# Dredging occurred at Plainwell Dam site.

# Data preparation --------------------------------------------------------
{
  # Create individual code for each site sampled
  site.numb <- kal$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  kal.tpcb <- cbind(factor(kal$SiteID), kal$Latitude, kal$Longitude,
                    as.matrix(kal$tPCB))
  # Add column names
  colnames(kal.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = kal.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/KalamazooRiverLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select Lake Washington data ---------------------------------------------------
lwa <- wdc[str_detect(wdc$LocationName, 'Lake Washington'),]

# Data preparation --------------------------------------------------------
{
  # Create individual code for each site sampled
  site.numb <- lwa$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  lwa.tpcb <- cbind(factor(lwa$SiteID), lwa$Latitude, lwa$Longitude,
                    as.matrix(lwa$tPCB))
  # Add column names
  colnames(lwa.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = lwa.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/LakeWashingtonLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select nbh River data ---------------------------------------------------
nbh <- wdc[str_detect(wdc$LocationName, 'New Bedford'),]

# Data preparation --------------------------------------------------------
{
  # Create individual code for each site sampled
  site.numb <- nbh$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  nbh.tpcb <- cbind(factor(nbh$SiteID), nbh$Latitude, nbh$Longitude,
                    as.matrix(nbh$tPCB))
  # Add column names
  colnames(nbh.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = nbh.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/NewBedfordHarborLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select Passaic River data ---------------------------------------------------
pas <- wdc[str_detect(wdc$LocationName, 'Passaic River'),]

# Data preparation --------------------------------------------------------
{
  # Create individual code for each site sampled
  site.numb <- pas$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  pas.tpcb <- cbind(factor(pas$SiteID), pas$Latitude, pas$Longitude,
                    as.matrix(pas$tPCB))
  # Add column names
  colnames(pas.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = pas.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/PassaicRiverLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select Portland Harbor data ---------------------------------------------------
por <- wdc[str_detect(wdc$LocationName, 'Portland Harbor'),]

# Data preparation --------------------------------------------------------
{
  # Create individual code for each site sampled
  site.numb <- por$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  por.tpcb <- cbind(factor(por$SiteID), por$Latitude, por$Longitude,
                    as.matrix(por$tPCB))
  # Add column names
  colnames(por.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = por.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/PortlandHarborLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select Richardson Hill Road Landfill data ---------------------------------------------------
rhr <- wdc[str_detect(wdc$LocationName, 'Richardson Hill Road Landfill'),]

# Data preparation --------------------------------------------------------
{
  # Create individual code for each site sampled
  site.numb <- rhr$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  rhr.tpcb <- cbind(factor(rhr$SiteID), rhr$Latitude, rhr$Longitude,
                    as.matrix(rhr$tPCB))
  # Add column names
  colnames(rhr.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = rhr.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/RichardsonHillLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select Spokane River data ---------------------------------------------------
spo <- wdc[str_detect(wdc$LocationName, 'Spokane River'),]

# Data preparation --------------------------------------------------------
{
  # Create individual code for each site sampled
  site.numb <- spo$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  spo.tpcb <- cbind(factor(spo$SiteID), spo$Latitude, spo$Longitude,
                    as.matrix(spo$tPCB))
  # Add column names
  colnames(spo.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = spo.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/SpokaneRiverLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Select tributary data from LMMB  data ---------------------------------------------------
# Just tributaries data
glt <- wdc[grepl("^Tributary", wdc$SiteName), ]

# Data preparation --------------------------------------------------------
{
  # Create individual code for each site sampled
  site.numb <- glt$SiteID %>% as.factor() %>% as.numeric
  # Create data frame
  glt.tpcb <- cbind(factor(glt$SiteID), glt$Latitude, glt$Longitude,
                    as.matrix(glt$tPCB))
  # Add column names
  colnames(glt.tpcb) <- c("SiteID", "Latitude", "Longitude", "tPCB")
}

# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = glt.tpcb, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Sites/GoogleEarth/TributariesLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)
