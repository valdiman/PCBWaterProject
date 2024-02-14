
# Install packages
install.packages("ggplot2")
install.packages("scales")
install.packages("RColorBrewer")

# Load libraries
{
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
}

# Read generated data for total PCB ---------------------------------------
{
  # 21 Mich
  mic <- read.csv("Output/Data/Sites/csv/21Mich/21MichRFPerformancetPCB.csv")
  mic <- mic[, 2]
  # Anacostia River
  anr <- read.csv("Output/Data/Sites/csv/AnacostiaRiver/AnacostiaRiverRFPerformancetPCB.csv")
  anr <- anr[, 2]
  # Bannister Federal Complex
  bfc <- read.csv("Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexRFPerformancetPCB.csv")
  bfc <- bfc[, 2]
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFPerformancetPCB.csv")
  che <- che[, 2]
  # Fox River
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverRFPerformancetPCB.csv")
  fox <- fox[, 2]
  # Great Lakes (Lake Michigan)
  lmi <- read.csv("Output/Data/Sites/csv/GreatLakes/GreatLakesRFPerformancetPCB.csv")
  lmi <- lmi[, 2]
  # Housatonic River
  hou <- read.csv("Output/Data/Sites/csv/HousatonicRiver/HousatonicRiverRFPerformancetPCB.csv")
  hou <- hou[, 2]
  # Hudson River
  hud <- read.csv("Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFPerformancetPCB.csv")
  hud <- hud[, 2]
  # Kalamazoo River
  kal <- read.csv("Output/Data/Sites/csv/KalamazooRiver/KalamazooRiverRFPerformancetPCB.csv")
  kal <- kal[, 2]
  # Lake Washington
  lwa <- read.csv("Output/Data/Sites/csv/LakeWashington/LakeWashingtonRFPerformancetPCB.csv")
  lwa <- lwa[, 2]
  # New Bedford Harbor
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHRFPerformancetPCB.csv")
  nbh <- nbh[, 2]
  # Passaic River
  pas <- read.csv("Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFPerformancetPCB.csv")
  pas <- pas[, 2]
  # Portland Harbor
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFPerformancetPCB.csv")
  por <- por[, 2]
  # Richardson Hill Road Landfill
  rhl <- read.csv("Output/Data/Sites/csv/Richardson/RichardsonRFPerformancetPCB.csv")
  rhl <- rhl[, 2]
  # Spokane River
  spo <- read.csv("Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFPerformancetPCB.csv")
  spo <- spo[, 2]
  # Combine the data frames
  combined_data <- rbind(anr, bfc,  che, mic, fox, lmi, hou, hud, kal, nbh,
                         lwa, pas, por, rhl, spo)
  colnames(combined_data) <- c("RMSE", "R2", "Factor2")
}

# see data
print(combined_data)

# Export results
write.csv(combined_data, file = "Output/Data/Sites/csv/AllRFPerformancetPCBV0.csv")

# Select locations with R2 > 0.4
combined_data.2 <- combined_data[combined_data[, 2] > 0.4, ]

# see data
print(combined_data.2)

# Export results
write.csv(combined_data.2, file = "Output/Data/Sites/csv/AllRFPerformancetPCB.csv")

# Read generated data for individual PCB ----------------------------------
{
  # 21 Mich
  mic <- read.csv("Output/Data/Sites/csv/21Mich/21MichRFPerformancePCB.csv")
  # Bannister Federal Complex
  bfc <- read.csv("Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexRFPerformancePCB.csv")
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFPerformancePCB.csv")
  # Fox River data
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverRFPerformancePCB.csv")
  # Great Lakes (Lake Michigan)
  grl <- read.csv("Output/Data/Sites/csv/GreatLakes/GreatLakesRFPerformancePCB.csv")
  # Hudson River
  hud <- read.csv("Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFPerformancePCB.csv")
  # Lake Washington
  lwa <- read.csv("Output/Data/Sites/csv/LakeWashington/LakeWashingtonRFPerformancePCB.csv")
  # New Bedford Harbor data
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHRFPerformancePCB.csv")
  # Passaic River
  pas <- read.csv("Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFPerformancePCB.csv")
  # Portland Harbor
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFPerformancePCB.csv")
  # Spokane River
  spo <- read.csv("Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFPerformancePCB.csv")
  
  # Combine the data frames
  combined_data <- rbind(mic, bfc, che, fox, grl, hud, lwa, nbh, pas, por, spo)
}

# Export results
write.csv(combined_data, file = "Output/Data/Sites/csv/AllRFPerformancePCBi.csv")



