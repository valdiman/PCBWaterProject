
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
  mic <- read.csv("Output/Data/Sites/csv/21Mich/21MichRFtPCB.csv")
  mic <- mic[, 2]
  # Anacostia River
  anr <- read.csv("Output/Data/Sites/csv/AnacostiaRiver/AnacostiaRiverRFtPCB.csv")
  anr <- anr[, 2]
  # Bannister Federal Complex
  bfc <- read.csv("Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexRFtPCB.csv")
  bfc <- bfc[, 2]
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFtPCB.csv")
  che <- che[, 2]
  # Fox River
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverRFtPCB.csv")
  fox <- fox[, 2]
  # Great Lakes (Lake Michigan)
  lmi <- read.csv("Output/Data/Sites/csv/GreatLakes/GreatLakesRFtPCB.csv")
  lmi <- lmi[, 2]
  # Housatonic River
  hou <- read.csv("Output/Data/Sites/csv/HousatonicRiver/HousatonicRiverRFtPCB.csv")
  hou <- hou[, 2]
  # Hudson River
  hud <- read.csv("Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFtPCB.csv")
  hud <- hud[, 2]
  # Kalamazoo River
  kal <- read.csv("Output/Data/Sites/csv/KalamazooRiver/KalamazooRiverRFtPCB.csv")
  kal <- kal[, 2]
  # Lake Washington
  lwa <- read.csv("Output/Data/Sites/csv/LakeWashington/LakeWashingtonRFtPCB.csv")
  lwa <- lwa[, 2]
  # New Bedford Harbor
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHRFtPCB.csv")
  nbh <- nbh[, 2]
  # Passaic River
  pas <- read.csv("Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFtPCB.csv")
  pas <- pas[, 2]
  # Portland Harbor
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFtPCB.csv")
  por <- por[, 2]
  # Richardson Hill Road Landfill
  rhl <- read.csv("Output/Data/Sites/csv/Richardson/RichardsonRFtPCB.csv")
  rhl <- rhl[, 2]
  # Spokane River
  spo <- read.csv("Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFtPCB.csv")
  spo <- spo[, 2]
  # Combine the data frames
  combined_data <- rbind(anr, bfc,  che, mic, fox, lmi, hou, hud, kal, nbh,
                         lwa, pas, por, rhl, spo)
  colnames(combined_data) <- c("RMSE", "R2", "Factor2")
}

# see data
print(combined_data)

# Export results
write.csv(combined_data, file = "Output/Data/Sites/csv/AllRFtPCB.csv")

# Read generated data for individual PCB ----------------------------------
{
  # 21 Mich
  mic <- read.csv("Output/Data/Sites/csv/21Mich/21MichRFPCB.csv")
  # Bannister Federal Complex
  bfc <- read.csv("Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexRFPCB.csv")
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFPCB.csv")
  # Fox River data
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverRFPCB.csv")
  # Great Lakes (Lake Michigan)
  grl <- read.csv("Output/Data/Sites/csv/GreatLakes/GreatLakesRFPCB.csv")
  # Hudson River
  hud <- read.csv("Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFPCB.csv")
  # Lake Washington
  lwa <- read.csv("Output/Data/Sites/csv/LakeWashington/LakeWashingtonRFPCB.csv")
  # New Bedford Harbor data
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHRFPCB.csv")
  # Passaic River
  pas <- read.csv("Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFPCB.csv")
  # Portland Harbor
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFPCB.csv")
  # Spokane River
  spo <- read.csv("Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFPCB.csv")
  
  # Combine the data frames
  combined_data <- rbind(mic, bfc, che, fox, grl, hud, lwa, nbh, pas, por, spo)
}

# Export results
write.csv(combined_data, file = "Output/Data/Sites/csv/AllRFPCB.csv")



