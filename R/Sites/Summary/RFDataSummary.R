# Summarize performance of random forest model analysis,
# for both tPCB and individual PCBs. RSME, R2, and factor or 2.
# For individual PCBs, only model that yielded R2 > 0 are included.

# Install packages
install.packages('dplyr')

# Load libraries
library(dplyr)

# Read generated data for total PCB ---------------------------------------
{
  # DEQ MI
  dmi <- read.csv("Output/Data/Sites/csv/DEQMichigan/DEQMIRFtPCB.csv")
  dmi <- dmi[, 2]
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
  combined_tPCB <- rbind(anr, bfc,  che, dmi, fox, lmi, hou, hud, kal, nbh,
                         lwa, pas, por, rhl, spo)
  colnames(combined_tPCB) <- c("RMSE", "R2", "Factor2")
}

# see data
print(combined_tPCB)

# Export results
write.csv(combined_tPCB, file = "Output/Data/Sites/csv/Summary/AllRFtPCB.csv")

# Read generated data for individual PCB ----------------------------------
{
  # DEQ MI
  dmi <- read.csv("Output/Data/Sites/csv/DEQMichigan/DEQMIRFPCB.csv")
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
  combined_PCB <- rbind(bfc, che, dmi, fox, grl, hud, lwa, nbh, pas, por, spo)
}

# see data
print(combined_PCB)

# Export results
write.csv(combined_PCB, file = "Output/Data/Sites/csv/Summary/AllRFPCB.csv")

# Summary of individual PCBs ----------------------------------------------
# Calculate the number of congeners per Location
congeners_per_location <- combined_PCB %>%
  group_by(Location) %>%
  summarize(NumberOfCongeners = n())

# View the results
print(congeners_per_location)

# Calculate the mean and standard deviation for RMSE, R_squared, and Factor2_Percentage per Location
stats_per_location <- combined_PCB %>%
  group_by(Location) %>%
  summarize(
    Mean_RMSE = mean(RMSE, na.rm = TRUE),
    SD_RMSE = sd(RMSE, na.rm = TRUE),
    Mean_R_squared = mean(R_squared, na.rm = TRUE),
    SD_R_squared = sd(R_squared, na.rm = TRUE),
    Mean_Factor2_Percentage = mean(Factor2_Percentage, na.rm = TRUE),
    SD_Factor2_Percentage = sd(Factor2_Percentage, na.rm = TRUE)
  )

# View the results
print(stats_per_location)

