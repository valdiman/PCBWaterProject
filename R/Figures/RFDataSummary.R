
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

# Read generated data for total PCB
{
  # 21 Mich
  mic <- read.csv("Output/Data/Sites/csv/21Mich/21MichRFPerformancetPCB.csv")
  mic <- mic[, 3]
  # Anacostia River
  anr <- read.csv("Output/Data/Sites/csv/AnacostiaRiver/AnacostiaRiverRFPerformancetPCB.csv")
  anr <- anr[, 3]
  # Bannister Federal Complex
  bfc <- read.csv("Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexRFPerformancetPCB.csv")
  bfc <- bfc[, 3]
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFPerformancetPCB.csv")
  che <- che[, 3]
  # Fox River data
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverRFPerformancetPCB.csv")
  fox <- fox[, 3]
  # Housatonic River
  hou <- read.csv("Output/Data/Sites/csv/HousatonicRiver/HousatonicRiverRFPerformancetPCB.csv")
  hou <- hou[, 3]
  # Hudson River
  hud <- read.csv("Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFPerformancetPCB.csv")
  hud <- hud[, 3]
  # Kalamazoo River
  kal <- read.csv("Output/Data/Sites/csv/KalamazooRiver/KalamazooRiverRFPerformancetPCB.csv")
  kal <- kal[, 3]
  # New Bedford Harbor data
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHRFPerformancetPCB.csv")
  nbh <- nbh[, 3]
  # Passaic River
  pas <- read.csv("Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFPerformancetPCB.csv")
  pas <- pas[, 3]
  # Portland Harbor data
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFPerformancetPCBV02.csv")
  por <- por[, 3]
  # Spokane River data
  spo <- read.csv("Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFPerformancetPCBV02.csv")
  spo <- spo[, 3]
  # Combine the data frames
  combined_data <- rbind(anr, bfc,  che, mic, fox, hou, hud, kal, nbh,
                         pas, por, spo)
  colnames(combined_data) <- c("RMSE", "R2", "Factor2")
}

# R2R
min(combined_data)*100
max(combined_data)*100
mean(combined_data)*100
sd(combined_data)*100

# t05
list_loc <- list(mic, bfc, che, fox, grl, glt, hud, lwa, nbh, pas, por, spo)

t05 <- list2env(setNames(lapply(list_loc, function(df) {
  df %>%
    filter(time.pv < 0.05) %>%
    select(Congeners, t05, t05.error)
}), c("mic_t05", "bfc_t05", "che_t05", "fox_t05", "grl_t05", "glt_t05",
      "hud_t05", "lwa_t05", "nbh_t05", "pas_t05", "por_t05", "spo_t05")))

# Combine the data frames into a single data frame
combined_t05 <- bind_rows(t05$bfc_t05, t05$che_t05)

# Combine all data frames into a single data frame
combined_t05 <- bind_rows(lapply(t05, function(df) as.data.frame(df)))


