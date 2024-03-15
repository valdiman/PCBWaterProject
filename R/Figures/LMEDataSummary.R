
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
# Data only with R2 > 0 and normality (p-value) > 0.05
{
  # Anacostia River
  anr <- read.csv("Output/Data/Sites/csv/AnacostiaRiver/AnacostiaRiverLmetPCB.csv")
  # Bannister Federal Complex
  bfc <- read.csv("Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexLmetPCB.csv")
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeLmetPCB.csv")
  # Fox River
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverLmetPCB.csv")
  # Kalamazoo River
  kal <- read.csv("Output/Data/Sites/csv/KalamazooRiver/KalamazooRiverLmetPCB.csv")
  # Lake Washington
  lwa <- read.csv("Output/Data/Sites/csv/LakeWashington/LakeWashingtonLmetPCB.csv")
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

# Read generated data from LME Model
{
  # 21 Mich
  mic <- read.csv("Output/Data/Sites/csv/21Mich/21MichLmePCB.csv")
  mic.r2r <- mic$R2R
  # Bannister Federal Complex
  bfc <- read.csv("Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexLmePCB.csv")
  bfc.r2r <- bfc$R2R
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeLmePCB.csv")
  che.r2r <- che$R2R
  # Fox River data
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverLmePCB.csv")
  fox.r2r <- fox$R2R
  # Great Lakes
  grl <- read.csv("Output/Data/Sites/csv/GreatLakes/GreatLakesLmePCB.csv")
  grl.r2r <- grl$R2R
  # Hudson River
  hud <- read.csv("Output/Data/Sites/csv/HudsonRiver/HudsonRiverLmePCB.csv")
  hud.r2r <- hud$R2R
  # Lake Washington
  lwa <- read.csv("Output/Data/Sites/csv/LakeWashington/LakeWashingtonLmePCB.csv")
  lwa.r2r <- lwa$R2R
  # Tributaries to Lake Michigan
  glt <- read.csv("Output/Data/Sites/csv/GreatLakes/TributariesLmePCB.csv")
  glt.r2r <- glt$R2R
  # New Bedford Harbor data
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHLmePCB.csv")
  nbh.r2r <- nbh$R2R
  # Passaic River
  pas <- read.csv("Output/Data/Sites/csv/PassaicRiver/PassaicLmePCB.csv")
  pas.r2r <- pas$R2R
  # Portland Harbord data
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/PortlandHarborLmePCB.csv")
  por.r2r <- por$R2R
  # Spokane River data
  spo <- read.csv("Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverLmePCB.csv")
  spo.r2r <- spo$R2R
  # Combine the data frames
  combined_data <- c(mic.r2r, bfc.r2r, che.r2r, fox.r2r, grl.r2r, glt.r2r,
                     hud.r2r, lwa.r2r, nbh.r2r, pas.r2r, por.r2r, spo.r2r)
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


