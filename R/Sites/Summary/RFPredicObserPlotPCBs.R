
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
# Data only with R2 > 0. See file AllRFtPCB.csv
{
  # Anacostia River
  anr <- read.csv("Output/Data/Sites/csv/AnacostiaRiver/AnacostiaRiverRFObsPredtPCB.csv")
  # Bannister Federal Complex
  bfc <- read.csv("Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexRFObsPredtPCB.csv")
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFObsPredtPCB.csv")
  # DEQ MI
  dmi <- read.csv("Output/Data/Sites/csv/DEQMichigan/DEQMIRFObsPredtPCB.csv")
  # Fox River data
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverRFObsPredtPCB.csv")
  # Housatonic River
  hou <- read.csv("Output/Data/Sites/csv/HousatonicRiver/HousatonicRiverRFObsPredtPCB.csv")
  # Hudson River
  hud <- read.csv("Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFObsPredtPCB.csv")
  # Kalamazoo River
  kal <- read.csv("Output/Data/Sites/csv/KalamazooRiver/KalamazooRiverRFObsPredtPCB.csv")
  # New Bedford Harbor data
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHRFObsPredtPCB.csv")
  # Passaic River
  pas <- read.csv("Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFObsPredtPCB.csv")
  # Portland Harbord data
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFObsPredtPCB.csv")
  # Spokane River data
  spo <- read.csv("Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFObsPredtPCB.csv")
  # Combine the data frames
  combined_data_tPCB <- rbind(anr, bfc, che, dmi, fox, hou, hud, kal, nbh, pas,
                              por, spo)
}

# Create a custom color palette with 9 different colors
n_colors <- 12
custom_palette <- colorRampPalette(brewer.pal(9, "Set1"))(n_colors)

# Add a new column with the number of rows for each Location
combined_data_tPCB$Location <- paste0(combined_data_tPCB$Location,
                                 " (n=", ave(combined_data_tPCB$Location,
                                             combined_data_tPCB$Location,
                                             FUN = length), ")")                                                                         

# Plot prediction vs. observations, 1:1 line
CombinePredObsPlot <- ggplot(combined_data_tPCB,
                             aes(x = 10^(Actual), y = 10^(Predicted),
                                 fill = Location)) +  # Use the new column for legend labels
  geom_point(shape = 21, size = 3, alpha = 0.5) +
  scale_y_log10(limits = c(1, 10^8), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^8), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values = custom_palette) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) +
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  theme_bw() +
  theme(aspect.ratio = 15/15,
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.position = c(0.25, 0.8),  # Adjust the position (values are between 0 and 1)
        legend.background = element_rect(fill = "transparent")) +  # Make the legend background transparent
  annotation_logticks(sides = "bl") +
  guides(fill = guide_legend(title = NULL))  # Remove the legend title

# Print plot
print(CombinePredObsPlot)

# Save plot
ggsave("Output/Plots/Sites/ObsPred/Summary/CombineRFObsPredtPCB.png",
       plot = CombinePredObsPlot, width = 18, height = 8, dpi = 2000)

# Read generated data for PCBs --------------------------------------------
# Data only with R2 > 0
{
  # Bannister Federal Complex
  bfc <- read.csv("Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexRFObsPredPCB.csv")
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFObsPredPCB.csv")
  # DEQ MI
  dmi <- read.csv("Output/Data/Sites/csv/DEQMichigan/DEQMIRFObsPredPCB.csv")
  # Fox River data
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverRFObsPredPCB.csv")
  # Great Lakes (Lake Michigan)
  grl <- read.csv("Output/Data/Sites/csv/GreatLakes/GreatLakesRFObsPredPCB.csv")
  # Hudson River
  hud <- read.csv("Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFObsPredPCB.csv")
  # Lake Washington
  lwa <- read.csv("Output/Data/Sites/csv/LakeWashington/LakeWashingtonRFObsPredPCB.csv")
  # New Bedford Harbor data
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHRFObsPredPCB.csv")
  # Passaic River
  pas <- read.csv("Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFObsPredPCB.csv")
  # Portland Harbord data
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFObsPredPCB.csv")
  # Spokane River data
  spo <- read.csv("Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFObsPredPCB.csv")
  # Combine the data frames
  combined_data_PCB <- rbind(bfc, che, dmi, fox, grl, hud, lwa, nbh, pas, por, spo)
}

# Create a custom color palette with 9 different colors
n_colors <- 11
custom_palette <- colorRampPalette(brewer.pal(9, "Set1"))(n_colors)

# Add a new column with the number of rows for each Location
combined_data_PCB$Location <- paste0(combined_data_PCB$Location,
                                 " (n=", ave(combined_data_PCB$Location,
                                             combined_data_PCB$Location,
                                             FUN = length), ")")                                                                         

# Plot prediction vs. observations, 1:1 line
CombinePredObsPlot <- ggplot(combined_data_PCB,
                             aes(x = 10^(Actual), y = 10^(Predicted),
                                 fill = Location)) +  # Use the new column for legend labels
  geom_point(shape = 21, size = 3, alpha = 0.5) +
  scale_y_log10(limits = c(0.001, 10^7), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.001, 10^7), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values = custom_palette) +
  xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
  ylab(expression(bold("Predicted concentration PCBi (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) +
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  theme_bw() +
  theme(aspect.ratio = 15/15,
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.position = c(0.25, 0.82),  # Adjust the position (values are between 0 and 1)
        legend.background = element_rect(fill = "transparent")) +  # Make the legend background transparent
  annotation_logticks(sides = "bl") +
  guides(fill = guide_legend(title = NULL))  # Remove the legend title

# Print plot
print(CombinePredObsPlot)

# Save plot
ggsave("Output/Plots/Sites/ObsPred/Summary/CombineRFObsPredPCB.png",
       plot = CombinePredObsPlot, width = 18, height = 8, dpi = 2000)

