
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

# Read generated data
{
  # 21 Mich
  mic <- read.csv("Output/Data/Sites/csv/21Mich/21MichObsPredPCB.csv")
  mic <- mic[, -c(1:2)]
  # Bannister Federal Complex
  bfc <- read.csv("Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexObsPredPCB.csv")
  bfc <- bfc[, -c(1:2)]
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeObsPredPCB.csv")
  che <- che[, -c(1:2)]
  # Fox River data
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverObsPredPCB.csv")
  fox <- fox[, -c(1:2)]
  # Great Lakes
  grl <- read.csv("Output/Data/Sites/csv/GreatLakes/GreatLakesObsPredPCB.csv")
  grl <- grl[, -1]
  # Hudson River
  hud <- read.csv("Output/Data/Sites/csv/HudsonRiver/HudsonRiverObsPredPCB.csv")
  hud <- hud[, -c(1:2)]
  # Lake Washington
  lwa <- read.csv("Output/Data/Sites/csv/LakeWashington/LakeWashingtonObsPredPCB.csv")
  lwa <- lwa[, -c(1:2)]
  # Tributaries to Lake Michigan
  glt <- read.csv("Output/Data/Sites/csv/GreatLakes/TributariesObsPredPCB.csv")
  glt <- glt[, -c(1:2)]
  # New Bedford Harbor data
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHObsPredPCB.csv")
  nbh <- nbh[, -c(1:2)]
  # Passaic River
  pas <- read.csv("Output/Data/Sites/csv/PassaicRiver/PassaicObsPredPCB.csv")
  pas <- pas[,-c(1:2)]
  # Portland Harbord data
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/PortlandHarborObsPredPCB.csv")
  por <- por[, -c(1:2)]
  # Spokane River data
  spo <- read.csv("Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverObsPredPCB.csv")
  spo <- spo[, -c(1:2)]
  # Combine the data frames
  combined_data <- rbind(mic, bfc, che, fox, grl, glt, hud, lwa, nbh, pas,
                         por, spo)
}

# Create a custom color palette with 12 different colors
n_colors <- 12
custom_palette <- colorRampPalette(brewer.pal(9, "Set1"))(n_colors)

# Add a new column with the number of rows for each LocationName
combined_data$Location <- paste0(combined_data$LocationName,
                                 " (n=", ave(combined_data$LocationName,
                                             combined_data$LocationName,
                                             FUN = length), ")")                                                                         

# Plot prediction vs. observations, 1:1 line
CombinePredObsPlot <- ggplot(combined_data,
                             aes(x = 10^(observed), y = 10^(predicted),
                                 fill = Location)) +  # Use the new column for legend labels
  geom_point(shape = 21, size = 1.2, alpha = 0.5) +
  scale_y_log10(limits = c(0.001, 10^7), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.001, 10^7), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values = custom_palette) +
  xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration PCBi (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) +
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  theme_bw() +
  theme(aspect.ratio = 15/15,
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  annotation_logticks(sides = "bl")

# Print plot
print(CombinePredObsPlot)

# Save plot
ggsave("Output/Figures/Sites/CombineObsPredPCBi.png",
       plot = CombinePredObsPlot, width = 8, height = 8, dpi = 500)

