
# Install packages
install.packages("ggplot2")
install.packages("scales")

# Load libraries
{
  library(ggplot2)
  library(scales) # function trans_breaks
}

# Read generated data
{
  # 21 Mich
  mic <- read.csv("Output/Data/Sites/csv/21Mich/ObsPred21MichPCB.csv")
  mic <- mic[, -c(1:2)]
  # Blue River
  blr <- read.csv("Output/Data/Sites/csv/BlueRiver/ObsPredBlueRiverPCB.csv")
  blr <- blr[, -c(1:2)]
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ObsPredChesapeakeBayPCB.csv")
  che <- che[, -c(1:2)]
  # Fox River data
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/ObsPredFoxRiverPCB.csv")
  fox <- fox[, -c(1:2)]
  # Great Lakes
  grl <- read.csv("Output/Data/Sites/csv/GreatLakes/ObsPredGreatLakesPCB.csv")
  grl <- grl[, -1]
  # Tributaries Great Lakes
  glt <- read.csv("Output/Data/Sites/csv/GreatLakes/ObsPredTributariesPCB.csv")
  glt <- glt[, -c(1:2)]
  # Hudson River
  hud <- read.csv("Output/Data/Sites/csv/HudsonRiver/ObsPredHudsonRiverPCB.csv")
  hud <- hud[, -c(1:2)]
  # New Bedford Harbor data
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/ObsPredNewBHPCB.csv")
  nbh <- nbh[, -c(1:2)]
  # Passaic River
  pas <- read.csv("Output/Data/Sites/csv/PassaicRiver/ObsPredPassaicPCB.csv")
  pas <- pas[,-c(1:2)]
  # Portland Harbord data
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/ObsPredPortlandHarborPCB.csv")
  por <- por[, -c(1:2)]
  # Spokane River data
  spo <- read.csv("Output/Data/Sites/csv/SpokaneRiver/ObsPredSpokaneRiverPCB.csv")
  spo <- spo[, -c(1:2)]
  # Combine the data frames
  combined_data <- rbind(mic, blr, che, fox, grl, glt, hud, nbh, pas, por, spo)
}

# Create a custom color palette with 11 different colors
custom_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                   "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                   "#999999", "#66C2A5", "#FFA07A")

# Add a new column with the number of rows for each LocationName
combined_data <- transform(combined_data,
                           Location = paste(LocationName,
                                            " (n =", ave(LocationName,
                                                         LocationName,
                                                         FUN = length), ")"))

# Plot prediction vs. observations, 1:1 line
CombinePredObsPlot <- ggplot(combined_data,
                             aes(x = 10^(observed), y = 10^(predicted),
                                 fill = Location)) +  # Use the new column for legend labels
  geom_point(shape = 21, size = 1.5, alpha = 0.5) +
  scale_y_log10(limits = c(0.001, 10^7), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.001, 10^7), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values = custom_colors) +
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
        legend.text = element_text(size = 12))  +
  annotation_logticks(sides = "bl")

# Print plot
print(CombinePredObsPlot)

# Save plot
ggsave("Output/Figures/Sites/CombineObsPredPCBi.png",
       plot = CombinePredObsPlot, width = 8, height = 8, dpi = 500)



