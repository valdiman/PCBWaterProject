
# Install packages
install.packages("ggplot2")
install.packages("scales")
install.packages("RColorBrewer")

# Load libraries
{
  library(ggplot2)
  library(scales) # function trans_breaks
  library(RColorBrewer)
}

# Read generated data
{
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakePredic_Obser.csv")
  che <- che[, -1]
  # Fox River data
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverPredic_Obser.csv")
  fox <- fox[, -1]
  # New Bedford Harbor data
  kal <- read.csv("Output/Data/Sites/csv/KalamazooRiver/KalamazooPredic_Obser.csv")
  kal <- kal[, -1]
  # New Bedford Harbor data
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHPredic_Obser.csv")
  nbh <- nbh[, -1]
  # Portland Harbord data
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/PortlandPredic_Obser.csv")
  por <- por[, -1]
  # Spokane River data
  spo <- read.csv("Output/Data/Sites/csv/SpokaneRiver/SpokanePredic_Obser.csv")
  spo <- spo[, -1]
  # Combine the data frames
  combined_data <- rbind(che, fox, kal, nbh, por, spo)
}

# Create a custom color palette with distinct colors for the 6 locations
custom_colors <- brewer.pal(6, "Set1")

# Plot prediction vs. observations, 1:1 line
CombinePredObsPlot <- ggplot(combined_data, aes(x = tPCB, y = predicted, fill = Location)) +
  geom_point(shape = 21, size = 2) +
  scale_y_log10(limits = c(1, 10^8), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^8), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values = custom_colors) +  # Use custom color palette
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line (factor of 2)
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
ggsave("Output/Figures/Sites/CombineObsPredtPCB.png",
       plot = CombinePredObsPlot, width = 8, height = 8, dpi = 500)

