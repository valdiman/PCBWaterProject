
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
  mic <- read.csv("Output/Data/Sites/csv/21Mich/21MichRFObsPredPCB.csv")
  # Bannister Federal Complex
  bfc <- read.csv("Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexRFObsPredPCB.csv")
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFObsPredPCB.csv")
  # Fox River data
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverRFObsPredPCB.csv")
  # Great Lakes
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
  combined_data <- rbind(mic, bfc, che, fox, grl, hud, lwa, nbh, pas, por, spo)
}

# Create a custom color palette with 9 different colors
n_colors <- 11
custom_palette <- colorRampPalette(brewer.pal(9, "Set1"))(n_colors)

# Add a new column with the number of rows for each Location
combined_data$Location <- paste0(combined_data$Location,
                                 " (n=", ave(combined_data$Location,
                                             combined_data$Location,
                                             FUN = length), ")")                                                                         

# Plot prediction vs. observations, 1:1 line
CombinePredObsPlot <- ggplot(combined_data,
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
        legend.text = element_text(size = 18),
        legend.position = c(1.28, 0.78),  # Adjust the position (values are between 0 and 1)
        legend.background = element_rect(fill = "transparent")) +  # Make the legend background transparent
  annotation_logticks(sides = "bl")

# Print plot
print(CombinePredObsPlot)

# Save plot
ggsave("Output/Figures/Sites/CombineRFObsPredPCBi.png",
       plot = CombinePredObsPlot, width = 18, height = 8, dpi = 500)
