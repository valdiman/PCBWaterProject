
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
  # Anacostia River
  anr <- read.csv("Output/Data/Sites/csv/AnacostiaRiver/AnacostiaRiverObsPredtPCB.csv")
  anr <- anr[, -1]
  # Bannister Federal Complex
  bfc <- read.csv("Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexObsPredtPCB.csv")
  bfc <- bfc[, -1]
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeObsPredtPCB.csv")
  che <- che[, -1]
  # Fox River data
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverObsPredtPCB.csv")
  fox <- fox[, -1]
  # Kalamazoo River data
  kal <- read.csv("Output/Data/Sites/csv/KalamazooRiver/KalamazooObsPredtPCB.csv")
  kal <- kal[, -1]
  # Lake Washington
  lwa <- read.csv("Output/Data/Sites/csv/LakeWashington/LakeWashingtonObsPredtPCB.csv")
  lwa <- lwa[, -1]
  # New Bedford Harbor data
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHObsPredtPCBV02.csv")
  nbh <- nbh[, -1]
  # Portland Harbord data
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/PortlandHarborObsPredtPCB.csv")
  por <- por[, -1]
  # Richardson Hill Road Landfill
  rhr <- read.csv("Output/Data/Sites/csv/Richardson/RichardsonObsPredtPCB.csv")
  rhr <- rhr[, -1]
  # Spokane River data
  spo <- read.csv("Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverObsPredtPCB.csv")
  spo <- spo[, -1]
  # Combine the data frames
  combined_data <- rbind(anr, bfc, che, fox, kal, lwa, nbh, por, rhr, spo)
}

# Create a custom color palette with 10 different colors
n_colors <- 10
custom_palette <- colorRampPalette(brewer.pal(9, "Set1"))(n_colors)

# Add a new column with the number of rows for each Location
combined_data$Location <- paste0(combined_data$Location,
                                 " (n=", ave(combined_data$Location,
                                             combined_data$Location,
                                             FUN = length), ")")   

# Plot prediction vs. observations, 1:1 line
CombinePredObsPlot <- ggplot(combined_data,
                             aes(x = tPCB, y = predicted, fill = Location)) +
  geom_point(shape = 21, size = 4, alpha = 0.5) +
  scale_y_log10(limits = c(1, 10^8), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^8), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values = custom_palette) +  # Use custom color palette
  xlab(expression(bold("Observed concentration " * Sigma * "PCB (pg/L)"))) +
  ylab(expression(bold("Predicted LME concentration " * Sigma * "PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15,
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = c(1.35, 0.8),  # Adjust the position (values are between 0 and 1)
        legend.background = element_rect(fill = "transparent")) +  # Make the legend background transparent
  annotation_logticks(sides = "bl")

# Print plot
print(CombinePredObsPlot)

# Save plot
ggsave("Output/Figures/Sites/CombineObsPredtPCBV02.png",
       plot = CombinePredObsPlot, width = 18, height = 8, dpi = 500)

