## Script to combine maps of avereage tPCB and temporal plots of tPCB.

# Install packages
install.packages("ggplot2")
install.packages("png")
install.packages("grid")

# Load libraries
{
  library(ggplot2)
  library(png)
  library(grid)
}

# Fox River ---------------------------------------------------------------
# Read the PNG images
mapFR <- readPNG("Output/Maps/Sites/maptPCBFoxRiverAveV01.png")
plotFR <- readPNG("Output/Plots/Sites/Temporal/plotFoxRiverTime.png")

# Create a plot with a smaller canvas
p <- ggplot() +
  xlim(0, 6) + ylim(0, 6) +  # Adjust the limits as needed
  theme_void()

# Add the images at specific positions with reduced size
p <- p +
  annotation_custom(rasterGrob(mapFR),
                    xmin = 0, xmax = 5, ymin = 1, ymax = 6) +
  annotation_custom(rasterGrob(plotFR, width = unit(6*0.5, "cm"),
                               height = unit(5*0.5, "cm")),
                    xmin = 3.25, xmax = 5.25, ymin = 1.5, ymax = 3.9)

# Print the plot
print(p)

# Save plot directly to a PNG file
ggsave("Output/Figures/Sites/FoxRiver.png",
       plot = p, width = 6, height = 4, dpi = 1000)

# Housatonic River --------------------------------------------------------
# Read the PNG images
mapHR <- readPNG("Output/Maps/Sites/maptPCBHousatonicRiverAveV02.png")
plotHR <- readPNG("Output/Plots/Sites/Temporal/plotHousRiverTimeV02.png")

# Create a plot with a smaller canvas
p <- ggplot() +
  xlim(0, 6) + ylim(0, 6) +  # Adjust the limits as needed
  theme_void()

# Add the images at specific positions with reduced size
p <- p +
  annotation_custom(rasterGrob(mapHR),
                    xmin = 0, xmax = 5, ymin = 1, ymax = 6) +
  annotation_custom(rasterGrob(plotHR, width = unit(18*0.38, "cm"),
                               height = unit(8*0.38, "cm")),
                    xmin = 3.64, xmax = 5.64, ymin = 1.5, ymax = 3.9)

# Print the plot
print(p)

# Save plot directly to a PNG file
ggsave("Output/Figures/Sites/HousatonicRiverV02.png",
       plot = p, width = 6, height = 4, dpi = 1000)

# Portland Harbor ---------------------------------------------------------
# Read the PNG images
mapPO <- readPNG("Output/Maps/Sites/maptPCBPortlandAveV01.png")
plotPO <- readPNG("Output/Plots/Sites/Temporal/plotPortlandTime.png")

# Create a plot with a smaller canvas
p <- ggplot() +
  xlim(0, 6) + ylim(0, 6) +  # Adjust the limits as needed
  theme_void()

# Add the images at specific positions with reduced size
p <- p +
  annotation_custom(rasterGrob(mapPO),
                    xmin = 0, xmax = 5, ymin = 1, ymax = 6) +
  annotation_custom(rasterGrob(plotPO, width = unit(6*0.6, "cm"),
                               height = unit(5*0.6, "cm")),
                    xmin = 3.2, xmax = 5.2, ymin = 1.5, ymax = 3.9)

# Print the plot
print(p)

# Save plot directly to a PNG file
ggsave("Output/Figures/Sites/PortlandHarbor.png",
       plot = p, width = 6, height = 4, dpi = 1000)

# Hudson River ------------------------------------------------------------
# Read the PNG images
mapHU <- readPNG("Output/Maps/Sites/maptPCBHudsonRiverAveV01.png")
plotHU <- readPNG("Output/Plots/Sites/Temporal/plotHudsonRiverTime.png")

# Create a plot with a smaller canvas
p <- ggplot() +
  xlim(0, 6) + ylim(0, 6) +  # Adjust the limits as needed
  theme_void()

# Add the images at specific positions with reduced size
p <- p +
  annotation_custom(rasterGrob(mapHU),
                    xmin = 0, xmax = 5, ymin = 1, ymax = 6) +
  annotation_custom(rasterGrob(plotHU, width = unit(8*0.65, "cm"),
                               height = unit(5*0.65, "cm")),
                    xmin = 2.75, xmax = 4.75, ymin = 0.9, ymax = 3.4)

# Print the plot
print(p)

# Save plot directly to a PNG file
ggsave("Output/Figures/Sites/HudsonRiver.png",
       plot = p, width = 6, height = 4, dpi = 1000)

# Kalamazoo River ---------------------------------------------------------
# Read the PNG images
mapKR <- readPNG("Output/Maps/Sites/maptPCBKalamazooAveV01.png")
plotKR <- readPNG("Output/Plots/Sites/Temporal/plotKalRiverTime.png")

# Create a plot with a smaller canvas
p <- ggplot() +
  xlim(0, 6) + ylim(0, 6) +  # Adjust the limits as needed
  theme_void()

# Add the images at specific positions with reduced size
p <- p +
  annotation_custom(rasterGrob(mapKR),
                    xmin = 0, xmax = 5, ymin = 1, ymax = 6) +
  annotation_custom(rasterGrob(plotKR, width = unit(7*0.4, "cm"),
                               height = unit(5*0.4, "cm")),
                    xmin = 3.2, xmax = 5.2, ymin = 1.7, ymax = 4.2)

# Print the plot
print(p)

# Save plot directly to a PNG file
ggsave("Output/Figures/Sites/KalamazooRiver.png",
       plot = p, width = 6, height = 4, dpi = 1000)

# New Bedford Harbor ------------------------------------------------------
# Read the PNG images
mapNBH <- readPNG("Output/Maps/Sites/maptPCBNBHAveV01.png")
plotNBH <- readPNG("Output/Plots/Sites/Temporal/plotNBHTime.png")

# Create a plot with a smaller canvas
p <- ggplot() +
  xlim(0, 6) + ylim(0, 6) +  # Adjust the limits as needed
  theme_void()

# Add the images at specific positions with reduced size
p <- p +
  annotation_custom(rasterGrob(mapNBH),
                    xmin = 0, xmax = 5, ymin = 1, ymax = 6) +
  annotation_custom(rasterGrob(plotNBH, width = unit(7*0.5, "cm"),
                               height = unit(6*0.5, "cm")),
                    xmin = 3.1, xmax = 5.1, ymin = 1.6, ymax = 4.1)

# Print the plot
print(p)

# Save plot directly to a PNG file
ggsave("Output/Figures/Sites/NewBedfordHarbor.png",
       plot = p, width = 6, height = 4, dpi = 1000)

# Spokane River -----------------------------------------------------------
# Read the PNG images
mapSPO <- readPNG("Output/Maps/Sites/maptPCBSpoAveV01.png")
plotSPO <- readPNG("Output/Plots/Sites/Temporal/plotSpokaneTime.png")

# Create a plot with a smaller canvas
p <- ggplot() +
  xlim(0, 6) + ylim(0, 6) +  # Adjust the limits as needed
  theme_void()

# Add the images at specific positions with reduced size
p <- p +
  annotation_custom(rasterGrob(mapSPO),
                    xmin = 0, xmax = 5, ymin = 1, ymax = 6) +
  annotation_custom(rasterGrob(plotSPO, width = unit(6*0.4, "cm"),
                               height = unit(5*0.4, "cm")),
                    xmin = 3.45, xmax = 5.45, ymin = 1.7, ymax = 4.2)

# Print the plot
print(p)

# Save plot directly to a PNG file
ggsave("Output/Figures/Sites/SpokaneRiver.png",
       plot = p, width = 6, height = 4, dpi = 1000)

# Chesapeake Bay ----------------------------------------------------------
# Read the PNG images
mapChe <- readPNG("Output/Maps/Sites/maptPCBCheAveV01.png")
plotChe <- readPNG("Output/Plots/Sites/Temporal/plotChesapeakeTime.png")

# Create a plot with a smaller canvas
p <- ggplot() +
  xlim(0, 6) + ylim(0, 6) +  # Adjust the limits as needed
  theme_void()

# Add the images at specific positions with reduced size
p <- p +
  annotation_custom(rasterGrob(mapChe),
                    xmin = 0, xmax = 5, ymin = 1, ymax = 6) +
  annotation_custom(rasterGrob(plotChe, width = unit(6*0.61, "cm"),
                               height = unit(5*0.61, "cm")),
                    xmin = 2.95, xmax = 4.95, ymin = 1.3, ymax = 3.8)

# Print the plot
print(p)

# Save plot directly to a PNG file
ggsave("Output/Figures/Sites/ChesapeakeBay.png",
       plot = p, width = 6, height = 4, dpi = 1000)
