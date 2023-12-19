## Water PCB concentrations plots per site
## Kalamazoo River
## Source: Allied Paper, Inc. Recycling of carbonless copy paper
## Aroclors 1242, 1254 and 1260, no congener analysis

# Install packages
install.packages("ggplot2")
install.packages("scales")
install.packages("stringr")
install.packages("zoo")

# Load libraries
{
  library(ggplot2)
  library(scales) # function trans_breaks
  library(stringr) # str_detect
  library(zoo) # yields seasons
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# Select Kalamazoo data ---------------------------------------------------
kal <- wdc[str_detect(wdc$LocationName, 'Kalamazoo River'),]
# Superfund site from Morrow Dam (Kalamazoo River) to Lake Michigan
# and 30 miles of Portage Creek (south), Cork St and Portage Creek Cork St sites
# Dredging occurred at Plainwell Dam site.

# Data preparation --------------------------------------------------------
{
  # Change date format
  kal$SampleDate <- as.Date(kal$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(kal$SampleDate) - min(as.Date(kal$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- kal$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(kal$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  kal.tpcb <- cbind(factor(kal$SiteID), kal$SampleDate, as.matrix(kal$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(kal.tpcb) <- c("SiteID", "date", "tPCB", "time",
                          "site.code", "season")
}

# General plots -------------------------------------------------------------------
# (1) Histograms
# tPCB
hist(kal.tpcb$tPCB)
hist(log10(kal.tpcb$tPCB))

# (2) Time trend plots
KRTime <- ggplot(kal.tpcb, aes(y = tPCB, x = format(date, '%Y'))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  xlab("") +
  scale_y_log10(
    breaks = c(100, 1000, 10000, 100000, 1000000),  # Specify the desired breaks
    labels = label_comma()(c(100, 1000, 10000, 100000, 1000000))  # Specify the desired labels
  ) +
  theme_classic() +
  ylab(expression(bold(Sigma*"PCB (pg/L)"))) +
  theme(
    axis.text.y = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 20, angle = 60, hjust = 1),
    axis.title.x = element_text(face = "bold", size = 17),
    plot.margin = margin(0, 0, 0, 0, unit = "cm"))

# Print the plot
print(KRTime)

# Save plot in folder
ggsave("Output/Plots/Sites/Temporal/KalamazooRiverTime.png",
       plot = KRTime, width = 7, height = 5, dpi = 500)

# (3) Seasonality
ggplot(kal.tpcb, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotate("text", x = 1.2, y = 10^6.5, label = "Kalamazoo River",
           size = 3)

# (4) Sites
ggplot(kal.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotate("text", x = 5, y = 10^6.5, label = "Kalamazoo River",
           size = 3)

# Remove site -------------------------------------------------------------
# Remove site PlainwellDam. Dredging = WCPCB-KAL023
kal.tpcb.1 <- subset(kal.tpcb, SiteID != c("WCPCB-KAL023"))

# 2nd plots -----------------------------------------------------------
# (1.1) tPCB
hist(kal.tpcb.1$tPCB)
hist(log10(kal.tpcb.1$tPCB))

# (2) Time trend plots
ggplot(kal.tpcb.1, aes(y = tPCB,
                       x = format(date,'%Y'))) +
  geom_point(shape = 21, fill = "#66ccff") +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9)) +
  annotate("text", x = 2, y = 100, label = "Kalamazoo River",
           size = 3)

