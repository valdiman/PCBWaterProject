## Water PCB concentrations plots per site
## Spokane River

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

# Select Spokane River data ---------------------------------------------------
spo <- wdc[str_detect(wdc$LocationName, 'Spokane River'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  spo$SampleDate <- as.Date(spo$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(spo$SampleDate) - min(as.Date(spo$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- spo$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(spo$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  spo.tpcb <- cbind(factor(spo$SiteID), spo$SampleDate, as.matrix(spo$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(spo.tpcb) <- c("SiteID", "date", "tPCB", "time",
                          "site.code", "season")
}

# General plots -------------------------------------------------------------------
# (1) Histograms
hist(spo.tpcb$tPCB)
hist(log10(spo.tpcb$tPCB))

# (2) Time trend plots
SPOTime <- ggplot(spo.tpcb, aes(y = tPCB, x = format(date, '%Y-%m'))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  xlab("") +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000, 10000),  # Specify the desired breaks
    labels = label_comma()(c(1, 10, 100, 1000, 10000))  # Specify the desired labels
  ) +
  theme_classic() +
  ylab(expression(bold(Sigma*"PCB (pg/L)"))) +
  theme(
    axis.text.y = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 20, angle = 60, hjust = 1),
    axis.title.x = element_text(face = "bold", size = 17),
    plot.margin = margin(0, 0, 0, 0, unit = "cm"))

# Print plot
print(SPOTime)

# Save plot in folder
ggsave("Output/Plots/Sites/Temporal/SpokaneTime.png",
       plot = SPOTime, width = 6, height = 5, dpi = 500)

# (3) Seasonality
ggplot(spo.tpcb, aes(x = season, y = tPCB)) +
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
              shape = 21, fill = "white") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotate("text", x = 1, y = 10^4.2, label = "Spokane River",
           size = 3)

# (4) Sites
ggplot(spo.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
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
              shape = 21, fill = "white") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotate("text", x = 2.5, y = 10^4.2, label = "Spokane River",
           size = 3)

# Remove site -------------------------------------------------------------
## Sample sites not located at the Spokane River
{
  spo.tpcb.1 <- subset(spo.tpcb, SiteID != c("WCPCB-SPR002")) # City of Spokane WRF
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR005")) # Regional WRF
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR006")) # Inland Empire paper
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR008")) # Kaiser Aluminum
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR010")) # Liberty Lake sewer
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR011")) # Post Falls WWTP
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR013")) # Coeur d'Alene WWTP
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR015")) # Hangman Creek
}

# (1) Histograms
hist(spo.tpcb.1$tPCB)
hist(log10(spo.tpcb.1$tPCB))

# (2) Time trend plots
SPOTimeV02 <- ggplot(spo.tpcb.1, aes(y = tPCB, x = format(date, '%Y-%m'))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  xlab("") +
  scale_y_log10(
    breaks = c(10, 100, 1000),  # Specify the desired breaks
    labels = label_comma()(c(10, 100, 1000))  # Specify the desired labels
  ) +
  theme_classic() +
  ylab(expression(bold(Sigma*"PCB (pg/L)"))) +
  theme(
    axis.text.y = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 20, angle = 60, hjust = 1),
    axis.title.x = element_text(face = "bold", size = 17),
    plot.margin = margin(0, 0, 0, 0, unit = "cm"))

# Print plot
print(SPOTimeV02)

# Save plot in folder
ggsave("Output/Plots/Sites/Temporal/SpokaneTimeV02.png",
       plot = SPOTimeV02, width = 6, height = 5, dpi = 500)

# (3) Sites
ggplot(spo.tpcb.1, aes(x = factor(SiteID), y = tPCB)) + 
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
  annotate("text", x = 6, y = 10^3, label = "Spokane River",
           size = 3)
