## Water PCB concentrations plots per site
## Newtown Creek

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

# Select 21 Michigan data ---------------------------------------------------
ntc <- wdc[str_detect(wdc$LocationName, 'Newtown Creek'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  ntc$SampleDate <- as.Date(ntc$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(ntc$SampleDate) - min(as.Date(ntc$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- ntc$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(ntc$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  ntc.tpcb <- cbind(factor(ntc$SiteID), ntc$SampleDate, as.matrix(ntc$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(ntc.tpcb) <- c("SiteID", "date", "tPCB", "time",
                          "site.code", "season")
}

# General plots -------------------------------------------------------------------
# (1) Histograms
hist(ntc.tpcb$tPCB)
hist(log10(ntc.tpcb$tPCB))

# (2) Time trend plots
NTCTime <- ggplot(ntc.tpcb, aes(y = tPCB, x = format(date, '%Y-%m'))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  xlab("") +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000, 10000, 100000),  # Specify the desired breaks
    labels = label_comma()(c(1, 10, 100, 1000, 10000, 100000))  # Specify the desired labels
  ) +
  theme_classic() +
  ylab(expression(bold(Sigma*"PCB (pg/L)"))) +
  theme(
    axis.text.y = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 20, angle = 60, hjust = 1),
    axis.title.x = element_text(face = "bold", size = 17),
    plot.margin = margin(0.1, 0, 0, 0, unit = "cm"))

# Print plot
print(NTCTime)

# Save plot in folder
ggsave("Output/Plots/Sites/Temporal/NewtownCreekTime.png",
       plot = NTCTime, width = 6, height = 5, dpi = 500)

# (3) Seasonality
ggplot(ntc.tpcb, aes(x = season, y = tPCB)) +
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
  annotate("text", x = 1, y = 500, label = "Newtown Creek",
           size = 3)

# (4) Sites
ggplot(ntc.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concetration",
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
  annotate("text", x = 11, y = 50000, label = "Newtown Creek",
           size = 3)

