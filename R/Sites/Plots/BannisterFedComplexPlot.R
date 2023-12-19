## Water PCB concentrations plots per site
## Blue River

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

# Select Bannister Federal Complex data ---------------------------------------------------
bfc <- wdc[str_detect(wdc$LocationName, 'Bannister Fed Complex'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  bfc$SampleDate <- as.Date(bfc$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(bfc$SampleDate) - min(as.Date(bfc$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- bfc$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(bfc$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  bfc.tpcb <- cbind(factor(bfc$SiteID), bfc$SampleDate, as.matrix(bfc$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(bfc.tpcb) <- c("SiteID", "date", "tPCB", "time",
                          "site.code", "season")
}

# General plots -------------------------------------------------------------------
# (1) Histograms
hist(bfc.tpcb$tPCB)
hist(log10(bfc.tpcb$tPCB))

# (2) Time trend plots
BCTime <- ggplot(bfc.tpcb, aes(y = tPCB, x = format(date, '%Y'))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  xlab("") +
  scale_y_log10(
    breaks = c(100, 1000, 10000, 100000),  # Specify the desired breaks
    labels = label_comma()(c(100, 1000, 10000, 100000))  # Specify the desired labels
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
print(BCTime)

# Save plot in folder
ggsave("Output/Plots/Sites/Temporal/BannisterFedComplexTime.png",
       plot = BCTime, width = 6, height = 5, dpi = 500)

# (4) Sites
ggplot(bfc.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
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
  annotate("text", x = 6.5, y = 10^5, label = "Blue River",
           size = 3)
