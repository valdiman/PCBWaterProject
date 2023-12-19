## Water PCB concentrations plots per site
## Richardson Hill Road Landfill
## Arolclor method

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

# Select Richardson Hill Road Landfill data ---------------------------------------------------
rhr <- wdc[str_detect(wdc$LocationName, 'Richardson Hill Road Landfill'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  rhr$SampleDate <- as.Date(rhr$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(rhr$SampleDate) - min(as.Date(rhr$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- rhr$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(rhr$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  rhr.tpcb <- cbind(factor(rhr$SiteID), rhr$SampleDate, as.matrix(rhr$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(rhr.tpcb) <- c("SiteID", "date", "tPCB", "time",
                          "site.code", "season")
}

# General plots -------------------------------------------------------------------
# (1) Histograms
hist(rhr.tpcb$tPCB)
hist(log10(rhr.tpcb$tPCB))

# (2) Time trend plots
RHTime <- ggplot(rhr.tpcb, aes(y = tPCB, x = format(date, '%Y-%m'))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  xlab("") +
  scale_y_log10(
    breaks = c(150000, 400000, 1000000, 3000000),  # Specify the desired breaks
    labels = label_comma()(c(150000, 400000, 1000000, 3000000))  # Specify the desired labels
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
print(RHTime)

# Save plot in folder
ggsave("Output/Plots/Sites/Temporal/RichardsonHillTime.png",
       plot = RHTime, width = 6, height = 5, dpi = 500)

# (4) Sites
ggplot(rhr.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
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
  annotate("text", x = 6.5, y = 10^4.5, label = "Richardson Hill Road Landfill",
           size = 3)
