## Water PCB concentrations data analysis per site
## Housatonic River
## Aroclors 1254 and 1260, no congener analysis
## GE facility map @https://semspub.epa.gov/work/01/574882.pdf

# Install packages
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("robustbase")
install.packages("dplyr")
install.packages("tibble")
install.packages("Matrix")
install.packages("lme4")
install.packages("MuMIn")
install.packages("lmerTest")
install.packages("zoo")
install.packages("dataRetrieval")
install.packages("reshape")
install.packages("scales")

# Load libraries
{
  library(ggplot2)
  library(scales) # function trans_breaks
  library(stringr) # str_detect
  library(robustbase) # function colMedians
  library(dplyr) # performs %>%
  library(tibble) # adds a column
  library(lme4) # performs lme
  library(MuMIn) # gets Rs from lme
  library(lmerTest) # gets the p-value from lme
  library(zoo) # yields seasons
  library(dataRetrieval) # read data from USGS
  library(reshape)
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# Select Housatonic River data ---------------------------------------------------
hou <- wdc[str_detect(wdc$LocationName, 'Housatonic River'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  hou$SampleDate <- as.Date(hou$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(hou$SampleDate) - min(as.Date(hou$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- hou$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hou$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  hou.tpcb <- cbind(factor(hou$SiteID), hou$SampleDate,
                    hou$Latitude, hou$Longitude, as.matrix(hou$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(hou.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Get coordinates per site to plot in Google Earth
hou.location <- hou.tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
hou.location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                          data = hou.location, mean)

# General plots -------------------------------------------------------------------
# (1) Histograms
hist(hou.tpcb$tPCB)
hist(log10(hou.tpcb$tPCB))

# (2) Time trend plots
HRTime <- ggplot(hou.tpcb, aes(y = tPCB, x = format(date, '%Y'))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  xlab("") +
  scale_y_log10(
    breaks = c(10000, 100000, 1000000, 10000000),  # Specify the desired breaks
    labels = label_comma()(c(10000, 100000, 1000000, 10000000))  # Specify the desired labels
  ) +
  theme_classic() +
  ylab(expression(bold(Sigma*"PCB (pg/L)"))) +
  theme(
    axis.text.y = element_text(face = "bold", size = 22),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.x = element_text(size = 24, angle = 60, hjust = 1),
    axis.title.x = element_text(face = "bold", size = 17),
    plot.margin = margin(0, 0, 0, 0, unit = "cm"))

# Pint plot
print(HRTime)

# Save plot in folder
ggsave("Output/Plots/Sites/Temporal/plotHousRiverTimeV02.png",
       plot = HRTime, width = 18, height = 8, dpi = 500)

# (3) Seasonality
ggplot(hou.tpcb, aes(x = season, y = tPCB)) +
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
  annotate("text", x = 4, y = 10^5.6, label = "Housotonic River",
           size = 3)

# (4) Sites
ggplot(hou.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
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
  annotate("text", x = 32, y = 10^5.6, label = "Housotonic River",
           size = 3)

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Housatonic River
{
  siteHouN1 <- "01197000" # EAST BRANCH HOUSATONIC RIVER AT COLTSVILLE, MA
  siteHouN2 <- "01197500" # HOUSATONIC RIVER NEAR GREAT BARRINGTON, MA
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  # Retrieve USGS data
  flow.1 <- readNWISdv(siteHouN1, paramflow,
                       min(hou.tpcb$date), max(hou.tpcb$date))
  flow.2 <- readNWISdv(siteHouN2, paramflow,
                       min(hou.tpcb$date), max(hou.tpcb$date))
  # Add USGS data to hou.tpcb, matching dates (m3/s, 0.03 conversion factor)
  hou.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(hou.tpcb$date,
                                                     flow.1$Date)]
  hou.tpcb$flow.2 <- 0.03*flow.2$X_00060_00003[match(hou.tpcb$date,
                                                     flow.2$Date)]
}

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- hou.tpcb$tPCB
time <- hou.tpcb$time
flow <- hou.tpcb$flow.1
site <- hou.tpcb$site.code
season <- hou.tpcb$season

# tPCB vs. time + flow + season + site
lme.hou.tpcb <- lmer(log10(tpcb) ~ 1 + time + flow + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lme.hou.tpcb)

# Look at residuals
{
  res.hou.tpcb <- resid(lme.hou.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.hou.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.hou.tpcb)
}
# Shapiro test
shapiro.test(res.hou.tpcb)
# Lme does not provide a good model for the data.
