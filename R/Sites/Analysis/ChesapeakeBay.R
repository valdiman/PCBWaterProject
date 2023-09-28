## Water PCB concentrations data analysis per site
## Chesapeake Bay & Delaware Canal

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
install.packages("tidyr")
install.packages('patchwork')
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
  library(tidyr) # function gather
  library(patchwork) # combine plots
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# Select Chesapeake Bay & Delaware Canal data ---------------------------------------------------
che <- wdc[str_detect(wdc$LocationName, 'Chesapeake Bay'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  che$SampleDate <- as.Date(che$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(che$SampleDate) - min(as.Date(che$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- che$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(che$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  che.tpcb <- cbind(factor(che$SiteID), che$SampleDate,
                    che$Latitude, che$Longitude, as.matrix(che$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(che.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Get coordinates per site to plot in Google Earth
che.location <- che.tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
che.location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                          data = che.location, mean)

# General plots -------------------------------------------------------------------
# (1) Histograms
hist(che.tpcb$tPCB)
hist(log10(che.tpcb$tPCB))

# (2) Time trend plots
CHTime <- ggplot(che.tpcb, aes(y = tPCB, x = format(date, '%Y'))) +
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
print(CHTime)

# Save plot in folder
ggsave("Output/Plots/Sites/Temporal/plotChesapeakeTime.png",
       plot = CHTime, width = 6, height = 5, dpi = 500)

# (3) Seasonality
ggplot(che.tpcb, aes(x = season, y = tPCB)) +
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
  annotate("text", x = 1, y = 20, label = "Chesapeake Bay",
           size = 3)

# (4) Sites
ggplot(che.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
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
  annotate("text", x = 5, y = 20, label = "Chesapeake Bay",
           size = 3)

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- che.tpcb$tPCB
time <- che.tpcb$time
site <- che.tpcb$site.code
season <- che.tpcb$season
# tPCB vs. time + season + flow + temp + site
lme.che.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.che.tpcb)
# Look at residuals
{
  res.che.tpcb <- resid(lme.che.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/ChesapeakeBayQ-QtPCB.pdf")
  qqnorm(res.che.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.che.tpcb)
  dev.off()
}
# Shapiro test
shapiro.test(resid(lme.che.tpcb))

# Remove minimum and bottom 2 values, ~ 10 pg/L
che.tpcb.1 <- subset(che.tpcb,tPCB > 10)

# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- che.tpcb.1$tPCB
time <- che.tpcb.1$time
site <- che.tpcb.1$site.code
season <- che.tpcb.1$season
# tPCB vs. time + season + flow + temp + site
lme.che.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.che.tpcb)
# Look at residuals
{
  res.che.tpcb <- resid(lme.che.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/ChesapeakeBayQ-QtPCBV02.pdf")
  qqnorm(res.che.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.che.tpcb)
  dev.off()
}
# Shapiro test
shapiro.test(resid(lme.che.tpcb))

# Create matrix to store results
{
  lme.tpcb <- matrix(nrow = 1, ncol = 21)
  lme.tpcb[1] <- fixef(lme.che.tpcb)[1] # intercept
  lme.tpcb[2] <- summary(lme.che.tpcb)$coef[1,"Std. Error"] # intercept error
  lme.tpcb[3] <- summary(lme.che.tpcb)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb[4] <- fixef(lme.che.tpcb)[2] # time
  lme.tpcb[5] <- summary(lme.che.tpcb)$coef[2,"Std. Error"] # time error
  lme.tpcb[6] <- summary(lme.che.tpcb)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb[7] <- fixef(lme.che.tpcb)[3] # season 1
  lme.tpcb[8] <- summary(lme.che.tpcb)$coef[3,"Std. Error"] # season 1 error
  lme.tpcb[9] <- summary(lme.che.tpcb)$coef[3,"Pr(>|t|)"] # season 1 p-value
  lme.tpcb[10] <- fixef(lme.che.tpcb)[4] # season 2
  lme.tpcb[11] <- summary(lme.che.tpcb)$coef[4,"Std. Error"] # season 2 error
  lme.tpcb[12] <- summary(lme.che.tpcb)$coef[4,"Pr(>|t|)"] # season 2 p-value
  lme.tpcb[13] <- fixef(lme.che.tpcb)[5] # season 3
  lme.tpcb[14] <- summary(lme.che.tpcb)$coef[5,"Std. Error"] # season 3 error
  lme.tpcb[15] <- summary(lme.che.tpcb)$coef[5,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb[16] <- -log(2)/lme.tpcb[4]/365 # t0.5
  lme.tpcb[17] <- abs(-log(2)/lme.tpcb[4]/365)*lme.tpcb[5]/abs(lme.tpcb[4]) # t0.5 error
  lme.tpcb[18] <- as.data.frame(VarCorr(lme.che.tpcb))[1,'sdcor']
  lme.tpcb[19] <- as.data.frame(r.squaredGLMM(lme.che.tpcb))[1, 'R2m']
  lme.tpcb[20] <- as.data.frame(r.squaredGLMM(lme.che.tpcb))[1, 'R2c']
  lme.tpcb[21] <- shapiro.test(resid(lme.che.tpcb))$p.value
}

# Just 3 significant figures
lme.tpcb <- formatC(signif(lme.tpcb, digits = 3))
# Add column names
colnames(lme.tpcb) <- c("Intercept", "Intercept.error",
                        "Intercept.pv", "time", "time.error", "time.pv",
                        "season1", "season1.error", "season1.pv", "season2",
                        "season2.error", "season2, pv", "season3",
                        "season3.error", "season3.pv", "t05", "t05.error",
                        "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality")

# Export results
write.csv(lme.tpcb, file = "Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeLmetPCB.csv")

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.che.tpcb <- as.data.frame(fitted(lme.che.tpcb))
# Add column name
colnames(fit.lme.values.che.tpcb) <- c("predicted")
# Add predicted values to data.frame
che.tpcb.1$predicted <- 10^(fit.lme.values.che.tpcb$predicted)
# Create overall plot prediction vs. observations
predic.obs <- data.frame(tPCB = che.tpcb.1$tPCB, predicted = che.tpcb.1$predicted)
predic.obs <- data.frame(Location = che$LocationName[1], predic.obs)
# Save new data
write.csv(predic.obs, "Output/Data/Sites/csv/ChesapeakeBay/ChesapeakePredic_Obser.csv")

# Plot prediction vs. observations, 1:1 line
p <- ggplot(che.tpcb.1, aes(x = tPCB, y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10, 10^5.5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10, 10^5.5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 75, y = 10^5,
           label = expression(atop(" Chesapeake Bay (R"^2*"= 0.49)",
                                   paste("t"[1/2]*" = 15 ± 4.6 (yr)"))),
           size = 4, fontface = 2)
# See plot
print(p)
# Save plot
ggsave("Output/Plots/Sites/ObsPred/ChesapeakeBay/ChesapeakeBayObsPredtPCB.png",
       plot = p, width = 8, height = 8, dpi = 500)

# Plot residuals vs. predictions
{
  # Open a PNG graphics device
  png("Output/Plots/Sites/Residual/res_plotlmeChesapeakeBayResidualtPCB.png", width = 800,
      height = 600)
  # Create your plot
  plot(che.tpcb.1$predicted, resid(lme.che.tpcb),
       points(che.tpcb.1$predicted, resid(lme.che.tpcb), pch = 16, 
              col = "white"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  # Add lines to the plot
  abline(0, 0)
  abline(h = c(-2, 2), col = "grey")
  abline(v = seq(0, 30000, 5000), col = "grey")
  # Close the PNG graphics device
  dev.off()
}

# Estimate a factor of 2 between observations and predictions
che.tpcb.1$factor2 <- che.tpcb.1$tPCB/che.tpcb.1$predicted
factor2.tpcb <- nrow(che.tpcb.1[che.tpcb.1$factor2 > 0.5 & che.tpcb.1$factor2 < 2,
                                ])/length(che.tpcb.1[,1])*100

# Individual PCB Analysis -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  che.pcb <- subset(che, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  che.pcb <- subset(che.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  che.pcb <- log10(che.pcb)
  # Replace -inf to NA
  che.pcb <- do.call(data.frame,
                     lapply(che.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  che.pcb.1 <- che.pcb[,
                       -which(colSums(is.na(che.pcb))/nrow(che.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(che$SiteID)
  # Change date format
  SampleDate <- as.Date(che$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Create individual code for each site sampled
  site.numb <- che$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(che$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to fox.pcb.1
  che.pcb.1 <- cbind(che.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     site.numb, season.s)
  # Remove metadata
  che.pcb.2 <- subset(che.pcb.1, select = -c(SiteID:season.s))
}

# LME for individual PCBs -------------------------------------------------
# Get covariates
time <- che.pcb.1$time
season <- che.pcb.1$season
site <- che.pcb.1$site.numb

# Create matrix to store results
lme.pcb <- matrix(nrow = length(che.pcb.2[1,]), ncol = 21)

# Perform LME
for (i in 1:length(che.pcb.2[1,])) {
  fit <- lmer(che.pcb.2[,i] ~ 1 + time + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"))
  lme.pcb[i,1] <- fixef(fit)[1] # intercept
  lme.pcb[i,2] <- summary(fit)$coef[1,"Std. Error"] # intercept error
  lme.pcb[i,3] <- summary(fit)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.pcb[i,4] <- fixef(fit)[2] # time
  lme.pcb[i,5] <- summary(fit)$coef[2,"Std. Error"] # time error
  lme.pcb[i,6] <- summary(fit)$coef[2,"Pr(>|t|)"] # time p-value
  lme.pcb[i,7] <- fixef(fit)[3] # # season 1
  lme.pcb[i,8] <- summary(fit)$coef[3,"Std. Error"] # season 1 error
  lme.pcb[i,9] <- summary(fit)$coef[3,"Pr(>|t|)"] # # season 1 p-value
  lme.pcb[i,10] <- fixef(fit)[4] # season 2
  lme.pcb[i,11] <- summary(fit)$coef[4,"Std. Error"] # season 2 error
  lme.pcb[i,12] <- summary(fit)$coef[4,"Pr(>|t|)"] # season 2 p-value
  lme.pcb[i,13] <- fixef(fit)[5] # season 3
  lme.pcb[i,14] <- summary(fit)$coef[5,"Std. Error"] # season 3 error
  lme.pcb[i,15] <- summary(fit)$coef[5,"Pr(>|t|)"] # season 3 p-value
  lme.pcb[i,16] <- -log(2)/lme.pcb[i,4]/365 # t0.5
  lme.pcb[i,17] <- abs(-log(2)/lme.pcb[i,4]/365)*lme.pcb[i,5]/abs(lme.pcb[i,4]) # t0.5 error
  lme.pcb[i,18] <- as.data.frame(VarCorr(fit))[1,'sdcor']
  lme.pcb[i,19] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
  lme.pcb[i,20] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
  lme.pcb[i,21] <- shapiro.test(resid(fit))$p.value
}

# Just 3 significant figures
lme.pcb <- formatC(signif(lme.pcb, digits = 3))
# Add congener names
congeners <- colnames(che.pcb.2)
lme.pcb <- as.data.frame(cbind(congeners, lme.pcb))
# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "season1", "season1.error", "season1.pv", "season2",
                       "season2.error", "season2, pv", "season3",
                       "season3.error", "season3.pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality")
# Remove congeners with no normal distribution
# Shapiro test p-value < 0.05
lme.pcb$Normality <- as.numeric(lme.pcb$Normality)
# Get the congeners that are not showing normality
lme.pcb.out <- lme.pcb[lme.pcb$Normality < 0.05, ]
lme.pcb <- lme.pcb[lme.pcb$Normality > 0.05, ]

# Export results
write.csv(lme.pcb, file = "Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeLmePCB.csv")

# Generate predictions
# Select congeners that are not showing normality to be remove from che.pcb.2
df <- data.frame(names_to_remove = lme.pcb.out$Congeners)
# Get column indices to remove
cols_to_remove <- which(names(che.pcb.2) %in% df$names_to_remove)
# Remove columns from che.pcb.2 with congeners that don't show normality
che.pcb.3 <- che.pcb.2[, -cols_to_remove]

# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(che.pcb.3[,1]),
                      ncol = length(che.pcb.3[1,]))

for (i in 1:length(che.pcb.3[1,])) {
  fit <- lmer(che.pcb.3[,i] ~ 1 + time + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  lme.fit.pcb[,i] <- fitted(fit)
}

# Estimate a factor of 2 between observations and predictions
factor2 <- 10^(che.pcb.3)/10^(lme.fit.pcb)
factor2.pcb <- sum(factor2 > 0.5 & factor2 < 2,
                   na.rm = TRUE)/(sum(!is.na(factor2)))*100

# Individual PCB congener plots -------------------------------------------
# (1)
# Plot 1:1 for all congeners
# Transform lme.fit.pcb to data.frame
lme.fit.pcb <- as.data.frame(lme.fit.pcb)
# Add congener names to lme.fit.pcb columns
colnames(lme.fit.pcb) <- colnames(che.pcb.3)
# Add code number to first column
df1 <- cbind(code = row.names(che.pcb.3), che.pcb.3)
df2 <- cbind(code = row.names(lme.fit.pcb), lme.fit.pcb)

for (i in 2:length(df1)) {
  col_name <- if (i == 1) {
    ""  # leave the name empty for the first plot
  } else {
    names(df1)[i] # use the column name for other plots
  }
  
  # create plot for each pair of columns
  p <- ggplot(data = data.frame(x = df1$code, y1 = 10^(df1[, i]), y2 = 10^(df2[, i])),
              aes(x = y1, y = y2)) +
    geom_point(shape = 21, size = 3, fill = "white") +
    scale_y_log10(limits = c(0.5, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(limits = c(0.5, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
    ylab(expression(bold("Predicted lme concentration PCBi (pg/L)"))) +
    theme_bw() +
    theme(aspect.ratio = 15/15) +
    annotation_logticks(sides = "bl") +
    geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
    geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
    geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
    annotate('text', x = 10^1, y = 10^4, label = gsub("\\.", "+", names(df1)[i]),
             size = 3, fontface = 2)
  # save plot
  ggsave(paste0("Output/Plots/Sites/ObsPred/ChesapeakeBay/", col_name, ".png"),
         plot = p, width = 6, height = 6, dpi = 500)
}

# (2)
# All plots in one page
# Create a list to store all the plots
plot_list <- list()

# Loop over the columns of df1 and df2
for (i in 2:length(df1)) {
  col_name <- paste(names(df1)[i], sep = "")  # use the column name for plot title
  # Create plot for each pair of columns and add to plot_list
  p <- ggplot(data = data.frame(x = df1$code, y1 = 10^(df1[, i]), y2 = 10^(df2[, i])),
              aes(x = y1, y = y2)) +
    geom_point(shape = 21, size = 3, fill = "white") +
    scale_y_log10(limits = c(0.5, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(limits = c(0.5, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
    ylab(expression(bold("Predicted lme concentration PCBi (pg/L)"))) +
    theme_bw() +
    theme(aspect.ratio = 15/15, 
          axis.title = element_text(size = 8)) +
    annotation_logticks(sides = "bl") +
    annotate('text', x = 25, y = 10^4, label = gsub("\\.", "+", col_name),
             size = 2.5, fontface = 2) +
    geom_abline(intercept = 0, slope = 1, col = "white", linewidth = 0.7) +
    geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
    geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7)
  
  plot_list[[i-1]] <- p  # add plot to list
}
# Combine all the plots using patchwork
combined_plot <- wrap_plots(plotlist = plot_list, ncol = 4)
# Save the combined plot
ggsave("Output/Plots/Sites/ObsPred/ChesapeakeBay/combined_plot.png", combined_plot,
       width = 15, height = 15, dpi = 500)

# (3)
# Create a list to store all the cleaned data frames
cleaned_df_list <- list()
# Loop over the columns of df1 and df2
for (i in 2:length(df1)) {
  # Create a new data frame by binding the columns of df1 and df2 for each pair of columns
  df_pair <- cbind(df1[,1], df1[,i], df2[,i])
  colnames(df_pair) <- c("code", "observed", "predicted")
  # Remove the rows with missing values
  cleaned_df_pair <- na.omit(df_pair)
  # Add the cleaned data frame to the list
  cleaned_df_list[[i-1]] <- cleaned_df_pair
}

{
  # Modify data to be plotted
  # Combine all the cleaned data frames using rbind
  combined_cleaned_df <- do.call(rbind, cleaned_df_list)
  # Convert the matrix to a data frame
  combined_cleaned_df <- as.data.frame(combined_cleaned_df)
  # Convert the code column to a factor
  combined_cleaned_df$code <- as.factor(combined_cleaned_df$code)
  # Convert the observed and predicted columns to numeric
  combined_cleaned_df[,2:3] <- apply(combined_cleaned_df[,2:3], 2, as.numeric)
}

# Plot all the pairs together
p <- ggplot(combined_cleaned_df, aes(x = 10^(observed), y = 10^(predicted))) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(0.1, 10^5), 
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.1, 10^5), 
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration PCBi (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 15/15, 
        axis.title = element_text(size = 10)) +
  annotation_logticks(sides = "bl") +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  annotate("text", x = 5, y = 10^4.5,
           label = expression(atop("Chesapeake Bay",
                                   paste("21 PCB congeners (n = 1090 pairs)"))),
           size = 4, fontface = 2)
# See plot
print(p)
# Save plot
ggsave("Output/Plots/Sites/ObsPred/ChesapeakeBay/ChesapeakeBayObsPredPCB.png",
       plot = p, width = 8, height = 8, dpi = 500)

         