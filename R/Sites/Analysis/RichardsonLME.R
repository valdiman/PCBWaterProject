## Water PCB concentrations data analysis per site
## Richardson Hill Road Landfill
## Arolclor method

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
  rhr.tpcb <- cbind(factor(rhr$SiteID), rhr$SampleDate,
                    rhr$Latitude, rhr$Longitude, as.matrix(rhr$tPCB),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(rhr.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- rhr.tpcb$tPCB
time <- rhr.tpcb$time
site <- rhr.tpcb$site.code
season <- rhr.tpcb$season
# tPCB vs. time + season + flow + temp + site
lme.rhr.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.rhr.tpcb)
# Look at residuals
{
  res.rhr.tpcb <- resid(lme.rhr.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/RichardsonQ-QtPCB.pdf")
  qqnorm(res.rhr.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.rhr.tpcb)
  dev.off()
}
# Shapiro test
shapiro.test(resid(lme.rhr.tpcb)) # p-value = 0.4767

# Create matrix to store results
{
  lme.tpcb <- matrix(nrow = 1, ncol = 16)
  lme.tpcb[1] <- fixef(lme.rhr.tpcb)[1] # intercept
  lme.tpcb[2] <- summary(lme.rhr.tpcb)$coef[1,"Std. Error"] # intercept error
  lme.tpcb[3] <- summary(lme.rhr.tpcb)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb[4] <- fixef(lme.rhr.tpcb)[2] # time
  lme.tpcb[5] <- summary(lme.rhr.tpcb)$coef[2,"Std. Error"] # time error
  lme.tpcb[6] <- summary(lme.rhr.tpcb)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb[7] <- fixef(lme.rhr.tpcb)[3] # season 1
  lme.tpcb[8] <- summary(lme.rhr.tpcb)$coef[3,"Std. Error"] # season 3 error
  lme.tpcb[9] <- summary(lme.rhr.tpcb)$coef[3,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb[10] <- -log(2)/lme.tpcb[4]/365 # t0.5
  lme.tpcb[11] <- abs(-log(2)/lme.tpcb[4]/365)*lme.tpcb[5]/abs(lme.tpcb[4]) # t0.5 error
  lme.tpcb[12] <- as.data.frame(VarCorr(lme.rhr.tpcb))[1,'sdcor']
  lme.tpcb[13] <- as.data.frame(r.squaredGLMM(lme.rhr.tpcb))[1, 'R2m']
  lme.tpcb[14] <- as.data.frame(r.squaredGLMM(lme.rhr.tpcb))[1, 'R2c']
  lme.tpcb[15] <- shapiro.test(resid(lme.rhr.tpcb))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(lme.rhr.tpcb)
  # Calculate residuals and RMSE
  residuals <- log10(tpcb) - predictions
  non_na_indices <- !is.na(residuals)
  lme.tpcb[16] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Obtain observations and predictions
# Get predicted values tpcb
fit.lme.values.rhr.tpcb <- as.data.frame(fitted(lme.rhr.tpcb))
# Add column name
colnames(fit.lme.values.rhr.tpcb) <- c("predicted")
# Add predicted values to data.frame
rhr.tpcb$predicted <- 10^(fit.lme.values.rhr.tpcb$predicted)
# Create overall plot prediction vs. observations
predic.obs <- data.frame(tPCB = rhr.tpcb$tPCB, predicted = rhr.tpcb$predicted)
predic.obs <- data.frame(Location = rhr$LocationName[1], predic.obs)
colnames(predic.obs) <- c("location", "observed", "predicted")
# Save new data
write.csv(predic.obs,
          "Output/Data/Sites/csv/Richardson/RichardsonLmeObsPredtPCB.csv",
          row.names = FALSE)

# Estimate a factor of 2 between observations and predictions
rhr.tpcb$factor2 <- rhr.tpcb$tPCB/rhr.tpcb$predicted
factor2.tpcb <- nrow(rhr.tpcb[rhr.tpcb$factor2 > 0.5 & rhr.tpcb$factor2 < 2,
                              ])/length(rhr.tpcb[,1])*100

# Transform lme.tpcb to data.frame so factor 2 can be included
lme.tpcb <- as.data.frame(lme.tpcb)

# Add factor 2 to lme.pcb data.frame
lme.tpcb$factor2 <- factor2.tpcb

# Change number format of factor 2 to 3 significant figures
lme.tpcb$factor2 <- formatC(signif(lme.tpcb$factor2, digits = 3))

# Add column names
colnames(lme.tpcb) <- c("Intercept", "Intercept.error",
                        "Intercept.pv", "time", "time.error", "time.pv",
                        "season3", "season3.error", "season3.pv", "t05",
                        "t05.error", "RandonEffectSiteStdDev", "R2nR", "R2R",
                        "Normality", "RMSE", "Factor2")

# Export results
write.csv(lme.tpcb,
          file = "Output/Data/Sites/csv/Richardson/RichardsonLmetPCB.csv",
          row.names = FALSE)

# Modeling plots
# Plot prediction vs. observations, 1:1 line
p <- ggplot(rhr.tpcb, aes(x = tPCB, y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(10^4, 10^7), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^4, 10^7), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# See plot
print(p)

# Save plot
ggsave("Output/Plots/Sites/ObsPred/Richardson/RichardsonLmeObsPredtPCB.png",
       plot = p, width = 8, height = 8, dpi = 500)

# Plot residuals vs. predictions
{
  # Open a PNG graphics device
  png("Output/Plots/Sites/Residual/res_plotlmeRichardsontPCB.png", width = 800,
      height = 600)
  # Create plot
  plot(rhr.tpcb$predicted, resid(lme.rhr.tpcb),
       points(rhr.tpcb$predicted, resid(lme.rhr.tpcb), pch = 16, 
              col = "white"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  # Add lines to the plot
  abline(0, 0)
  abline(h = c(-1, 1), col = "grey")
  abline(v = seq(0, 1400000, 100000), col = "grey")
  # Close the PNG graphics device
  dev.off()
}

