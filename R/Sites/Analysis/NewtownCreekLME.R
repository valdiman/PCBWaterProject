## Water PCB concentrations data analysis per site
## Newtown Creek
## Data do not work for any of the models.

# Install packages
install.packages("tidyverse")
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
install.packages("scales")

# Load libraries
{
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
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# Select Newtown Creek data ---------------------------------------------------
ntc <- wdc[str_detect(wdc$LocationName, 'Newtown Creek'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  ntc$SampleDate <- as.Date(ntc$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(ntc$SampleDate) - min(as.Date(ntc$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(ntc$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  ntc.tpcb <- cbind(factor(ntc$SiteID), as.matrix(ntc$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(ntc.tpcb) <- c("SiteID", "tPCB", "time", "season")
}

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- ntc.tpcb$tPCB
time <- ntc.tpcb$time
site <- ntc.tpcb$SiteID
season <- ntc.tpcb$season

# tPCB vs. time + season + flow + temp + site
lme.ntc.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# Something not working!
# See results
summary(lme.ntc.tpcb)
# Look at residuals
{
  res.ntc.tpcb <- resid(lme.ntc.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/TributariesQ-QtPCB.pdf")
  qqnorm(res.ntc.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.ntc.tpcb)
  dev.off()
}
# Shapiro test
shapiro.test(resid(lme.ntc.tpcb)) # p-value = 0.2819

# Create matrix to store results from lme analysis
{
  lme.tpcb <- matrix(nrow = 1, ncol = 16)
  lme.tpcb[1] <- fixef(lme.ntc.tpcb)[1] # intercept
  lme.tpcb[2] <- summary(lme.ntc.tpcb)$coef[1,"Std. Error"] # intercept error
  lme.tpcb[3] <- summary(lme.ntc.tpcb)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb[4] <- fixef(lme.ntc.tpcb)[2] # time
  lme.tpcb[5] <- summary(lme.ntc.tpcb)$coef[2,"Std. Error"] # time error
  lme.tpcb[6] <- summary(lme.ntc.tpcb)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb[7] <- fixef(lme.ntc.tpcb)[3] # season 3
  lme.tpcb[8] <- summary(lme.ntc.tpcb)$coef[3,"Std. Error"] # season 3 error
  lme.tpcb[9] <- summary(lme.ntc.tpcb)$coef[3,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb[10] <- -log(2)/lme.tpcb[4]/365 # t0.5
  lme.tpcb[11] <- abs(-log(2)/lme.tpcb[4]/365)*lme.tpcb[5]/abs(lme.tpcb[4]) # t0.5 error
  lme.tpcb[12] <- as.data.frame(VarCorr(lme.ntc.tpcb))[1,'sdcor']
  lme.tpcb[13] <- as.data.frame(r.squaredGLMM(lme.ntc.tpcb))[1, 'R2m']
  lme.tpcb[14] <- as.data.frame(r.squaredGLMM(lme.ntc.tpcb))[1, 'R2c']
  lme.tpcb[15] <- shapiro.test(resid(lme.ntc.tpcb))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(lme.ntc.tpcb)
  # Calculate residuals and RMSE
  residuals <- log10(tpcb) - predictions
  non_na_indices <- !is.na(residuals)
  lme.tpcb[16] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Obtain observations and predictions
# Get predicted values tpcb
fit.lme.values.ntc.tpcb <- as.data.frame(fitted(lme.ntc.tpcb))
# Add column name
colnames(fit.lme.values.ntc.tpcb) <- c("predicted")
# Add predicted values to data.frame
ntc.tpcb$predicted <- 10^(fit.lme.values.ntc.tpcb$predicted)
# Something odd with the model. Predictions = Observations!

# Estimate a factor of 2 between observations and predictions
ntc.tpcb$factor2 <- ntc.tpcb$tPCB/ntc.tpcb$predicted
factor2.tpcb <- nrow(ntc.tpcb[ntc.tpcb$factor2 > 0.5 & ntc.tpcb$factor2 < 2,
                              ])/length(ntc.tpcb[,1])*100

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

# Add Location Name
lme.tpcb <- cbind(LocationName = rep("NewTown Creek",
                                     nrow(lme.tpcb)), lme.tpcb)
# Select relevant columns
lme.tpcb.t <- lme.tpcb[, c("LocationName", "t05", "t05.error",
                           "R2R", "RMSE", "Factor2")]
# Time coefficient not significant

# LME for individual PCBs -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  ntc.pcb <- subset(ntc, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  ntc.pcb <- subset(ntc.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  ntc.pcb <- log10(ntc.pcb)
  # Replace -inf to NA
  ntc.pcb <- do.call(data.frame,
                     lapply(ntc.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Check if there are any NA values in the dataframe
  if (anyNA(ntc.pcb)) {
    # Remove columns with more than 70% NA values
    ntc.pcb.1 <- ntc.pcb[, -which(colSums(is.na(ntc.pcb))/nrow(ntc.pcb) > 0.7)]
  } else {
    # If there are no NA values, you might want to handle it differently,
    # such as keeping all columns
    ntc.pcb.1 <- ntc.pcb
  }
  # Add site ID
  SiteID <- factor(ntc$SiteID)
  # Change date format
  SampleDate <- as.Date(ntc$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(ntc$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to ntc.pcb.1
  ntc.pcb.1 <- cbind(ntc.pcb.1, SiteID, data.frame(time.day), season.s)
  # Remove metadata
  ntc.pcb.2 <- subset(ntc.pcb.1, select = -c(SiteID:season.s))
}

# Get covariates
time <- ntc.pcb.1$time
season <- ntc.pcb.1$season
site <- ntc.pcb.1$SiteID

# Create matrix to store results
lme.pcb <- matrix(nrow = length(ntc.pcb.2[1,]), ncol = 16)

# Perform LME
for (i in 1:length(ntc.pcb.2[1,])) {
  fit <- lmer(ntc.pcb.2[,i] ~ 1 + time + season + (1|site),
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
  lme.pcb[i,7] <- fixef(fit)[3] # # season 3
  lme.pcb[i,8] <- summary(fit)$coef[3,"Std. Error"] # season 3 error
  lme.pcb[i,9] <- summary(fit)$coef[3,"Pr(>|t|)"] # # season 3 p-value
  lme.pcb[i,10] <- -log(2)/lme.pcb[i,4]/365 # t0.5
  lme.pcb[i,11] <- abs(-log(2)/lme.pcb[i,4]/365)*lme.pcb[i,5]/abs(lme.pcb[i,4]) # t0.5 error
  lme.pcb[i,12] <- as.data.frame(VarCorr(fit))[1,'sdcor']
  lme.pcb[i,13] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
  lme.pcb[i,14] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
  lme.pcb[i,15] <- shapiro.test(resid(fit))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(fit)
  # Calculate residuals and RMSE
  residuals <- ntc.pcb.2[, i] - predictions
  non_na_indices <- !is.na(residuals)
  lme.pcb[i, 16] <- sqrt(mean(residuals[non_na_indices]^2))
}

# There is something with the data. Stop analysis.
