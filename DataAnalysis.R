## Water PCB concentrations analysis.
# Data were obtained from EPA and contractors from PCB Superfund
# sites in USA

# Install packages
install.packages("ggpubr")
install.packages("ggpmisc")
install.packages("tidyverse")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("robustbase")
install.packages("dplyr")
install.packages("tibble")
install.packages("Matrix")
install.packages("lme4")
install.packages("MuMIn")
install.packages("lmerTest")

# Load libraries
library(ggplot2)
#library(ggpubr)
#library(ggpmisc)
library(scales) # function trans_breaks
#library(gridExtra)
#library(tidyverse)
#library(reshape2)
library(stringr) # str_detect
library(robustbase) #function colMedians
library(dplyr) # %>%
library(tibble) # function add a column
library(lme4) # perform lme
library(MuMIn) # get Rs from lme
library(lmerTest) # get the p-value from lme

# Data in pg/L
wdc.0 <- read.csv("WaterDataCongenerAroclor08052022.csv")

# Summary plots -----------------------------------------------------------

# Box plots with Aroclor dataset
# Remove samples (rows) with total PCBs  = 0
wdc.1 <- wdc.0[!(rowSums(wdc.0[, c(12:115)], na.rm = TRUE)==0),]
# Remove metadata
wdc.1 <- subset(wdc.1, select = -c(ID:AroclorCongener))
# Remove Aroclor data
wdc.1 <- subset(wdc.1, select = -c(A1016:A1260))

# Create a frequency detection plot
wdc.freq <- colSums(! is.na(wdc.1) & (wdc.1 !=0))/nrow(wdc.1)
wdc.freq <- data.frame(wdc.freq)
colnames(wdc.freq) <- c("PCB.frequency")
congener <- row.names(wdc.freq)
wdc.freq <- cbind(congener, wdc.freq$PCB.frequency)
colnames(wdc.freq) <- c("congener", "PCB.frequency")
wdc.freq <- data.frame(wdc.freq)
wdc.freq$congener <- as.character(wdc.freq$congener)
wdc.freq$congener <- gsub('\\.', '+', wdc.freq$congener) # replace dot for +
wdc.freq$PCB.frequency <- as.numeric(as.character(wdc.freq$PCB.frequency))
wdc.freq$congener <- factor(wdc.freq$congener,
                            levels = rev(wdc.freq$congener)) # change the order to be plotted.

# Summary statistic of frequency of detection
summary(wdc.freq$PCB.frequency)

# Frequency detection plot
ggplot(wdc.freq, aes(x = 100*PCB.frequency, y = congener)) +
  geom_bar(stat = "identity", fill = "#66ccff") +
  ylab("") +
  theme_bw() +
  xlim(c(0,100)) +
  theme(aspect.ratio = 20/5) +
  xlab(expression(bold("Frequency detection (%)"))) +
  theme(axis.text.x = element_text(face = "bold", size = 9),
        axis.title.x = element_text(face = "bold", size = 10)) +
  theme(axis.text.y = element_text(face = "bold", size = 5))

# Summary statistic of total PCBs in pg/L
summary(rowSums(wdc.1, na.rm = T))

# Total PCBs in 1 box plot
ggplot(wdc.1, aes(x = "", y = rowSums(wdc.1, na.rm = T))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  theme(aspect.ratio = 14/2) +
  xlab(expression(bold(Sigma*"PCB (n = 5798)")))+
  ylab(expression(bold("Water Concentration 1990 - 2020 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10,
                                    angle = 45, hjust = 1.8,
                                    vjust = 2)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotation_logticks(sides = "l")

# Individual congeners
# Summary statistic of individual congeners in pg/L
summary(wdc.1, na.rm = T, zero = T)
# Obtain the median for each individual congener
cong.median <- as.numeric(sub('.*:',
                              '', summary(wdc.1, na.rm = T,
                                          zero = T)[3,]))

# Individual PCB boxplot
ggplot(stack(wdc.1), aes(x = ind, y = values)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot(width = 0.6, outlier.colour = "#66ccff", col = "#66ccff",
               outlier.shape = 1) +
  scale_x_discrete(labels = wdc.freq$congener) + # use to change the "." to "+"
  theme_bw() +
  theme(aspect.ratio = 25/135) +
  xlab(expression("")) +
  ylab(expression(bold("PCB congener concentration (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8,
                                   color = "black"),
        axis.title.y = element_text(face = "bold", size = 8,
                                    color = "black")) +
  theme(axis.text.x = element_text(face = "bold", size = 6,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.6, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# Spatial plots and analysis ----------------------------------------------
# Modify x-axis
# States
sites <- c("CA", "DE", "ID", "IN", "MA", "MI", "MO",
           "MT", "NM", "NY", "OR", "TX", "WA", "WI")

# Total PCBs
ggplot(wdc.0, aes(x = factor(StateSampled, levels = sites),
                y = rowSums(wdc.0[, c(12:115)],  na.rm = T))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2019 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  geom_hline(yintercept = 63, color = "#cc0000") # median from line 78

# Selected StateSampled and individual PCB congener
wdc.pcb.sp <- subset(wdc.0, select = c(StateSampled, PCB4.10))

# Plot
ggplot(wdc.pcb.sp, aes(x = factor(StateSampled, levels = sites),
                   y = PCB4.10)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Water Conncetration PCB 4+10 1990 - 2019 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  geom_hline(yintercept = 0.03, color = "#cc0000") # median 

# Temporal analysis and plots---------------------------------------------
# Analysis
# Needs to add a individual number for each site name
# to perform regression analysis, linear and Linear Mixed-Effects Model (LME)
site.numb <- wdc.0$SiteName %>% as.factor() %>% as.numeric
wdc.site <- add_column(wdc.0,
                          site.numb, .after = "AroclorCongener")
wdc.site$SampleDate <- as.Date(wdc.site$SampleDate,
                                  format = "%m/%d/%y")
# Calculate total PCB per sample
tpcb.tmp <- rowSums(wdc.site[, c(13:116)], na.rm = T)
time.day <- data.frame(as.Date(wdc.site$SampleDate) - min(as.Date(wdc.site$SampleDate)))
# Generate data.frame for analysis and plots
tpcb.tmp <- cbind(data.frame(time.day), as.matrix(tpcb.tmp),
                  wdc.site$SampleDate)
colnames(tpcb.tmp) <- c("time", "tPCB", "date")

# Perform linear regression (lr)
# (i) + 1
lr.tpcb <- lm(log10(tPCB + 1) ~ time, data = tpcb.tmp)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (ii) tPCB^(1/4)
lr.tpcb.2 <- lm(1/(tPCB + 1)^0.5 ~ time, data = tpcb.tmp)
# See results
summary(lr.tpcb.2)
# Look at residuals
res.2 <- resid(lr.tpcb.2) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.2)
# Add a straight diagonal line to the plot
qqline(res.2)
# Shapiro test
shapiro.test(res.2)
# One-sample Kolmogorov-Smirnov test
ks.test(res.2, 'pnorm')

# Extract intercept and slope (std error and p-value too) values
inter.lr <- summary(lr.tpcb)$coef[1,"Estimate"]
slope.lr <- summary(lr.tpcb)$coef[2,"Estimate"]
slope.lr.error <- summary(lr.tpcb)$coef[2,"Std. Error"]
slope.lr.pvalue <- summary(lr.tpcb)$coef[2,"Pr(>|t|)"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5.lr <- -log(2)/slope.lr/365 # half-life tPCB in yr = -log(2)/slope/365
# Calculate error
t0.5.lr.error <- abs(t0.5.lr)*slope.lr.error/abs(slope.lr)
# Extract R2 adjust
R2.adj <- summary(lr.tpcb)$adj.r.squared

# Perform Linear Mixed-Effects Model (LME)
time <- tpcb.tmp$time
site <- wdc.site$site.numb

# (i) + 1
lmem.tpcb <- lmer(log10(tpcb.tmp$tPCB/4 + 1) ~ 1 + time + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.tpcb)
# Look at residuals
res <- resid(lmem.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res, main = "log10(C + 1)")
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (ii) tPCB^(1/4)
lmem.tpcb.2 <- lmer((tpcb.tmp$tPCB)^(1/4) ~ 1 + time + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.tpcb.2)
# Look at residuals
res.2 <- resid(lmem.tpcb.2) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.2, main = "C^(1/4)")
# Add a straight diagonal line to the plot
qqline(res.2)
# Shapiro test
shapiro.test(res.2)
# One-sample Kolmogorov-Smirnov test
ks.test(res.2, 'pnorm')

# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.tpcb))[1, 'R2c']
# Extract intercept and slope (std error and p-value too) values
inter.lmem <- summary(lmem.tpcb)$coef[1, 1]
slope.lmem <- summary(lmem.tpcb)$coef[2, 1]
slope.lmem.error <- summary(lmem.tpcb)$coef[2, 2]
slope.lmem.pvalue <- summary(lmem.tpcb)$coef[ ,5]
# Std. Dev. of random effect
as.data.frame(VarCorr(lmem.tpcb))[1,'sdcor']
# Half-life in yr
t0.5.lmem <- -log(2)/slope.lmem/365
t0.5.lmem.error <- abs(t0.5.lmem)*slope.lmem.error/abs(slope.lmem)

# Final plot tPCB, including lr and lmem
ggplot(tpcb.tmp, aes(y = tPCB,
                       x = format(date,'%Y'))) +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(Sigma*"PCB concentration (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 7,
                                   angle = 45, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 10)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)
  geom_abline(intercept = inter.lr,
              slope = slope.lr*365, color = "#fc8d62", size = 0.8) +
  geom_abline(intercept = inter.lmem,
              slope = slope.lmem*365, color = "#7570b3", size = 0.8)

# Congeners
# Format data for selected PCBs
# Linear regression per congener
pcbi.tmp.lr <- subset(d.cong.site, select = -c(ID:site.numb))
pcbi.tmp.lr <- subset(pcbi.tmp.lr, select = -c(A1016:A1260))
# log10
log.pcbi.tmp.lr <- log10(pcbi.tmp.lr)
# Remove infinite values
log.pcbi.tmp.lr[log.pcbi.tmp.lr == -Inf] = NA

# Create matrix to store data
lr.pcbi <- matrix(nrow = length(log.pcbi.tmp.lr), ncol = 7)

for(i in 1:length(log.pcbi.tmp.lr)) {
  fit.lr <- lm(log.pcbi.tmp.lr[,i] ~ tpcb.tmp$time)
  lr.pcbi[i,1] <- summary(fit.lr)$coef[1,"Estimate"] # intercept
  lr.pcbi[i,2] <- summary(fit.lr)$coef[2,"Estimate"] # slope
  lr.pcbi[i,3] <- summary(fit.lr)$coef[2,"Std. Error"] # slope error
  lr.pcbi[i,4] <- summary(fit.lr)$coef[2,"Pr(>|t|)"] # slope p-value
  lr.pcbi[i,5] <- -log(2)/lr.pcbi[i,2]/365 # t0.5
  lr.pcbi[i,6] <- abs(-log(2)/lr.pcbi[i,2]/365)*lr.pcbi[i,3]/abs(lr.pcbi[i,2]) # t0.5 error
  lr.pcbi[i,7] <- summary(fit.lr)$adj.r.squared # R2 adj
}

# Add column names
colnames(lr.pcbi) <- c("intercept", "slope", "slope.error",
                       "p-value", "half-life", "half-life.error", "R2.adj")
# Just 3 significant figures
lr.pcbi <- formatC(signif(lr.pcbi, digits = 3))
# Add PCB congener names
congener <- colnames(pcbi.tmp.lr)
lr.pcbi2 <- cbind(congener, lr.pcbi)
# Export results
write.csv(lr.pcbi, file = "linearRegressionPCBi.csv")

# Linear Mixed-effects model per individual congeners
# Create matrix to store data
lmem.pcbi <- matrix(nrow = length(log.pcbi.tmp.lr), ncol = 8)

for(i in 1:length(log.pcbi.tmp.lr)) {
  fit.lmem <- lmer(log.pcbi.tmp.lr[,i] ~ 1 + time + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"))
  lmem.pcbi[i,1] <- summary(fit.lmem)$coef[1,1] # intercept
  lmem.pcbi[i,2] <- summary(fit.lmem)$coef[2,1] # slope
  lmem.pcbi[i,3] <- summary(fit.lmem)$coef[2,"Std. Error"] # slope error
  lmem.pcbi[i,4] <- summary(fit.lmem)$coef[2,"Pr(>|t|)"] # slope p-value
  lmem.pcbi[i,5] <- -log(2)/lmem.pcbi[i,2]/365 # t0.5
  lmem.pcbi[i,6] <- abs(-log(2)/lmem.pcbi[i,2]/365)*lmem.pcbi[i,3]/abs(lmem.pcbi[i,2]) # t0.5 error
  lmem.pcbi[i,7] <- as.data.frame(r.squaredGLMM(fit.lmem))[1, 'R2m'] # R2 w/o random effect
  lmem.pcbi[i,8] <- as.data.frame(r.squaredGLMM(fit.lmem))[1, 'R2c'] # R2 w/ random effect
}

# Add column names
colnames(lmem.pcbi) <- c("intercept", "slope", "slope.error",
                       "p-value", "half-life", "half-life.error", "R2.nre", "R2.re")
# Just 3 significant figures
lmem.pcbi <- formatC(signif(lmem.pcbi, digits = 3))
# Add PCB congener names
congener <- colnames(pcbi.tmp.lr)
lmem.pcbi <- cbind(congener, lmem.pcbi)
# Export results
write.csv(lmem.pcbi, file = "LmemPCBi.csv")

# Individual PCB plots
# Include sampling dates
plot.pcbi.tmp <- cbind(tpcb.tmp$date, pcbi.tmp.lr)
# change date name
colnames(plot.pcbi.tmp)[colnames(plot.pcbi.tmp) == "tpcb.tmp$date"] <- "date"

# Plot
ggplot(plot.pcbi.tmp, aes(y = PCB1,
                          x = format(date,'%Y'))) +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("PCBi concentration (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 7,
                                   angle = 45, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 10)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotation_logticks(sides = "l") +
  geom_smooth(method = lm, se = TRUE, formula = y ~ x,
              aes(group = 1), colour = "#ff6611", size = 0.5)


# Plot and analysis per sites ---------------------------------------------
# (i) Fox River, WI
# Select Fox River only
Fox.River <- d.cong.0[str_detect(d.cong.0$SiteName, 'FoxRiver'),] 
# Remove samples (rows) with total PCBs  = 0
Fox.River <- Fox.River[!(rowSums(Fox.River[, c(12:115)], na.rm = TRUE)==0),]
# Remove metadata
Fox.River.1 <- subset(Fox.River, select = -c(ID:AroclorCongener))
# Remove Aroclor data
Fox.River.1 <- subset(Fox.River.1, select = -c(A1016:A1260))

# Summary statistic of total PCBs
summary(rowSums(Fox.River.1, na.rm = T))

# Total PCBs in 1 box plot
ggplot(Fox.River.1, aes(x = "", y = rowSums(Fox.River.1, na.rm = T))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)(c(0.1, 1e4)),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  theme(aspect.ratio = 14/2) +
  xlab(expression(bold(Sigma*"PCB (n = 151)")))+
  ylab(expression(bold("Fox River 2010 - 2018 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10,
                                    angle = 45, hjust = 1.8,
                                    vjust = 2)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotation_logticks(sides = "l")

# Individual congeners
# Summary statistic of individual congeners
summary(Fox.River.1, na.rm = T, zero = T)
# Obtain the median for each individual congener
Fox.River.cong.median <- as.numeric(sub('.*:',
                                        '', summary(Fox.River.1,
                                                    na.rm = T, zero = T)[3,]))

# Individual PCB boxplot
ggplot(stack(Fox.River.1), aes(x = ind, y = values)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot(width = 0.6, outlier.colour = "#66ccff", col = "#66ccff",
               outlier.shape = 1) +
  scale_x_discrete(labels = d.cong.freq$congener) + # This is from the frequency detection plot
  theme_bw() +
  theme(aspect.ratio = 25/135) +
  xlab(expression("")) +
  ylab(expression(bold("PCB congener (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8,
                                   color = "black"),
        axis.title.y = element_text(face = "bold", size = 8,
                                    color = "black")) +
  theme(axis.text.x = element_text(face = "bold", size = 6,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.6, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# Spatial plots and analysis
# Modify x-axis
sites.FR <- c("LakeWinnebago", "OperableUnit1", "OperableUnit2A",
              "OperableUnit2B", "OperableUnit2C", "OperableUnit3")

# Total PCBs
ggplot(Fox.River, aes(x = factor(SiteSampled, levels = sites.FR),
                   y = rowSums(Fox.River[, c(12:115)],  na.rm = T))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)(c(0.1, 1e6)),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(Sigma*"PCB 2010 - 2018 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# Selected PCB congeners
# Format data for selected PCBs
# Remove samples with = 0 and NA
FR.pcbi.sp <- subset(Fox.River,
                        Fox.River$PCB4.10 != 0 & Fox.River$PCB4.10 != "NA")

# Plot
ggplot(FR.pcbi.sp, aes(x = factor(SiteSampled, levels = sites.FR),
                          y = PCB4.10)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)(c(0.1, 1e6)),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("PCB 4+10 2010 - 2018 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# Temporal plot
# Total PCBs
ggplot(Fox.River, aes(y = rowSums(Fox.River[, c(12:115)],  na.rm = T),
                   x = format(SampleDate,'%Y'))) +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(Sigma*"PCB concentration (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 7,
                                   angle = 45, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 10)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotation_logticks(sides = "l")

# (ii) Hudson River, NY
# Select Kalamazoo River only
Hud.River <- d.cong[str_detect(d.cong$SiteName, 'HudsonRiver'),] 
# Remove samples (rows) with total PCBs  = 0
Hud.River <- Hud.River[!(rowSums(Hud.River[, c(12:115)], na.rm = TRUE)==0),]
# Remove metadata
Hud.River.1 <- subset(Hud.River, select = -c(ID:AroclorCongener))
# Remove Aroclor data
Hud.River.1 <- subset(Hud.River.1, select = -c(A1016:A1260))

# Summary statistic of total PCBs
summary(rowSums(Hud.River.1, na.rm = T))

# Total PCBs in 1 box plot
ggplot(Hud.River.1, aes(x = "", y = rowSums(Hud.River.1, na.rm = T))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  theme(aspect.ratio = 14/2) +
  xlab(expression(bold(Sigma*"PCB (n = 453)")))+
  ylab(expression(bold("Hudson River 2005 - 2017 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10,
                                    angle = 45, hjust = 1.8,
                                    vjust = 2)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotation_logticks(sides = "l")

# Individual congeners
# Summary statistic of individual congeners
summary(Hud.River.1, na.rm = T, zero = T)
# Obtain the median for each individual congener
Hud.River.cong.median <- as.numeric(sub('.*:', '', summary(Hud.River.1,
                                                           na.rm = T, zero = T)[3,]))

# Individual PCB boxplot
ggplot(stack(Hud.River.1), aes(x = ind, y = values)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot(width = 0.6, outlier.colour = "#66ccff", col = "#66ccff",
               outlier.shape = 1) +
  scale_x_discrete(labels = d.cong.freq$congener) + # This is from the frequency detection plot
  theme_bw() +
  theme(aspect.ratio = 25/135) +
  xlab(expression("")) +
  ylab(expression(bold("PCB congener (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8,
                                   color = "black"),
        axis.title.y = element_text(face = "bold", size = 8,
                                    color = "black")) +
  theme(axis.text.x = element_text(face = "bold", size = 6,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.6, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# Spatial plots and analysis
# Modify x-axis
sites.HR <- c("Albany", "BakersFalls", "DixBridge", "Lock1",
              "MohawkRiverAtCohoes", "Poughkeepsie", "RogersIsland",
              "Schuylerville", "Stillwater", "ThompsonIsland",
              "ThompsonIslandDam", "Waterford")

# Total PCBs
ggplot(Hud.River, aes(x = factor(SiteSampled, levels = sites.HR),
                      y = rowSums(Hud.River[, c(12:115)],  na.rm = T))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(Sigma*"PCB 2005 - 2017 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# Selected PCB congeners
# Format data for selected PCBs
# Remove samples with = 0 and NA
HR.pcbi.sp <- subset(Hud.River,
                     Hud.River$PCB4.10 != 0 & Hud.River$PCB4.10 != "NA")

# Plot
ggplot(HR.pcbi.sp, aes(x = factor(SiteSampled, levels = sites.HR),
                       y = PCB4.10)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)(c(0.1, 1e4)),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("PCB 4+10 2005 - 2017 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# Temporal plot
# Total PCBs
ggplot(Hud.River, aes(y = rowSums(Hud.River[, c(12:115)],  na.rm = T),
                      x = format(SampleDate,'%Y'))) +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(Sigma*"PCB concentration (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 7,
                                   angle = 45, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 10)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotation_logticks(sides = "l")


