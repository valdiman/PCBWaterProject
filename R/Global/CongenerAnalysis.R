## Water PCB concentrations analysis.
## Data were obtained from EPA and contractors from PCB Superfund
## sites in USA. Data only AroclorCongener = Congener

# Install packages
install.packages("ggplot2")
install.packages('FactoMineR')
install.packages('factoextra')

# Load libraries
{
  library(ggplot2)
  library(FactoMineR) # Perform PCA
  library(factoextra) # Plot result from PCA
}

# Read data
wdc <- read.csv("Data/WaterDataCongenerAroclor08052022.csv")

# Data preparation
cong <- subset(wdc, AroclorCongener == "Congener")
cong <- cong[!(rowSums(cong[, 14:117], na.rm = TRUE) == 0), ]
cong.1 <- subset(cong, !LocationName == "New Bedford Harbor") # Only 18 congeners
cong.1 <- subset(cong, select = -c(SampleID:AroclorCongener))
cong.1 <- subset(cong.1, select = -c(A1016:A1260))

# Create an average PCB profile distribution
tmp <- rowSums(cong.1, na.rm = TRUE)
prof <- sweep(cong.1, 1, tmp, FUN = "/")
prof.ave <- colMeans(prof, na.rm = TRUE)
prof.sd <- apply(prof, 2, sd, na.rm = TRUE)
congener <- as.character(colnames(prof))
prof.ave <- data.frame(congener, mean = as.numeric(prof.ave), sd = as.numeric(prof.sd))
prof.ave$congener <- factor(prof.ave$congener, levels = unique(congener))

# Plot profiles -----------------------------------------------------------
# Plot average PCB profile
ggplot(prof.ave, aes(x = congener, y = mean)) +
  geom_bar(position = position_dodge(), stat = "identity",
           fill = "black") +
  geom_errorbar(aes(ymin = mean, ymax = (mean+sd)), width = 0.9,
                position = position_dodge(0.9)) +
  xlab("") +
  ylim(0, 0.5) +
  theme_bw() +
  theme(aspect.ratio = 4/12) +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  annotate("text", x = 5.4, y = 0.4, label = "PCBs 4+10", size = 2,
           fontface = 1, angle = 90) +
  annotate("text", x = 8, y = 0.36, label = "PCB 11", size = 2,
           fontface = 1, angle = 90) +
  annotate("text", x = 14, y = 0.32, label = "PCBs 18+30", size = 2,
           fontface = 1, angle = 90) +
  annotate("text", x = 18, y = 0.33, label = "PCBs 20+21+28+\n31+33+50+53",
           size = 2, fontface = 1, angle = 90) +
  annotate("text", x = 30.5, y = 0.25, label = "PCBs 43+49+\n52+69+73",
           size = 2, fontface = 1, angle = 90) +
  annotate("text", x = 39.3, y = 0.30, label = "PCBs 61+66+70+74+\n76+93+95+98+100+102",
           size = 2, fontface = 1, angle = 90) +
  annotate("text", x = 70.5, y = 0.33, label = "PCBs 132+153+\n161+168",
           size = 2, fontface = 1, angle = 90)

# PCB congener analysis ---------------------------------------------------
# Prepare data for PCA
# Add sample names to first column
prof <- cbind(cong$SampleID, prof)
# (1) All samples
# Subset with samples with more than 75% congeners
prof.1 <- prof[rowMeans(!is.na(prof)) >= 0.70, ]
# Remove congeners with < 75% detection frequency
prof.2 <- prof.1[, colMeans(!is.na(prof.1)) >= 0.75]
# Perform PCA
PCA.1 <- PCA(prof.2[, -1], graph = FALSE)
fviz_eig(PCA.1, addlabels = TRUE, ylim = c(0, 100))
fviz_pca_var(PCA.1, col.var = "cos2",
             repel = TRUE) 
fviz_pca_ind(PCA.1, geom.ind = "point", pointshape = 21, 
             pointsize = 2, col.ind = "black", palette = "jco", 
             addEllipses = TRUE, label = "var",
             col.var = "black", repel = TRUE)

# (2) Samples with Method 1668
# Data preparation
cong.1668 <- subset(wdc, EPAMethod == "M1668")
cong.1668 <- cong.1668[!(rowSums(cong.1668[, c(14:117)], na.rm = TRUE)==0), ]
sampleID <- cong.1668$SampleID
cong.1668 <- subset(cong.1668, select = -c(SampleID:AroclorCongener))
cong.1668 <- subset(cong.1668, select = -c(A1016:A1260))
tmp <- rowSums(cong.1668, na.rm = TRUE)
prof.1668 <- sweep(cong.1668, 1, tmp, FUN = "/")
# Add sample names to first column
prof.1668 <- cbind(sampleID, prof.1668)
# Subset data.frame to rows with >= 75% non-NA values
prof.1668.1 <- prof.1668[rowMeans(!is.na(prof.1668)) >= 0.6, ]
# Remove congeners with < 75% detection frequency
prof.1668.2 <- prof.1668[, colMeans(!is.na(prof.1668.1)) >= 0.75]
# Perform PCA
PCA.2 <- PCA(prof.1668.2[,-1], graph = FALSE)
fviz_eig(PCA.2, addlabels = TRUE, ylim = c(0, 100))
fviz_pca_ind(PCA.2, geom.ind = "point", pointshape = 21, 
             pointsize = 2, col.ind = "black", palette = "jco", 
             addEllipses = TRUE, label = "var",
             col.var = "black", repel = TRUE)

# Cosine theta analysis ---------------------------------------------------
# Samples with 100% congeners only
prof.cos.1 <- prof[rowMeans(!is.na(prof)) >= 0.9, ]
# Transpose and remove sample names
prof.cos.2 <- t(prof[,-1])
# Create matrix to storage results
costheta <- matrix(nrow = length(prof.cos.2[1,]),
                   ncol = length(prof.cos.2[1,]))
# Perform Cosine Theta
for (i in 1:length(prof.cos.2[1,])) {
  for (j in 1:length(prof.cos.2[1,])) {
    m1 <- prof.cos.2[,i]
    m2 <- prof.cos.2[,j]
    costheta[i,j] <- sum(m1*m2)/(sum(m1^2)*sum(m2^2))^0.5
  }
}
# Just 3 significant figures
costheta <- formatC(signif(costheta, digits = 3))
# Remove upper diagonal values
costheta[upper.tri(costheta)] <- NA
# Add name to columns
colnames(costheta) <- prof.cos.1[,1]
# Add names to rows
rownames(costheta) <- prof.cos.1[,1]
# Export data
write.csv(costheta, file = "Output/Data/csv/costheta.csv")

