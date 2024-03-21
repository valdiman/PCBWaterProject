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
wdc <- read.csv("Data/WaterDataCongenerAroclor09072023.csv")

# Data preparation
cong <- subset(wdc, AroclorCongener == "Congener")
cong <- cong[!(rowSums(cong[, 14:117], na.rm = TRUE) == 0), ]
cong.1 <- subset(cong, !LocationName == "New Bedford Harbor") # Only 18 congeners
cong.1 <- subset(cong, select = -c(SampleID:AroclorCongener))
cong.1 <- subset(cong.1, select = -c(A1016:tPCB))

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
  geom_bar(stat = "identity", fill = "black") +
  geom_errorbar(aes(ymin = mean, ymax = mean + sd), width = 0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw() +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 7,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

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
cong.1668 <- subset(cong.1668, select = -c(A1016:tPCB))
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

