# Code to combine all Passaic River data

# Read data ---------------------------------------------------------------
# Data in pg/L
pass.1 <- read.csv("Data/PassaicRiver/2010-2016_Diamond-Alkali_OU3_EPA-NBSA-Split-samplesFV.csv")
pass.2 <- read.csv("Data/PassaicRiver/2011 CPG CWCM Sampling  - Round 1FV.csv")
pass.3 <- read.csv("Data/PassaicRiver/2012 CPG CWCM Sampling - Low FlowFV.csv")
pass.4 <- read.csv("Data/PassaicRiver/2012 CPG CWCM Sampling - Round 2FV.csv")
pass.5 <- read.csv("Data/PassaicRiver/2012 CPG CWCM Sampling - Round 4FV.csv")
pass.6 <- read.csv("Data/PassaicRiver/2012 CPG CWCM Sampling - Round 5FV.csv")
pass.7 <- read.csv("Data/PassaicRiver/2013 CPG CWCM Sampling - High Flow 1FV.csv")
pass.8 <- read.csv("Data/PassaicRiver/2013 CPG CWCM Sampling - High Flow 2FV.csv")
pass.9 <- read.csv("Data/PassaicRiver/2017-2019_OU2_PDI_Water_Column_20210924FV.csv")

merged_pass <- rbind(pass.1, pass.2, pass.3, pass.4, pass.5, pass.6, pass.7,
                   pass.8, pass.9)
