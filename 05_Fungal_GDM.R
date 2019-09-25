# Load packages
library(gdm)
library(tidyverse)
library(Hmisc)
library(corrplot)

# Load phyloseq object (reduced number of variables...from "./04_fungal_gdm_pre-processing_relabund.R")
ps <- readRDS("./output/WWP_ITS_reduced_pres-abs.RDS")

# Set OTU Table and Sample Data
otu <- as.data.frame(as(otu_table(ps),"matrix"))
sd <- as(sample_data(ps),"data.frame")
sd$X <- row.names(sd) # Add column "X" for sample IDs
otu$X <- row.names(sd) # Add column "X" for sample IDs
class(sd);class(otu)


#Full model variables with low ext nutrients removed
vars <- c("X","lon","lat","pH","EC","GWC","per_clay","per_sand","per_silt",
          "Al","Ca","Co","Cr","Cu","Fe","K","Mg","Mn","Na","Ni","P","S",
          "Ti","Zn","per_som","PerN","PerC","CNRatio","Ele")
sd.full <- sd[,vars]


#change the predData to the subset you want
gdmTab <- formatsitepair(bioData = otu, bioFormat = 1, XColumn = "lon", YColumn = "lat", 
                         siteColumn = "X", predData = sd.full, abundance = FALSE,
                         dist = "jaccard")

gdmTab[1:3,]


# Run a model to test for var significance and contribution to model.
# Significance testing of environmental variables will be done by a 
# combination of Monte Carlo sampling and stepwise backward elimination 
# as implemented in the gdm.varImp function. This is set for 150 permutations 
# per step until only significant variables remain in the model. 

# Try with sampleSitePairs = .01, nPerm = 50
# Re-run on cluster with ncores=12, nPerm = 99, and sampleSitePairs = 1
all.vars <- gdm.varImp(gdmTab, geo=TRUE, splines = NULL, knots = NULL, fullModelOnly = TRUE, 
                       nPerm = 50, sampleSites = 1, 
                       sampleSitePairs = 1, outFile = "GDM_VarImp_Jaccard")

warnings()
# That step takes about forever minutes to run on 4 cores! ... Save the output somewhere handy:
saveRDS(all.vars, "./output/GDM_Var_Imp.RDS")

# inspect GDM variable importance
summary(all.vars)

#Significant varibles only
sd.sig <- sd[,c(1:4,8,10,13,18,32,33)]
names(sd[,c(1:4,8,10,13,18,32,33)])


#change the predData to the subset you want
gdmTab.sig <- formatsitepair(bioData = otu, bioFormat = 1, XColumn = "lon", YColumn = "lat", 
                             siteColumn = "X", predData = sd, abundance = FALSE)
gdmTab.sig[1:3,]

# Run GDM
gdm.sig <- gdm(gdmTab, geo=T)

# Summary
summary(gdm.sig)






# Remove "X" column from OTU table
otu$X <- NULL
