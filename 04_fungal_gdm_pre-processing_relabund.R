# Prepare phyloseq object of fungal ITS1 for GDM analysis ####

# Load required libraries ####
library(phyloseq)


##Load the phyloseq object ####
ps <- readRDS("./output/WWP_ITS_phyloseq_relabund.RDS")
ps2 <- readRDS("./output/WWP_ITS_phyloseq_pres-abs.RDS")
#View the phyloseq object
ps


#Create a reduced sample data file with no categorical variables ####
sample_data(ps) <- sample_data(ps)[,!names(sample_data(ps)) %in% 
                  c("SampleID","cluster","num","topo","Elevation","Longitude","Latitude")]
sample_data(ps2) <- sample_data(ps2)[,!names(sample_data(ps2)) %in% 
                  c("SampleID","cluster","num","topo","Elevation","Longitude","Latitude")]

# Exmine ps otu tables
plot(rowSums(otu_table(ps)))
plot(rowSums(otu_table(ps2)))

#Save reduced sample data phyloseq objects for later use (like GDM) ####
saveRDS(ps,"./output/WWP_ITS_reduced_relabund.RDS")
saveRDS(ps2,"./output/WWP_ITS_reduced_pres-abs.RDS")
