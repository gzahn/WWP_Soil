# load packages
library(phyloseq)
library(tidyverse)

# "not in" operator
'%ni%' = Negate('%in%')

#Load phyloseq objects
Bact = readRDS("./output/phyloseq_object_16S_noncontam.RDS")
Fung = readRDS("./output/phyloseq_object_ITS_noncontam.RDS")

# Subset to WWP samples
Bact = subset_samples(Bact, Project_Name == "Woolsey_Wet_Prairie")
Fung = subset_samples(Fung, Project_Name == "Woolsey_Wet_Prairie")

# remove unwnated columns
bactcols = names(sample_data(Bact)) %in% c("SampleID","Latitude","Longitude","Elevation")
fungcols = names(sample_data(Fung)) %in% c("SampleID","Latitude","Longitude","Elevation")
sample_data(Bact) <- sample_data(Bact)[,bactcols]
sample_data(Fung) <- sample_data(Fung)[,fungcols]

# remove PCR postive control in Fung
Fung = subset_samples(Fung, SampleID != "Fungal-POS")

# Add soil metadata and merge all into one phyloseq object
sample_names(Bact) <- str_remove(sample_names(Bact),"-B")
sample_names(Fung) <- str_remove(sample_names(Fung),"-F")

meta = read.csv("./soils_data_w2_parse.csv")
meta$SampleID <- paste0("WS-",meta$cluster,"-",meta$num)
row.names(meta) <- meta$SampleID

wwp = merge_phyloseq(Fung,Bact)
wwp = merge_phyloseq(wwp,sample_data(meta))

# remove negative controls
wwp = subset_samples(wwp, sample_names(wwp) %ni% c("WS-Neg-1","WS-Neg-2","WS-Neg-3","WS-Pos-1"))

# Correct mistakes in SOM data entry (sample 3.2 = 4.83; sample 15.5 = 2.6725)
sample_data(wwp)["WS-3-2","per_som"] <- 4.83
sample_data(wwp)["WS-15-5","per_som"] <- 2.6725

# Re-separate fungi and bacteria for downstream
Bact = subset_taxa(wwp, Kingdom == "Bacteria")
Fung = subset_taxa(wwp, Kingdom == "k__Fungi")

# Remove singletons and doublets
Bact = subset_taxa(Bact, taxa_sums(Bact) > 2)
Fung = subset_taxa(Fung, taxa_sums(Fung) > 2)
wwp = subset_taxa(wwp, taxa_sums(wwp) > 2)

# Remove empty samples (if any)
Bact = subset_samples(Bact, sample_sums(Bact) > 0)
Fung = subset_samples(Fung, sample_sums(Fung) > 0)
wwp = subset_samples(wwp, sample_sums(wwp) > 0)


# convert wwp OTU table to presence-absence so ITS and 16S can be compared
wwp.pa <- transform_sample_counts(wwp, function(abund) 1*(abund>0))

# Save phyloseq objects as files
saveRDS(Bact, "./output/WWP_16S_phyloseq.RDS")
saveRDS(Fung, "./output/WWP_ITS_phyloseq.RDS")
saveRDS(wwp.pa, "./output/WWP_Full_phyloseq_object_presence-absence.RDS")
saveRDS(wwp,"output/WWP_Full_phyloseq_object.RDS")


