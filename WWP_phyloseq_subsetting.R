library(phyloseq)

#Load phyloseq objects
Bact = readRDS("./16S/output/phyloseq_object_16S_noncontam.RDS")
Fung = readRDS("./ITS/output/phyloseq_object_ITS_noncontam.RDS")

# Subset to WWP samples
Bact = subset_samples(Bact, Project_Name == "Woolsey_Wet_Prairie")
Fung = subset_samples(Fung, Project_Name == "Woolsey_Wet_Prairie")

# Remove singletons and doublets
Bact = subset_taxa(Bact, taxa_sums(Bact) > 2)
Fung = subset_taxa(Fung, taxa_sums(Fung) > 2)

# Remove empty samples (if any)
Bact = subset_samples(Bact, sample_sums(Bact) > 0)
Fung = subset_samples(Fung, sample_sums(Fung) > 0)

saveRDS(Bact, "./WWP_16S_phyloseq.RDS")
saveRDS(Fung, "./WWP_ITS_phyloseq.RDS")
