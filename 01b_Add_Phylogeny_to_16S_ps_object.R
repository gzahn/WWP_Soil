library(tidyverse)
library(phyloseq)
library(vegan)
library(phangorn)
library(msa)
library(microbiome)
library(ape)

# Read original ps object ####
ps = readRDS("output/phyloseq_object_16S_noncontam.RDS")

ps@sam_data$Project_Name

# subset to Woolsey Wet Prairie samples and remove singletons #
ps = subset_samples(ps, Project_Name == "Woolsey_Wet_Prairie")
ps = subset_taxa(ps, taxa_sums(ps) > 1)

# Build Phylogenetic Tree ####
seqs <- rownames(tax_table(ps))
names(seqs) <- seqs # This propagates to the tip labels of the tree

# alignment (took 37 hours to complete!)
alignment <- msa(seqs,method = "Muscle", type = "dna",verbose = TRUE,order = "input",maxiters = 12)

# save progress 
saveRDS(alignment,"./output/WWP_dna_alignment_muscle.RDS")


# re-load point
alignment <- readRDS("./output/WWP_dna_alignment_muscle.RDS")
phang.align = as.phyDat(alignment, type = "DNA")

# distance max likelihood
dm <- dist.ml(phang.align)

#save
saveRDS(dm,"./output/ML_Distance.RDS")

# neighbor-joining tree
treeNJ <- NJ(dm) # Note, tip order != sequence order
treeNJ$tip.label <- seqs

#save
saveRDS(treeNJ, "./output/treeNJ.RDS")
treeNJ <- readRDS("./output/treeNJ.RDS")


fit = pml(treeNJ, data=phang.align)

#save
saveRDS(fit,"./output/fit_treeNJ.RDS")

## negative edges length changed to 0!

fitJC <- optim.pml(fit, TRUE)

# save
saveRDS(fitJC, "./output/WWP_fitJC.RDS") # This is the new tree using optim.pml
write.tree(fitJC$tree, file = "./output/WWP_bact_tree_JC.nwk")

# reload point
fitJC <- readRDS("./output/WWP_fitJC.RDS")

(unique(names(seqs)))

# Bootstrap 
bs = bootstrap.pml(fitJC, bs=100, control = pml.control(trace = 0),multicore = TRUE, mc.cores = 10)
?bootstrap.pml
# check tree with plot
plotBS(fitJC$tree,type="p")



#### redo with GTR model?
# fitGTR <- update(fit, k=4, inv=0.2)
# fitGTR <- phangorn::optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                               control = phangorn::pml.control(trace = 1L),rearrangement = "stochastic")


# saveRDS(fitGTR, "./output/WWP_fitGTR2.RDS") # This is the new tree using optim.pml
# fitGTR = readRDS("./output/WWP_fitGTR2.RDS") # loading new tree. does it work??
# write.tree(fitGTR$tree, file = "./output/WWP_bact_tree.nwk")

# bs = bootstrap.pml(fitGTR, bs=100, optNni=TRUE, multicore=TRUE)

detach("package:phangorn", unload=TRUE)

# add tree to phyloseq object ####
ps2 = phyloseq(tax_table(tax_table(ps)),
               otu_table(otu_table(ps)),
               sample_data(sample_data(ps)),
               phy_tree(fitJC$tree))
saveRDS(ps2, "./output/WWP_phyloseq_object_noncontam_w_ML-JCtree.RDS")

