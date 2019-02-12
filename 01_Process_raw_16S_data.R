# Process Raw 16S reads ####

# Load packages ####
library(dada2); packageVersion("dada2")
library(purrr)
library(tidyr)
library(ggplot2)
library(readxl)
library(decontam)
library(phyloseq)

# Find raw fastq files and prepare workspace ####
path <- "./16S" 
list.files(path)

# Parse fwd and rev reads
fnFs <- sort(list.files(path, pattern="_R1_", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_", full.names = TRUE))

# Get Sample Names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Peek at quality profiles
plotQualityProfile(fnFs[c(1,2,190,200)]) # fwd reads
plotQualityProfile(fnRs[c(1,2,190,200)]) # rev reads

# Make filtered outfile names
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Filter and trim ####
# cut fwd reads at 290 and rev reads at 160
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(290,160),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

# reassign filts for lost samples
filtFs <- list.files("./16S/filtered", pattern = "_F_filt", full.names = TRUE)
filtRs <- list.files("./16S/filtered", pattern = "_R_filt", full.names = TRUE)



# Learn error rates ####
errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e8)
errR <- learnErrors(filtRs, multithread=TRUE, nbases = 1e8)

# plot error rates for sanity
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplicate ####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepRs) <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)

# Sample Inferrence ####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge fwd and rev reads ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, trimOverhang = TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Save progress
saveRDS(seqtab.nochim,"./16S/output/seqtab.nochim.RDS")

## Track reads through pipeline ####
# getN <- function(x) sum(getUniques(x))
# track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
# rownames(track) <- sample.names
# 
# # Export results
# write.csv(track, file = "./16S/output/tracked_reads.csv", quote = FALSE)
# long.track = gather(as.data.frame(track), Step, Reads, 1:6)
# 
# # looks like merging is where a lot of reads are lost. Poor quality rev reads!  Just go with fwd reads
# ggplot(long.track, aes(x=Step,y=Reads)) + geom_violin() +
# theme_bw()  
# ggsave("./16S/output/tracked_reads.png")

# Generate new seqtable from just fwd reads ####
seqtab.fwd <- makeSequenceTable(dadaFs)
seqtab.fwd.nochim <- removeBimeraDenovo(seqtab.fwd, method="consensus", multithread=TRUE, verbose=TRUE)

# Save Data object
saveRDS(seqtab.fwd.nochim, "./16S/output/seqtab.fwd.nochim.RDS")

# Re-load, if necessary
seqtab.fwd.nochim <- readRDS("./16S/output/seqtab.fwd.nochim.RDS")

# remove all ASVs that don't have at least 10 hits
seqtab.fwd.nochim <- seqtab.fwd.nochim[,(colSums(seqtab.fwd.nochim) > 9)]

# Assign taxonomy - Silva v132 exact match / 80% bootstrap min ####
taxa <- assignTaxonomy(seqtab.fwd.nochim,"./taxonomy/silva_nr_v132_train_set.fa.gz", minBoot = 80,multithread = TRUE)
saveRDS(taxa,"./16S/output/taxa.RDS")
taxa <- addSpecies(taxa, "./taxonomy/silva_species_assignment_v132.fa.gz")
saveRDS(taxa,"./16S/output/taxa_exact.RDS")

# rename seqtab object samples
seqtab.df <- as.data.frame(seqtab.fwd.nochim)
row.names(seqtab.df) <- map(strsplit(row.names(seqtab.df), "_"),1)

# Evaluate accuracy against mock community ####
unqs.mock <- seqtab.df["WS-Pos-1-B",]
mockhits = which(unqs.mock > 1)
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
mock.genus = as.character(as.data.frame(taxa)[mockhits,"Genus"])
unqs.mock = unqs.mock[which(mock.genus != "NA")]
mock.genus = mock.genus[which(mock.genus != "NA")]
mock.results = data.frame(Genus = mock.genus, Abund = unqs.mock)

ggplot(mock.results, aes(x=Genus,y=Abund/sum(mock.results$Abund))) +
  geom_bar(stat = "identity")
ggsave("./16S/output/mock_community_observed.png")

cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.seqs = names(seqtab.df[,mockhits])
mock.ref <- getSequences(file.path("./taxonomy/mock/mock.fasta"))
match.ref <- sum(sapply(mock.seqs, function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")





# Create PhyloSeq object ####

# read in metadata
meta = read_xlsx("./metadata.xlsx")
meta$controls <- meta$SampleSource
meta$controls[meta$controls == "Extraction Negative"] <- TRUE
meta$controls[meta$controls != TRUE] <- FALSE
meta$controls = as.logical(meta$controls)
  
# subset meta to just remaining 16S samples
in.meta = which(names(seqtab.fwd.nochim[,1]) %in% meta$SampleID == TRUE)
seqtab.fwd.nochim = seqtab.fwd.nochim[in.meta,]
dim(seqtab.fwd.nochim)
in.seqtab = which(meta$SampleID %in% names(seqtab.fwd.nochim[,1]))
meta = meta[in.seqtab,]

#re-order
meta = meta[order(meta$SampleID),]
row.names(meta) <- meta$SampleID

# df = as.data.frame(seqtab.fwd.nochim, row.names = 1)
# row.names(df) <- meta$SampleID
identical(row.names(seqtab.fwd.nochim), meta$SampleID)
dim(taxa)
dim(seqtab.fwd.nochim)
dim(meta)

# make phyloseq object 
ps <- phyloseq(otu_table(seqtab.fwd.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa))
# save it
saveRDS(ps, "./16S/output/phyloseq_object_16S.RDS")

ps@sam_data$controls

# identify and remove contaminants ####
contamdf.prev <- isContaminant(ps, neg=c(325:327), threshold = 0.1)
table(contamdf.prev$contaminant)

png(filename = "./16S/output/neg_control_asv_prevalence.png")
plot(colSums(otu_table(ps)[,!meta$controls]), main = "ASV Prevalence with Extraction Negatives")
points(colSums(otu_table(ps)[,meta$controls]), col = "Red")
dev.off()

#remove them
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
identical(ps, ps.noncontam)

# save contam-free phyloseq object
saveRDS(ps.noncontam, "./16S/output/phyloseq_object_16S_noncontam.RDS")
