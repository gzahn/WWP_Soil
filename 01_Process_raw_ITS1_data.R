# Load packages ####
library(dada2); packageVersion("dada2")
library(purrr)
library(tidyr)
library(ggplot2)
library(readxl)
library(phyloseq)
library(decontam)


# Find raw fastq files and prepare workspace ####
path <- "./ITS/fwd" 
list.files(path)

# Parse fwd and rev reads
fnFs <- sort(list.files(path, pattern=".ITS1.fastq.gz", full.names = TRUE))

# Get Sample Names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Peek at quality profiles
plotQualityProfile(fnFs[c(80)]) # fwd reads

# Make filtered outfile names
filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))

# Filter and trim ####
# cut fwd reads at 290 and rev reads at 160
out <- filterAndTrim(fnFs, filtFs,
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

# reassign filts for lost samples
filtFs <- list.files("./ITS/fwd/filtered", pattern = "filt.fastq.gz", full.names = TRUE)

# Learn error rates ####
errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e8)

# plot error rates for sanity
plotErrors(errF, nominalQ=TRUE)

# Dereplicate ####
derepFs <- derepFastq(filtFs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sapply(strsplit(basename(names(derepFs)), "_"), `[`, 1)

# Sample Inferrence ####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

# Construct sequence table ####
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

# Inspect distribution of sequence lengths and remove short (<100bp) reads
table(nchar(getSequences(seqtab)))
longseqs = which(nchar(getSequences(seqtab)) >= 100)
seqtab = seqtab[,longseqs]

# Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Save progress
saveRDS(seqtab.nochim,"./ITS/output/seqtab.nochim.RDS")

# Track reads through pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out)
colnames(track) <- c("input", "filtered")
# rownames(track) <- sample.names

# Export results
write.csv(track, file = "./ITS/output/tracked_reads.csv", quote = FALSE)

seqtab.nochim <- readRDS("./ITS/output/seqtab.nochim.RDS")
# Assign taxonomy - Custom database including LOTS of outgroups ####
taxa <- assignTaxonomy(seqtab.nochim,"./taxonomy/custom_RDP_ITS1_database.fasta.gz", multithread = TRUE)

# save it
saveRDS(taxa,"./ITS/output/taxa.RDS")

# rename seqtab object samples
seqtab.df <- as.data.frame(seqtab.nochim)
rownames(seqtab.df)
#++++++++++++++++++++++++++++++++++

# Evaluate accuracy against mock community ####
unqs.mock <- seqtab.df["WS-Pos-1-F",]
mockhits = which(unqs.mock > 1)
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
mock.genus = as.character(as.data.frame(taxa)[mockhits,"Genus"])
unqs.mock = unqs.mock[which(mock.genus != "NA")]
mock.genus = mock.genus[which(mock.genus != "NA")]
mock.results = data.frame(Genus = mock.genus, Abund = unqs.mock)

ggplot(mock.results, aes(x=Genus,y=Abund/sum(mock.results$Abund))) +
  geom_bar(stat = "identity")
ggsave("./ITS/output/mock_community_observed.png")

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

# subset meta to just remaining ITS samples
in.meta = which(names(seqtab.nochim[,1]) %in% meta$SampleID == TRUE)
seqtab.fwd.nochim = seqtab.nochim[in.meta,]
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
saveRDS(ps, "./ITS/output/phyloseq_object_ITS.RDS")

ps@sam_data$controls

# identify and remove contaminants ####
contamdf.prev <- isContaminant(ps, neg=c(which(ps@sam_data$controls==TRUE)), threshold = 0.1)
table(contamdf.prev$contaminant)

png(filename = "./ITS/output/neg_control_asv_prevalence.png")
plot(colSums(otu_table(ps)[,!meta$controls]), main = "ASV Prevalence with Extraction Negatives")
points(colSums(otu_table(ps)[,meta$controls]), col = "Red")
dev.off()

#remove them
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
identical(ps, ps.noncontam)

# save contam-free phyloseq object
saveRDS(ps.noncontam, "./ITS/output/phyloseq_object_ITS_noncontam.RDS")
