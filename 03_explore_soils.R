# Explore soils data

# load packages ####
library(tidyverse)
library(ggplot2)
library(ggtern)
library(plyr)
library(vegan)
library(ade4)

# load phyloseq object (Full WWP data set) ####
wwp = readRDS("./output/WWP_Full_phyloseq_object.RDS")
wwp.pa = readRDS("./output/WWP_Full_phyloseq_object_presence-absence.RDS")

# bact = readRDS("./output/WWP_16S_phyloseq.RDS")
# fung = readRDS("./output/WWP_ITS_phyloseq.RDS")

# remove any empty samples or ASVs
wwp.pa = subset_taxa(wwp.pa, sample_sums(wwp.pa) > 0)
bact = subset_taxa(wwp, sample_sums(bact) > 0)
fung = subset_taxa(wwp, sample_sums(fung) > 0)

wwp.pa = subset_samples(wwp.pa, taxa_sums(wwp.pa) > 0)
bact = subset_samples(wwp, taxa_sums(bact) > 0)
fung = subset_samples(wwp, taxa_sums(fung) > 0)


# Examine metadata ####
glimpse(as.data.frame(sample_data(wwp.pa)))
glimpse(as.data.frame(sample_data(fung)))

# Species richness (bact + fungi)
richness.total = sample_sums(wwp.pa)
richness.bact = sample_sums(transform_sample_counts(bact, function(abund) 1*(abund>0)))
richness.fung = sample_sums(transform_sample_counts(fung, function(abund) 1*(abund>0)))

# Sand-Silt-Clay Plot ####

# make sand-silt-clay df
ssc = data.frame(Sand = sample_data(wwp.pa)$per_sand, Silt = sample_data(wwp.pa)$per_silt, Clay = sample_data(wwp.pa)$per_clay)
ssc.bact = data.frame(Sand = sample_data(bact)$per_sand, Silt = sample_data(bact)$per_silt, Clay = sample_data(bact)$per_clay)
ssc.fung = data.frame(Sand = sample_data(fung)$per_sand, Silt = sample_data(fung)$per_silt, Clay = sample_data(fung)$per_clay)


# Make base plot of soil textures chart
data(USDA)

# Put tile labels at the midpoint of each tile.
USDA.LAB = ddply(USDA, 'Label', function(df) {
  apply(df[, 1:3], 2, mean)
})

# Tweak
USDA.LAB$Angle = 0
USDA.LAB$Angle[which(USDA.LAB$Label == 'Loamy Sand')] = -35

base = ggplot(data = USDA, aes(y=Clay, x=Sand, z=Silt)) +
  coord_tern(L="x",T="y",R="z") +
  geom_polygon(alpha = 0.25, size = 0.5, color = 'black',aes(color=Label,fill=Label)) +
  geom_text(data = USDA.LAB,
            aes(label = Label, angle = Angle),
            color = 'black',
            size = 3.5) +
  theme_rgbw() +
  theme_showsecondary() +
  theme_showarrows() +
  custom_percent("Percent") +
  theme(legend.justification = c(0, 1),
        legend.position      = c(0, 1)) +
  labs(title = 'USDA Textural Classification Chart',
       fill  = 'Textural Class',
       color = 'Textural Class')

base + geom_point(data=ssc,alpha = .25, aes(size=richness.total)) + labs(size = "Total Richness")
ggsave("./output/figs/Soil_Textures_w_Total_Richness.png", dpi=300,width = 12,height = 12)

base + geom_point(data=ssc.bact,alpha = .25, aes(size=richness.bact)) + labs(size = "Bacterial Richness")
ggsave("./output/figs/Soil_Textures_w_Bacteria_Richness.png", dpi=300,width = 12,height = 12)

base + geom_point(data=ssc.fung,alpha = .25, aes(size=richness.fung)) + labs(size = "Fungal Richness")
ggsave("./output/figs/Soil_Textures_w_Fungal_Richness.png", dpi=300,width = 12,height = 12)


# explore soil metadata distributions ####

# convert metadata to "vanilla" data frame
meta <- as(sample_data(wwp.pa), "data.frame")
element.names=c("Al","As","B","Ca","Cd","Co","Cr","Cu","Fe","K","Mg","Mn","Mo","Na","Ni","P","S","Ti","Zn")

# add variable showing whether each sample is "mound" or "depression"

??????

# gather element values to make plotting easier
elements = gather(meta,key = Element, value = Value, element.names)

ggplot(elements, aes(x=Value, fill = Element)) +
  geom_histogram() + facet_wrap(~Element, scales = "free") + theme_bw()
ggsave("./output/figs/Soil_Element_Distributions.png", dpi=300, width = 12, height = 12)

ggplot(elements, aes(x=Value, fill = factor(cluster))) +
  geom_histogram(alpha=.5) + facet_wrap(~Element, scales = "free") + theme_bw()
ggsave("./output/figs/Soil_Element_Distributions_w_cluster.png", dpi=300, width = 12, height = 12)


# Mantel tests
# Mantel Test ####
spatial.dist.full = dist(cbind(wwp.pa@sam_data$lon, wwp.pa@sam_data$lat), method = "euclidean")
comm.dist.full = dist(as.matrix(wwp.pa@otu_table), method = "binary")
mantel.test.full = mantel.rtest(spatial.dist.full, comm.dist.full, nrepet = 999)

ggplot(mapping = aes(x=jitter(spatial.dist,amount=1), y=comm.dist.full)) +
  geom_point(alpha=.05) + stat_smooth(method = "lm") + 
  labs(x="Spatial Distance",y="Full Community Distance") + theme_bw()
ggsave("./output/figs/full_community_mantel.png", dpi=300)

spatial.dist.bact = dist(cbind(bact@sam_data$Longitude, bact@sam_data$Latitude), method = "euclidean")
comm.dist.bact = dist(as.matrix(bact@otu_table), method = "binary")
mantel.test.bact = mantel.rtest(spatial.dist.bact, comm.dist.bact, nrepet = 999)






# WCMD Ordination ####
# Try WCMD
otus = dist(otu_table(wwp.pa))
wcmd.full = wcmdscale(otus,k=2, eig=TRUE)
wcmd.df = data.frame(Dim1 = wcmd.full$points[,1], Dim2 = wcmd.full$points[,2], cluster = meta$cluster)

ggplot(wcmd.df, aes(x=Dim1,y=Dim2, color = factor(cluster))) +
  geom_point()



biplot(wcmd.full)