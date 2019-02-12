# Explore soils data

# load packages ####
library(tidyverse)
library(ggplot2)
library(ggtern)
library(plyr)

# load phyloseq object (Full WWP data set) ####
wwp = readRDS("./output/WWP_Full_phyloseq_object.RDS")
wwp.pa = readRDS("./output/WWP_Full_phyloseq_object_presence-absence.RDS")

bact = readRDS("./output/WWP_16S_phyloseq.RDS")
fung = readRDS("./output/WWP_ITS_phyloseq.RDS")

# remove any empty samples or ASVs
wwp.pa = subset_taxa(wwp.pa, sample_sums(wwp.pa) > 0)
bact = subset_taxa(bact, sample_sums(bact) > 0)
fung = subset_taxa(fung, sample_sums(fung) > 0)

wwp.pa = subset_samples(wwp.pa, taxa_sums(wwp.pa) > 0)
bact = subset_samples(bact, taxa_sums(bact) > 0)
fung = subset_samples(fung, taxa_sums(fung) > 0)


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




