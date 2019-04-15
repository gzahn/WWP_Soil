# Explore soils data

# load packages ####
library(tidyverse)
library(ggplot2)
library(ggtern)
library(plyr)
library(vegan)
library(ade4)
library(phyloseq)
library(colorblindr)
library(indicspecies)
library(cooccur)

# Custom palette
pal = c("#c4a113","#c1593c","#643d91","#38302f","#477887","#688e52","#12aa91","#705f36","#8997b2","#33312c","#3c3e44","#b3bf2d")
palette_plot(pal)
colorblindr::cvd_grid(palette_plot(pal))


# load and clean phyloseq object (Full WWP data set) ####
wwp = readRDS("./output/WWP_Full_phyloseq_object.RDS")

# remove non-fungal and non-bacterial reads
pa = subset_taxa(wwp,Kingdom != "Eukaryota")
pa = subset_taxa(pa,Kingdom != "k__Viridiplantae")
pa = subset_taxa(pa,Kingdom != "k__Metazoa")
pa = subset_taxa(pa,Kingdom != "k__Rhizaria")
pa = subset_taxa(pa,Kingdom != "Eukaryota")
pa = subset_taxa(pa,Kingdom != "NA")
wwp <- pa
rm(pa)

# Convert to presence-absence
wwp.pa <- transform_sample_counts(wwp, function(abund) 1*(abund>0))

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
  theme_rgbw() +
  theme_showsecondary() +
  theme_showarrows() +
  custom_percent("Percent") +
  theme(legend.justification = c(0, 1),
        legend.position      = c(0, 1)) +
  labs(title = 'USDA Textural Classification Chart',
       fill  = 'Textural Class',
       color = 'Textural Class')

base + geom_point(data=ssc,alpha = .25, aes(size=richness.total)) + labs(size = "Total Richness") +
  scale_fill_manual(values = pal) + 
  geom_text(data = USDA.LAB,
            aes(label = Label, angle = Angle),
            color = 'black',
            size = 3.5, face = "bold")
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
sample_data(wwp.pa)$Position = sample_data(wwp.pa)$Ele >= 375.2
sample_data(wwp.pa)$Position[sample_data(wwp.pa)$Position == TRUE] <- "Mound"
sample_data(wwp.pa)$Position[sample_data(wwp.pa)$Position == FALSE] <- "InterMound"

sample_data(fung)$Position = sample_data(fung)$Ele >= 375.2
sample_data(fung)$Position[sample_data(fung)$Position == TRUE] <- "Mound"
sample_data(fung)$Position[sample_data(fung)$Position == FALSE] <- "InterMound"

sample_data(bact)$Position = sample_data(bact)$Ele >= 375.2
sample_data(bact)$Position[sample_data(bact)$Position == TRUE] <- "Mound"
sample_data(bact)$Position[sample_data(bact)$Position == FALSE] <- "InterMound"


# gather element values to make plotting easier
elements = gather(meta,key = Element, value = Value, element.names)
elements$Position = sample_data(wwp.pa)$Position


ggplot(elements, aes(x=Value, fill = Element)) +
  geom_histogram() + facet_wrap(~Element, scales = "free") + theme_bw()
ggsave("./output/figs/Soil_Element_Distributions.png", dpi=300, width = 12, height = 12)

ggplot(elements, aes(x=Value, fill = Position)) +
  geom_histogram(alpha=1) + facet_wrap(~Element, scales = "free") + theme_bw() +
  scale_fill_manual(values = c(pal[1],pal[5]))
ggsave("./output/figs/Soil_Element_Distributions_w_position.png", dpi=300, width = 12, height = 12)


# Boron looks bimodal. What's the reason? Check other soil metadata for possible correlations
B = filter(elements, Element == "B")
for(i in names(B)[c(2:19)]){
  plot(B[,"Value"] ~ B[,i], xlab=i,ylab="Boron")
} # Don't see any correlations between B and other metadata...perhaps microbe community can explain?

# Add categorical B variable
B.clust = kmeans(x = sample_data(wwp.pa)$B,2)$cluster
B.clust[B.clust == 2] <- "High"
B.clust[B.clust == 1] <- "Low"

sample_data(wwp.pa)$B.cat <- B.clust

  
  
  
# Mantel Test ####
spatial.dist.full = dist(cbind(wwp.pa@sam_data$lon, wwp.pa@sam_data$lat), method = "euclidean")
comm.dist.full = dist(as.matrix(wwp.pa@otu_table), method = "binary")
mantel.test.full = mantel.rtest(spatial.dist.full, comm.dist.full, nrepet = 999)

ggplot(mapping = aes(x=jitter(spatial.dist.full,amount=1), y=comm.dist.full)) +
  geom_point(alpha=.05) + stat_smooth(method = "lm") + 
  labs(x="Spatial Distance",y="Full Community Distance") + theme_bw()
ggsave("./output/figs/full_community_mantel.png", dpi=300)

spatial.dist.bact = dist(cbind(bact@sam_data$Longitude, bact@sam_data$Latitude), method = "euclidean")
comm.dist.bact = dist(as.matrix(bact@otu_table), method = "binary")
spatial.dist.fung = dist(cbind(fung@sam_data$Longitude, fung@sam_data$Latitude), method = "euclidean")
comm.dist.fungi = dist(as.matrix(fung@otu_table), method = "binary")

plot(comm.dist.fungi)
mantel.test.bact = mantel.rtest(spatial.dist.full, comm.dist.bact, nrepet = 999)
mantel.test.fung = mantel.rtest(spatial.dist.fung, comm.dist.fungi, nrepet = 999)


# Compare community distance with soil composition distance
metals = c("Al","As","B","Ca","Cd","Co","Cr","Cu","Fe","K","Mg","Mn","Mo","Na","Ni","P","S","Ti","Zn")
metals.dist = dist(meta[,metals])
mantel.test.metals = mantel.rtest(metals.dist,comm.dist.full,nrepet = 999)
plot(mantel.test.metals)

ggplot(mapping = aes(x=metals.dist, y=comm.dist.full)) +
  geom_point(alpha=.05) + stat_smooth(method = "lm") + 
  labs(x="Elemental Composistion Distance",y="Full Community Distance") + theme_bw() +
  ylim(c(.5,1.01))

ggplot(mapping = aes(x=metals.dist, y=comm.dist.bact)) +
  geom_point(alpha=.05) + stat_smooth(method = "lm") + 
  labs(x="Elemental Composistion Distance",y="Bacterial Community Distance") + theme_bw() +
  ylim(c(.5,1.01))

metals.dist.fung = dist(sample_data(fung)[,metals])
ggplot(mapping = aes(x=metals.dist.fung, y=comm.dist.fungi)) +
  geom_point(alpha=.05) + stat_smooth(method = "lm") + 
  labs(x="Elemental Composistion Distance",y="Fungal Community Distance") + theme_bw() +
  ylim(c(.5,1.01))


# Heatmap of soil property relative values ####
metals.table = meta[,metals]
soil.vars = c(metals,"pH","EC","GWC","per_clay","per_sand","per_silt","per_som","PerN","PerC","CNRatio","Ele")
soil.vars.table = meta[,soil.vars]
soil.vars.norm = as.matrix(decostand(soil.vars.table,"normalize",2))

# create color vectors for grouping
mound.colors = plyr::mapvalues(sample_data(wwp.pa)$Position, from=c("InterMound","Mound"), to = c(pal[1],pal[5]))
# color by bact vs fung dominated?

fungal.pa = otu_table(subset_taxa(wwp.pa, Kingdom == "k__Fungi"))
bacterial.pa = otu_table(subset_taxa(wwp.pa, Kingdom == "Bacteria"))
fungal.presence = rowSums(fungal.pa)/ntaxa(fungal.pa)
bacterial.presence = rowSums(bacterial.pa)/ntaxa(bacterial.pa)
plot(fungal.presence, bacterial.presence)
predominance = (fungal.presence - bacterial.presence)
predominance[predominance > 0] <- "Fungal"
predominance[predominance != "Fungal"] <- "Bacterial"
group.colors = plyr::mapvalues(predominance, from = c(unique(predominance)), to= c(pal[3],pal[6]))

# heatmaps 
# colored by mound/inter-mound
heatmap(as.matrix(metals.table), col = gray.colors(100), Colv = NA, RowSideColors = mound.colors)
# colored by fungal/bacterial dominance
heatmap(as.matrix(metals.table), col = gray.colors(100), Colv = NA, RowSideColors = group.colors)


# WCMD Ordination ####
# Try WCMD
otus = dist(otu_table(wwp.pa))
wcmd.full = wcmdscale(otus,k=2, eig=TRUE)
wcmd.df = data.frame(Dim1 = wcmd.full$points[,1], Dim2 = wcmd.full$points[,2], 
                     Position = sample_data(wwp.pa)$Position, 
                     Cluster = sample_data(wwp.pa)$cluster)

ggplot(wcmd.df, aes(x=Dim1,y=Dim2, color = Position)) +
  geom_point() + stat_ellipse() + theme_bw()
ggsave("./output/figs/NMDS_Mound-Intermound.png",dpi=300)



# remove empty sample from wwp.pa
empty = names(which(sample_sums(otu_table(wwp.pa))==0))
wwp.pa2 = subset_samples(wwp.pa,sample_names(wwp.pa) != empty)


# PermANOVA ####
mod.adonis = adonis(otu_table(wwp.pa2) ~ sample_data(wwp.pa2)$cluster + sample_data(wwp.pa2)$PerN + 
                      sample_data(wwp.pa2)$PerC + sample_data(wwp.pa2)$Position,method = "jaccard")
sink(file = "./output/PermANOVA_Table.txt")
mod.adonis
sink(NULL)

# indicator species ####

# meta2 = as(sample_data(wwp.pa),"data.frame")
# otu2 = as(otu_table(wwp.pa),"matrix")
# otu2 = as.data.frame(otu2)
# 
# 
# full.indval = multipatt(otu2,cluster = meta2$Position, control = how(nperm=999))
# mound.indicators = summary(full.indval, alpha = .005)
# full.indval.sign = full.indval$sign
# full.indval.sign = full.indval.sign[full.indval.sign$p.value <=0.001,]
# full.indval.sign = full.indval.sign[complete.cases(full.indval.sign),]
# 
# mound.indicators = row.names(full.indval.sign[full.indval.sign$s.Mound == 1,])
# mound.indicators = tax_table(wwp.pa)[mound.indicators]
# table(mound.indicators)

# Species co-occurance ####

# spp matrix

tax = as(tax_table(wwp.pa),"matrix")
taxids = paste(tax[,1],tax[,2],tax[,3],tax[,4],tax[,5],tax[,6],tax[,7])
mat = as(otu_table(pa),"matrix")
mat = as.data.frame(mat)
names(mat) <- taxids

unique(tax_table(pa)[,1])

# remove taxa found only from two or fewer sample
mat <- mat[,(colSums(mat) > 2)]

# make site mask where all species can be present in all samples
sitemask = mat
sitemask[sitemask==0] <- 1
sitemask = (as.matrix(sitemask))
dim(sitemask)
dim(mat)
?cooccur

# commented out because it takes a long time to run... just reload the RDS object on subsequent runs!
# co = cooccur(mat=t(as.matrix(mat)), spp_names = TRUE, site_mask = t(sitemask),thresh = TRUE)
# saveRDS(co,"./output/co-occurrance_obj.RDS")


co = readRDS("./output/co-occurrance_obj.RDS")
ggsave("./output/figs/Co-Occurrance_obs_v_exp.png")

co
summary(co)
plot(co) + theme_void()
ggsave("./output")
print(co)
pair.profile(co)

# find effect sizes ... reload from RDS
# eff = effect.sizes(co)
# saveRDS(eff,"./output/co-occurrances_effect_sizes.RDS")

eff = readRDS("./output/co-occurrances_effect_sizes.RDS")
eff.big = eff[which(abs(eff$effects) > 0.1),]
eff.big = eff.big[!str_detect(eff.big$sp1,"k__Fungi NA NA NA NA NA NA"),]
eff.big = eff.big[!str_detect(eff.big$sp2,"k__Fungi NA NA NA NA NA NA"),]
plot(eff.big$effects)
