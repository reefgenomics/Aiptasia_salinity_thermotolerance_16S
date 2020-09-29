#### adonis and Phyloseq

#setwd("~/Documents/Bioinformatics_scripts/R_scripts/Janna/")
met=read.csv("metadata_2812OTUs.csv")
otu=t(read.table("OTU_table_2812.txt", sep = "\t", header = T, row.names = 1))

library(pairwiseAdonis)
otu=t(read.table("~/Documents/Bioinformatics_scripts/R_scripts/Aiptasia_salinity/Input_files/OTU_table_2812.txt", sep = "\t", header = T, row.names = 1))
groups=ifelse(rownames(otu) %like% "H2", "H2", ifelse(rownames(otu) %like% "CC7", "CC7", "Water")) 
pairwise.adonis(otu, groups , p.adjust.m ='fdr', sim.method = 'bray')

#phyloseq
map= read.table("~/Documents/Bioinformatics_scripts/R_scripts/Aiptasia_salinity/Input_files/metadata_phylo.txt", header = TRUE, row.names = 1) #lookup table

otu.t= otu_table(otu, taxa_are_rows=FALSE)
sam.t= sample_data(map)
#tax.t= tax_table(tax.noConta)
#tre.t=read_tree("R_data/RAxML_bestTree.CBASS84_noConta_ML.tre") UniFracs did not provide good resolution 
#tre.t$tip.label=gsub("/.*", "", tre.t$tip.label)
phy.all= phyloseq(otu.t,  sam.t)
P4=c("#823021","#0F8029", "#5DE1FE")
phy.t=microbiome::transform(phy.all, transform = "compositional", target = "OTU", shift = 0, scale = 1)

PCOA_br = ordinate(phy.t, method = "PCoA", distance = "bray")
pdf("~/Documents/Bioinformatics_scripts/R_scripts/Aiptasia_salinity/outputs/Aiptasia_PCoA_hosts.pdf", width = 5, height = 4, pointsize = 12)
plot_ordination(phy.t, PCOA_br, color = "host", shape = "host")  + geom_point(size = 3, alpha = 1) + theme_bw() + scale_colour_manual(values=P4)  + theme(plot.title = element_text(hjust = 0.5)) + theme_classic() 
dev.off()

