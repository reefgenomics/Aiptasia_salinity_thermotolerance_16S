library(vegan)

#setwd("~/Documents/Bioinformatics_scripts/R_scripts/Janna/")
met=read.csv("metadata_2812OTUs.csv")
otu=t(read.table("OTU_table_2812.txt", sep = "\t", header = T, row.names = 1))
grous=met$ID4[match(rownames(otu), met$ID)]
cha=specpool(otu,rownames(otu))


shann=diversity(otu, index = "shannon", MARGIN = 1, base = exp(1))
sim=diversity(otu, index = "simpson", MARGIN = 1, base = exp(1))
inv.sim=diversity(otu, index = "invsimpson", MARGIN = 1, base = exp(1))
div=as.data.frame(cbind(shann, sim, inv.sim))
spe=as.data.frame(specnumber(otu, MARGIN = 1), row.names = FALSE)
div.1=merge(div,spe, by.x = "row.names", by.y = "row.names")
colnames(div.1) = c("old", "shann","Simpson.D", "Inv.Simpson", "Species")
div.1$Simpson.E=(1/div.1$Inv.Simpson)/div.1$Species ### adding eveness. (1/Inv.simp)/No.species
div.1$ID=map$ID2[match(map$ID,div.1$old)]

div.s=ddply(div.1, .(ID), summarize,sha_mean=mean(shann), sha_std_dev=sd(shann), sim_mean=mean(Simpson.D), sim_std_dev=sd(Simpson.D), invsim_mean=mean(Inv.Simpson), invsim_std_dev=sd(Inv.Simpson),spec_mean=mean(Species), spec_std_dev=sd(Species), even_mean=mean(Simpson.E), even_std_dev=sd(Simpson.E))
rownames(div.s)=div.s$ID

all=merge(cha,div.s,by="row.names",all.x=TRUE)
all.1=all[,c(1:4,12:17,20,21,10)]
write.csv(all.1, file = "Aiptasia_AlphaDiversity2.csv", quote = FALSE, row.names=FALSE)


##### test normality
shapiro.test(div.1$shann) #W = 0.97725, p-value = 0.4167
shapiro.test(div.1$Simpson.D) # W = 0.89145, p-value = 0.0001921
shapiro.test(div.1$Inv.Simpson) # W = 0.9445, p-value = 0.01713
shapiro.test(div.1$Species) # W = 0.97551, p-value = 0.3565
shapiro.test(div.1$Simpson.E) # W = 0.86493, p-value = 2.933e-05

# non-parametric 

kruskal.test(shann ~ ID, data = div.1) # Kruskal-Wallis chi-squared = 11.319, df = 11, p-value = 0.4169
kruskal.test(Simpson.D ~ ID, data = div.1) # Kruskal-Wallis chi-squared = 10.418, df = 11, p-value = 0.4932
kruskal.test(Inv.Simpson ~ ID, data = div.1) # Kruskal-Wallis chi-squared = 10.418, df = 11, p-value = 0.49
kruskal.test(Simpson.E ~ ID, data = div.1) # Kruskal-Wallis chi-squared = 11.384, df = 11, p-value = 0.4116
kruskal.test(Species ~ ID, data = div.1) # Kruskal-Wallis chi-squared = 21.325, df = 11, p-value = 0.03015
