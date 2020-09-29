library("phyloseq")
library("ggplot2")
library("gridExtra")
library(reshape)
library(tidyr)
library(dplyr)
library(scales)

setwd("~/Documents/Projects/Janna/")

otu = t(read.table("Janna.final.OTU_table", header = FALSE))
### Using column 2 as the title of the columns (Aiptasia_001, Aiptasia_002....) ###
colnames(otu) = otu[2, ]
### removing the first 3 rows from the table now ##
otu=otu[-c(1,2,3), ] 


## Reads the taxonomy csv file, first line is considered as header, each tab is seperated, so in this case separated into OTU, Size, Taxonomy ####
tax = read.csv("Janna.final.taxonomy", header = TRUE, sep = "\t")
## remove all the numbers in "", i.e in the bootstrap ##
tax$Taxonomy = gsub("[(0-9)]", "", tax$Taxonomy) 
## replaces the letter that are before the underscore ####
tax$Taxonomy = gsub("[a-z]__", "", tax$Taxonomy) 
## split based on semicolons ###
tax[c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')] <- colsplit(tax$Taxonomy,';',c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
### This removes taxonomy which is in column three ###
tax=tax[, -c(3)]
## Merged the two tables together
all=merge(otu,tax, by.x= "Group", by.y= "OTU")


## makes columns numeric ####
all[,2:61] = as.numeric(as.character(unlist(all[,2:61]))) 

### to check that it is numeric ##
class(all$CC7_36_t25_5_R1)
## Suming the 3 negatives ###
all$sumNeg=all$N_1+all$N_2+all$N_3
## Determining the contamination factor ###
all$Conta=(all$sumNeg/all$Size)*(100)
### Removing Contamination.. keeps everything <10 ####
noConta=subset(all, all$Conta<10)
conta=subset(all, all$Conta>10)


mat.fin=noConta[, c(1:55,59:61)]
#write.csv(mat.fin, "Aiptasia_diversity_NoConta.csv")

###nuevo##
##OTUs####
fam.wid=noConta[, c(1,2:33,35:51,53:55,59:61,67)] # col 35 and 52 are outliers 
fam.wid.agg=aggregate(fam.wid[, 2:56], by = list(fam.wid[, 57]), FUN =  sum)
fam.wid.agg$sum = rowSums(fam.wid.agg[,2:55], na.rm=TRUE) # excluding water

fam.top=subset(fam.wid.agg, fam.wid.agg$sum>=27200)#cutoff for top 20
fam.bot=subset(fam.wid.agg, fam.wid.agg$sum<=27200)# will be aggregated
fam.bot$sum=gsub(".*","Others", fam.bot$sum) # add tag to aggregate, all that are over 20 are named others
others=aggregate(fam.bot[, 2:56], by = list(fam.bot[, 57]), FUN =  sum)


fam.top1=fam.top[, -c(57)] #remove sum column
colnames(others)[1] <- "Family"
colnames(fam.top1)[1] <- "Family"
all.2 =rbind(fam.top1, others)


library(limma)
means=as.data.frame("avearrays"(all.2, ID=colnames(all.2), weights=NULL))
##jannas request
means[,2:ncol(means)] = as.numeric(as.character(unlist(means[,2:ncol(means)]))) 
means.2=as.matrix(means[,-c(1)])
test=sweep(means.2,2,colSums(means.2),`/`) 
row.names(test)=means[,c(1)]
write.table(test, "percetages_per_family.txt", sep = "\t")


all.l=melt(means, id.vars=c("Family"), variable.name = "Family", value.name = "Abundance")
colnames(all.l)=c("Family","Sample","Abundance")
all.l$Abundance=as.numeric(as.character(all.l$Abundance))
all.l$Family= factor(all.l$Family, as.character(all.l$Family))

P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")

a=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Family), data = all.l, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA reads", x="Samples") + scale_fill_manual(values=P21)
svg(filename = "Aiptasia_diversity_top20.svg",  width = 10, height = 5, pointsize = 12)
plot(a)
dev.off()

#non averaged barplots
#remove H2 39 t25 6 R2 and H2 42 t34 5 R1
all.l=melt(all.2, id.vars=c("Family"), variable.name = "Family", value.name = "Abundance")
colnames(all.l)=c("Family","Sample","Abundance")
all.l$Abundance=as.numeric(as.character(all.l$Abundance))
all.l$Family= factor(all.l$Family, as.character(all.l$Family))
map= read.csv("lookup.csv", header = TRUE) #lookup table
map.1=map[c(-58),]
unwanted=c("H2_t25_39_6_R2", "H2_t34_42_5_R1")
all.sub=subset(all.l, !Sample %in% unwanted)
all.sub$Sample=gsub("H2_", "aH2_", all.sub$Sample)


P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")

a=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Family), data = all.sub, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA reads", x="Samples") + scale_fill_manual(values=P21)
svg(filename = "Aiptasia_diversity_no_Averages.svg",  width = 10, height = 5, pointsize = 12)
plot(a)
dev.off()

