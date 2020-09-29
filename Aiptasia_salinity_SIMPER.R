
######SIMPER############
library("vegan")
library("pracma")

mat=t(read.table("OTU_noConta", row.names = 1))

mat.n=sweep(mat.m,2,colSums(mat.m),`/`)
sqrtm(mat.n)

map= read.csv("lookup.csv", header = TRUE) #lookup table
groups = map[-c(58:59), c(2)] # assign groups

simp=simper(mat.n, groups, permutations = 0, trace = FALSE)

results=summary(simp)

write.table(results, file = "SIMPER_results.txt")




#########SIMPER##########
library("vegan")
library("pracma")


h2.t25.met=subset(h2.groups, Temperature=="25")
h2.t34.met=subset(h2.groups, Temperature=="34")
cc7.t25.met=subset(cc7.groups, Temperature=="25")
cc7.t34.met=subset(cc7.groups, Temperature=="34")

h2.t25=subset(h2, rownames(h2) %in% h2.t25.met$ID)
h2.t34=subset(h2, rownames(h2) %in% h2.t34.met$ID)
cc7.t25=subset(cc7, rownames(cc7) %in% cc7.t25.met$ID)
cc7.t34=subset(cc7, rownames(cc7) %in% cc7.t34.met$ID)

##pretty simper

simper.pretty = function(x, metrics, interesting, perc_cutoff, low_cutoff, low_val, output_name){
  library(vegan)
  #handling otu tables for taxa levels
  save=list(0)
  if(grepl("Otu", colnames(x)[1])!=TRUE){
    #converts output from A.Neumann Taxonomy script
    save=list(58)
    x=as.data.frame(t(x))
    orig_names=colnames(x)
    new_names=list()
    l=1
    for(n in colnames(x)){
      ifelse((l<10), (colnames(x)[l]=c(paste0('Otu000',c(l)))), (colnames(x)[l]=c(paste0('Otu00',c(l))))) 
      new_names=append(new_names, colnames(x)[l])
      l=l+1
    }
    orig_names=as.data.frame(orig_names, row.names = as.character(new_names))
  }
  #running simper
  for(variables in interesting){
    test_1=with(metrics, simper(x, metrics[[variables]]))
    #parsing through simper output, saving desired info to table
    for(name in names(test_1)){
      testmx=matrix(ncol=length(interesting))
      testmx=cbind(test_1[[name]]$ord,test_1[[name]]$cusum)
      sorted=testmx[order(testmx[,1]),]
      sorted=cbind(sorted,test_1[[name]]$species)
      sorted=sorted[order(sorted[,2]),]
      t=matrix(sorted[sorted[,2]<=perc_cutoff,],ncol=3)
      i=nrow(t)
      #converting percents to percent of whole
      while(i>1){
        t[i,2]=as.character(as.numeric(t[i,2])-as.numeric(t[i-1,2]))
        i=i-1
      }
      t[,1]=name
      write.table(t,file=paste(output_name,'_simper.csv',sep=""), append=TRUE, sep=",", col.names = FALSE)
    }}
  y=read.table(paste(output_name,'_simper.csv',sep=""), header=FALSE,sep=",",fill = TRUE,row.names = NULL)
  file.remove(paste(output_name,'_simper.csv',sep = ""))
  y=y[-c(1)]
  colnames(y) = c("Comparison", "SIMPER", "OTU")
  #removing results lower than low cutoff
  if(low_cutoff=='y'){
    y=y[!(as.numeric(as.character(y$SIMPER))<low_val),]
  }
  #prevents changing of colnames if OTU table
  if(58 %in% save){
    y$OTU=orig_names[match(y$OTU, rownames(orig_names)),1]
  }
  write.csv(y,file=paste(output_name,'_clean_simper.csv', sep=''))
}


simper.pretty(h2.t25, h2.t25.met, c('Salinity'), perc_cutoff=0.5, low_cutoff = 'y', low_val=0.01, 'h2_t25')
h2_t25_sim.otus=read.csv("h2_t25_clean_simper.csv")

simper.pretty(h2.t34, h2.t34.met, c('Salinity'), perc_cutoff=0.5, low_cutoff = 'y', low_val=0.01, 'h2_t34')
h2_t34_sim.otus=read.csv("h2_t34_clean_simper.csv")

simper.pretty(cc7.t25, cc7.t25.met, c('Salinity'), perc_cutoff=0.5, low_cutoff = 'y', low_val=0.01, 'cc7_t25')
cc7_t25_sim.otus=read.csv("cc7_t25_clean_simper.csv")

simper.pretty(cc7.t34, cc7.t34.met, c('Salinity'), perc_cutoff=0.5, low_cutoff = 'y', low_val=0.01, 'cc7_t34')
cc7_t34_sim.otus=read.csv("cc7_t34_clean_simper.csv")


#Heatmaps with selected species

h2.t25.simp.otus.m=as.matrix(t(subset(h2.t25, select=c(colnames(h2.t25) %in% h2_t25_sim.otus$OTU))))
h2.t25.simp.otus.r = sqrt(h2.t25.simp.otus.m)
h2.t25.simp.otus.r = log(h2.t25.simp.otus.m+1)
colnames(h2.t25.simp.otus.r) = h2.t25.met$Salinity[match(colnames(h2.t25.simp.otus.r),h2.t25.met$ID)] # change colnames by lookup table
salinities = h2.t25.met$Salinity[match(colnames((h2.t25.simp.otus.r),h2.t25.met$Temperature))]
colorCode = as.matrix(h2.t25.met[,c(6,7)])
colorCode = as.matrix(h2.t25.met$Salinity_col)
## add a mapping heatmap showing SIMPER contribution on row variables (OTUs)
# h2_t25_otus_wid=reshape(h2_t25_sim.otus[,-c(1)], idvar = "OTU", timevar = "Comparison", direction = "wide")
# h2_t25_otus_wid[is.na(h2_t25_otus_wid)] = 0
# contribution =h2_t25_otus_wid[match(rownames(h2.t25.simp.otus.r), h2_t25_otus_wid$OTU),]
# contribution = h2.t25.met$Salinity[match(rownames((h2.t25.simp.otus.r),h2_t25_sim.otus$OTU))]
# colorCode2=as.matrix(contribution[,-c(1)])


# library("gplots") # for the heatmap.2 function
# library("devtools")
# par(mfrow=c(1,1))
# 
# #Define custom dist and hclust functions for use with heatmaps
# mydist=function(c) {dist(c,method="euclidian")} 
# 
# 
# #Load latest version of heatmap.3 function
# source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
# 
# require(gtools)
# require(RColorBrewer)
# cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
# colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# distCor <- function(x) as.dist(1-cor(t(x)))
# hclustAvg <- function(x) hclust(x, method="average")
# hclustAvg <-function(x)  hclust(as.dist(1-cor(t(x), method="pearson")), method="complete")
# #svg(filename = "H2_t25_hm1.svg",  width = 7, height = 5, pointsize = 12)
# png(filename = "H2_t25_hm1.png", width = 1500, height = 1500, units = "px", pointsize = 12, bg = "NA",  res = 300)
# heatmap.3(h2.t25.simp.otus.r, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", dendrogram="row", margins=c(5,10),
#           Rowv=TRUE, Colv=FALSE, ColSideColors=colorCode,  
#           density.info="none", trace="none", labCol=FALSE,  cexRow=1, col = colorRampPalette( rev(brewer.pal(3, "RdBu")) )(255), ColSideColorsSize=2, RowSideColorsSize=2, KeyValueName="Abundance", main = "H2_t25")
# dev.off()

h2.t34.simp.otus.m=as.matrix(t(subset(h2.t34, select=c(colnames(h2.t34) %in% h2_t34_sim.otus$OTU))))
#h2.t34.simp.otus.r = sqrt(h2.t34.simp.otus.m)
h2.t34.simp.otus.r = log(h2.t34.simp.otus.m+1)
colnames(h2.t34.simp.otus.r) = h2.t34.met$Salinity[match(colnames(h2.t34.simp.otus.r),h2.t34.met$ID)] # change colnames by lookup table
#salinities = h2.t34.met$Salinity[match(colnames((h2.t34.simp.otus.r),h2.t34.met$Salinity))]
colorCode = as.matrix(h2.t34.met[,c(6,7)])
colorCode = as.matrix(h2.t34.met$Salinity_col)
png(filename = "H2_t34_hm1.png",   width = 1500, height = 1500, units = "px", pointsize = 12, bg = "NA",  res = 300)
heatmap.3(h2.t34.simp.otus.r, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", dendrogram="row", margins=c(5,10),
          Rowv=TRUE, Colv=FALSE, ColSideColors=colorCode,  
          density.info="none", trace="none", labCol=FALSE,  cexRow=1, col = colorRampPalette( rev(brewer.pal(3, "RdBu")) )(255), ColSideColorsSize=2, RowSideColorsSize=2, KeyValueName="Abundance", main = "H2_t34")
dev.off()

cc7.t34.simp.otus.m=as.matrix(t(subset(cc7.t34, select=c(colnames(cc7.t34) %in% cc7_t34_sim.otus$OTU))))
cc7.t34.simp.otus.r = sqrt(cc7.t34.simp.otus.m)
cc7.t34.simp.otus.r = log(cc7.t34.simp.otus.m+1)
colnames(cc7.t34.simp.otus.r) = cc7.t34.met$Salinity[match(colnames(cc7.t34.simp.otus.r),cc7.t34.met$ID)] # change colnames by lookup table
#salinities = cc7.t34.met$Salinity[match(colnames((cc7.t34.simp.otus.r),cc7.t34.met$Temperature))]
colorCode = as.matrix(cc7.t34.met[,c(6,7)])
colorCode = as.matrix(cc7.t34.met$Salinity_col)
png(filename = "cc7_t34_hm1.png",  width = 1500, height = 1500, units = "px", pointsize = 12, bg = "NA",  res = 300)
heatmap.3(cc7.t34.simp.otus.r, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", dendrogram="row", margins=c(5,10),
          Rowv=TRUE, Colv=FALSE, ColSideColors=colorCode,  
          density.info="none", trace="none", labCol=FALSE,  cexRow=1, col = colorRampPalette( rev(brewer.pal(3, "RdBu")) )(255), ColSideColorsSize=2, RowSideColorsSize=2, KeyValueName="Abundance", main = "CC7_t34")
dev.off()

cc7.t25.simp.otus.m=as.matrix(t(subset(cc7.t25, select=c(colnames(cc7.t25) %in% cc7_t25_sim.otus$OTU))))
cc7.t25.simp.otus.r = sqrt(cc7.t25.simp.otus.m)
cc7.t25.simp.otus.r = log(cc7.t25.simp.otus.m+1)
colnames(cc7.t25.simp.otus.r) = cc7.t25.met$Salinity[match(colnames(cc7.t25.simp.otus.r),cc7.t25.met$ID)] # change colnames by lookup table
#salinities = cc7.t25.met$Salinity[match(colnames((cc7.t25.simp.otus.r),cc7.t25.met$Temperature))]
colorCode = as.matrix(cc7.t25.met[,c(6,7)])
colorCode = as.matrix(cc7.t25.met$Salinity_col)
png(filename = "cc7_t25_hm1.png",   width = 1500, height = 1500, units = "px", pointsize = 12, bg = "NA",  res = 300)
heatmap.3(cc7.t25.simp.otus.r, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", dendrogram="row", margins=c(5,10),
          Rowv=TRUE, Colv=FALSE, ColSideColors=colorCode,  
          density.info="none", trace="none", labCol=FALSE,  cexRow=1, col = colorRampPalette( rev(brewer.pal(3, "RdBu")) )(255), ColSideColorsSize=2, RowSideColorsSize=2, KeyValueName="Abundance", main = "CC7_t25")
dev.off()


#means

h2.t25.simp.otus.means=as.matrix("avearrays"(h2.t25.simp.otus.r, ID=colnames(h2.t25.simp.otus.r), weights=NULL))
# #png(filename = "H2_t25_averages.png", width = 1500, height = 1500, units = "px", pointsize = 12, bg = "NA",  res = 300)
# heatmap.3(h2.t25.simp.otus.means, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", dendrogram="row", margins=c(5,10),
#           Rowv=TRUE, Colv=FALSE, density.info="none", trace="none", labCol=FALSE,  cexRow=1, col = colorRampPalette( rev(brewer.pal(3, "RdBu")) )(255), ColSideColorsSize=2, RowSideColorsSize=2, KeyValueName="Abundance", main = "H2_t25_ave")
# dev.off()

h2.t34.simp.otus.means=as.matrix("avearrays"(h2.t34.simp.otus.r, ID=colnames(h2.t34.simp.otus.r), weights=NULL))
# png(filename = "H2_t34_averages.png", width = 1500, height = 1500, units = "px", pointsize = 12, bg = "NA",  res = 300)
# heatmap.3(h2.t34.simp.otus.means, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", dendrogram="row", margins=c(5,10),
#           Rowv=TRUE, Colv=FALSE, density.info="none", trace="none", labCol=FALSE,  cexRow=1, col = colorRampPalette( rev(brewer.pal(3, "RdBu")) )(255), ColSideColorsSize=2, RowSideColorsSize=2, KeyValueName="Abundance", main = "H2_t34_ave")
# dev.off()

cc7.t25.simp.otus.means=as.matrix("avearrays"(cc7.t25.simp.otus.r, ID=colnames(cc7.t25.simp.otus.r), weights=NULL))
# png(filename = "cc7_t25_averages.png", width = 1500, height = 1500, units = "px", pointsize = 12, bg = "NA",  res = 300)
# heatmap.3(cc7.t25.simp.otus.means, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", dendrogram="row", margins=c(5,10),
#           Rowv=TRUE, Colv=FALSE, density.info="none", trace="none", labCol=FALSE,  cexRow=1, col = colorRampPalette( rev(brewer.pal(3, "RdBu")) )(255), ColSideColorsSize=2, RowSideColorsSize=2, KeyValueName="Abundance", main = "CC7_t25_ave")
# dev.off()

cc7.t34.simp.otus.means=as.matrix("avearrays"(cc7.t34.simp.otus.r, ID=colnames(cc7.t34.simp.otus.r), weights=NULL))
# png(filename = "CC7_t34_averages.png", width = 1500, height = 1500, units = "px", pointsize = 12, bg = "NA",  res = 300)
# heatmap.3(cc7.t34.simp.otus.means, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", dendrogram="row", margins=c(5,10),
#           Rowv=TRUE, Colv=FALSE, density.info="none", trace="none", labCol=FALSE,  cexRow=1, col = colorRampPalette( rev(brewer.pal(3, "RdBu")) )(255), ColSideColorsSize=2, RowSideColorsSize=2, KeyValueName="Abundance", main = "CC7_t34_ave")



# Clustering using pearson correlation 

mycol <- colorpanel(50, "black", "white", "green") # or try redgreen(75)
mycol2=colorRampPalette( rev(brewer.pal(3, "BrBG")) )(255)
mycol3 <- colorpanel(25, "black", "white", "green")
## Plot heatmaps

hr_h2_25 <- hclust(as.dist(1-cor(t(h2.t25.simp.otus.means), method="pearson")), method="complete")
mycl <- cutree(hr_h2_25, h=max(hr_h2_25$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
svg(filename = "H2_t25_ave_pearson.svg",  width = 5, height = 7, pointsize = 12)
#png(filename = "H2_t25_ave_pearson.png", width = 1500, height = 1500, units = "px", pointsize = 12, bg = "NA",  res = 300)
heatmap.2(h2.t25.simp.otus.means, Rowv=as.dendrogram(hr_h2_25 ),  Colv=TRUE, col=mycol, scale="row", density.info="none", trace="none",  main = "H2 25ºC")
dev.off()

pdf(file = "H2_t25_ave_pearson.pdf",  width = 5, height = 7, pointsize = 12)

hr_h2_34 <- hclust(as.dist(1-cor(t(h2.t34.simp.otus.means), method="pearson")), method="complete")
mycl <- cutree(hr_h2_34, h=max(hr_h2_34$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
svg(filename = "H2_t34_ave_pearson.svg",  width = 5, height = 7, pointsize = 12)
heatmap.2(h2.t34.simp.otus.means, Rowv=as.dendrogram(hr_h2_34),  Colv=TRUE, col=mycol, scale="row", density.info="none", trace="none",  main = "H2 34ºC")
dev.off()

hr_cc7_25 <- hclust(as.dist(1-cor(t(cc7.t25.simp.otus.means), method="pearson")), method="complete")
mycl <- cutree(hr_h2_25, h=max(hr_cc7_25$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
svg(filename = "CC7_t25_ave_pearson.svg",  width =5, height = 7.5, pointsize = 12)
heatmap.2(cc7.t25.simp.otus.means, Rowv=as.dendrogram(hr_cc7_25),  Colv=TRUE, col=mycol, scale="row", density.info="none", trace="none",  main = "CC7 25ºC")
dev.off()

hr_cc7_34 <- hclust(as.dist(1-cor(t(cc7.t34.simp.otus.means), method="pearson")), method="complete")
mycl <- cutree(hr_cc7_34, h=max(hr_cc7_34 $height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
svg(filename = "CC7_t34_ave_pearson.svg",  width = 5, height = 8.81, pointsize = 12)
heatmap.2(cc7.t34.simp.otus.means, Rowv=as.dendrogram(hr_cc7_34),  Colv=TRUE, col=mycol3, scale="row", density.info="none", trace="none",  main = "CC7 34ºC")
dev.off()


library(gridGraphics)
library(grid)
library(gridSVG)

grid.newpage() 
pushViewport(viewport(0,.6, .3, .5, just=c("left", "bottom"), name="A"))
grid.echo(function() { heatmap.2(h2.t25.simp.otus.means, Rowv=as.dendrogram(hr_h2_25 ),  Colv=FALSE, col=mycol, scale="row", density.info="none", trace="none",  main = "", key = FALSE, cexRow=0.7, cexCol =0.7)
}, newpage=FALSE) 
popViewport() 
pushViewport(viewport(.6, .6, .3, .5, just=c("left", "bottom"))) 
grid.echo(function() { heatmap.2(h2.t34.simp.otus.means, Rowv=as.dendrogram(hr_h2_34),  Colv=FALSE, col=mycol, scale="row", density.info="none", trace="none",  main = "", key = FALSE, cexRow=0.7, cexCol =0.7)
}, newpage=FALSE) 
popViewport() 
pushViewport(viewport(0, .05, .3, .53, just=c("left", "bottom"))) 
grid.echo(function() { heatmap.2(cc7.t25.simp.otus.means, Rowv=as.dendrogram(hr_cc7_25),  Colv=FALSE, col=mycol, scale="row", density.info="none", trace="none",  main = "", key = FALSE, cexRow=0.7, cexCol =0.7)
}, newpage=FALSE) 
popViewport() 
pushViewport(viewport(.6, 0, .3, .61, just=c("left", "bottom"))) 
grid.echo(function() { heatmap.2(cc7.t34.simp.otus.means, Rowv=as.dendrogram(hr_cc7_34),  Colv=FALSE, col=mycol, scale="row", density.info="none", trace="none",  main = "", key = FALSE, cexRow=0.7, cexCol =0.7)
}, newpage=FALSE) 
popViewport() 
#grid.export("heatmap_grid.svg")


#Venn Diagrams
venn_h2_25 = rownames(h2.t25.simp.otus.means)
write.table(venn_h2_25, quote = FALSE, file ="venn_diagram_H2_t25", col.names=FALSE, row.names = FALSE)

venn_h2_34 = rownames(h2.t34.simp.otus.means)
write.table(venn_h2_34, quote = FALSE, file ="venn_diagram_H2_t34", col.names=FALSE, row.names = FALSE)

venn_cc7_25 = rownames(cc7.t25.simp.otus.means)
write.table(venn_cc7_25, quote = FALSE, file ="venn_diagram_cc7_t25", col.names=FALSE, row.names = FALSE)

venn_cc7_34 = rownames(cc7.t34.simp.otus.means)
write.table(venn_cc7_34, quote = FALSE, file ="venn_diagram_cc7_t34", col.names=FALSE, row.names = FALSE)


#Table taxonomy

library(reshape2)
## Reads the taxonomy csv file,
tax = read.csv("../janna.final.taxonomy", header = TRUE, sep = "\t")
tax$Taxonomy = gsub("[(0-9)]", "", tax$Taxonomy) 
tax$Taxonomy = gsub("[a-z]__", "", tax$Taxonomy) 
tax[c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')] <- colsplit(tax$Taxonomy,';',c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
tax=tax[, -c(3)]

all.otu.simp=rbind(h2_t25_sim.otus, h2_t34_sim.otus, cc7_t25_sim.otus, cc7_t34_sim.otus)
simper.taxa=subset(tax, tax$OTU %in% all.otu.simp$OTU)
write.table(simper.taxa, quote = FALSE,  row.names = FALSE, file = "simper_taxa.txt", sep = "\t")

fast = read.csv("../janna.final.fasta2", header = FALSE, sep = "\t")
simper.fast=subset(fast, fast$V1 %in% all.otu.simp$OTU)
write.table(simper.fast, quote = FALSE,  row.names = FALSE,  col.names =  FALSE, file = "simper_fasta.txt", sep = "\t")
