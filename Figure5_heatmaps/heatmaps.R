####!/usr/bin/env Rscript
##############################################################################
### Heatmaps and hierarchical clustering on metabolites and phenotypic classes
# author: marc galland
# contact: m.galland@uva.nl
# version 1.0 created on February 13, 2017
#####################################################################


################
# load libraries
################
library(RColorBrewer)
library(gplots)
library(ggplot2)
library("superheat") # download it from github

# result directory (use the working directory and the today's date)
outdir = file.path(getwd(),"Figure5_heatmaps/")

############################################
# load input dataframe and convert to logs  
############################################

terpenoids = read.delim("Figure5_heatmaps/pheno_terpenoids.tsv",header=T,stringsAsFactors = F,check.names = F)
acylsugars = read.delim("Figure5_heatmaps/pheno_acylsugars.tsv",header=T,stringsAsFactors = F,check.names = F)

# add row names
row.names(terpenoids)=terpenoids$sample;terpenoids$sample <- NULL
row.names(acylsugars)=acylsugars$sample;acylsugars$sample <- NULL

# I order the accessions based on whitefly or thrips median survival time
# add 1 to the zeros
# calculate log2 values
terpenoids.wf.ordered = terpenoids[order(terpenoids$wf),]
mat.terpenoids.wf.ordered = as.matrix(terpenoids.wf.ordered[,3:ncol(terpenoids.wf.ordered)])
mat.terpenoids.wf.ordered[mat.terpenoids.wf.ordered == 0] <- 1
mat.terpenoids.wf.ordered = log10(mat.terpenoids.wf.ordered)

acylsugars.wf.Ordered = acylsugars[order(acylsugars$wf),]
mat.acylsugars.wf.Ordered = as.matrix(acylsugars.wf.Ordered[,3:ncol(acylsugars.wf.Ordered)])
mat.acylsugars.wf.Ordered[mat.acylsugars.wf.Ordered == 0] <- 1
mat.acylsugars.wf.Ordered = log10(mat.acylsugars.wf.Ordered)

terpenoids.thripsOrdered = terpenoids[order(terpenoids$thrips),]
mat.terpenoids.thripsOrdered = as.matrix(terpenoids.thripsOrdered[,3:ncol(terpenoids.thripsOrdered)])
mat.terpenoids.thripsOrdered[mat.terpenoids.thripsOrdered == 0] <- 1
mat.terpenoids.thripsOrdered = log10(mat.terpenoids.thripsOrdered)

acylsugars.thripsOrdered = acylsugars[order(acylsugars$thrips),]
mat.acylsugars.thripsOrdered = as.matrix(acylsugars.thripsOrdered[,5:ncol(acylsugars.thripsOrdered)])
mat.acylsugars.thripsOrdered[mat.acylsugars.thripsOrdered == 0] <- 1
mat.acylsugars.thripsOrdered = log10(mat.acylsugars.thripsOrdered)

###################################
# Heatmaps with heatmap2 + side plots
###################################

## whitefly / terpenoids heatmap
pdf(file = "Figure5_heatmaps/terpenoids.wf.pdf",width = 7,height = 5)
heatmap.2(x = mat.terpenoids.wf.ordered,Rowv = F,Colv = F,scale = "none",revC = T,col=brewer.pal(9,"OrRd"),dendrogram = "none",margins = c(10,10),density.info = "none",cexCol = 0.3,cexRow = 0.8,key.title = "Scaled metabolite abundance (AU)",
          rowsep = 1:nrow(mat.terpenoids.wf.ordered),
          trace="none",
          sepcolor = "black",sepwidth = c(0.2,0))
dev.off()

## whitefly / acylsugars
pdf(file = "Figure5_heatmaps/acylsugars.wf.pdf",width = 7,height = 5)
heatmap.2(x = mat.acylsugars.wf.Ordered,Rowv = F,Colv = F,scale = "none",revC = T,col=brewer.pal(9,"BuPu"),
          trace = "none",dendrogram = "none",margins = c(10,10),density.info = "none",cexCol = 0.3,cexRow = 0.8,key.title = "Scaled metabolite abundance (AU)",
          rowsep = 1:nrow(mat.acylsugars.wf.Ordered),
          sepcolor = "black",
          sepwidth = c(0.2,0)
          )
dev.off()

# line plot whitefly median survival
genos = row.names(acylsugars.wf.Ordered)
wf = acylsugars.wf.Ordered[,"wf"]
tmp.df = data.frame(genos=genos,wf=wf,stringsAsFactors = F)
tmp.df$genos = factor(tmp.df$genos,levels = tmp.df$genos)
tmp.df = tmp.df[order(tmp.df$wf,decreasing = F),]
pdf(file = "Figure5_heatmaps/lineplot.wf.pdf",width = 4,height = 8)
ggplot(tmp.df,aes(x = genos,y = wf)) +
  geom_point(size=6) +
  coord_flip() +
  geom_line(group=1,size=2) +
  scale_y_continuous(limits = c(0,100)) +
  theme_bw()
dev.off()

########## thrips / volatiles
pdf(file = "Figure5_heatmaps/terpenoids.thrips.pdf",width = 7,height = 5)
heatmap.2(x = mat.terpenoids.thripsOrdered,Rowv = F,Colv = F,scale = "none",revC = T,col=brewer.pal(9,"OrRd"),dendrogram = "none",margins = c(10,10),density.info = "none",cexCol = 0.3,cexRow = 0.8,key.title = "Scaled metabolite abundance (AU)",
          rowsep = 1:nrow(mat.terpenoids.thripsOrdered),
          trace="none",
          sepcolor = "black",sepwidth = c(0.2,0))
dev.off()

########## thrips / acylsugars
pdf(file = "Figure5_heatmaps/acylsugars.thrips.pdf",width = 7,height = 5)
heatmap.2(x = mat.acylsugars.thripsOrdered,Rowv = F,Colv = F,scale = "none",revC = T,col=brewer.pal(9,"BuPu"),
          trace = "none",dendrogram = "none",margins = c(10,10),density.info = "none",cexCol = 0.3,cexRow = 0.8,key.title = "Scaled metabolite abundance (AU)",
          rowsep = 1:nrow(mat.acylsugars.wf.Ordered),
          sepcolor = "black",
          sepwidth = c(0.2,0)
)
dev.off()

# line plot thrips median survival
genos = row.names(acylsugars.thripsOrdered)
thrips = acylsugars.thripsOrdered[,"thrips"]
tmp.df = data.frame(genos=genos,thrips=thrips,stringsAsFactors = F)
tmp.df$genos = factor(tmp.df$genos,levels = tmp.df$genos)
tmp.df = tmp.df[order(tmp.df$thrips,decreasing = F),]
pdf(file = "Figure5_heatmaps/lineplot.thrips.pdf",width = 4,height = 8)
ggplot(tmp.df,aes(x = genos,y = thrips)) +
  geom_point(size=6) +
  coord_flip() +
  geom_line(group=1,size=2) +
  scale_y_continuous(limits = c(0,15)) +
  theme_classic()
dev.off()


#####################################################
# Write used matrices to make the supplemental tables
#####################################################
write.table(x = mat.terpenoids.wf.ordered,   file = "Figure5_heatmaps/tables/TableS1_terpenoids_wf.tsv",quote = F,row.names = T,sep = "\t")
write.table(x = mat.terpenoids.thripsOrdered,file = "Figure5_heatmaps/tables/TableS2_terpenoids_thrips.tsv",quote = F,row.names = T,sep = "\t")
write.table(x = mat.acylsugars.wf.Ordered,   file = "Figure5_heatmaps/tables/TableS3_acylsugars_wf.tsv",quote = F,row.names = T,sep = "\t")
write.table(x = mat.acylsugars.thripsOrdered,file = "Figure5_heatmaps/tables/TableS4_acylsugars_thrips.tsv",quote = F,row.names = T,sep = "\t")


