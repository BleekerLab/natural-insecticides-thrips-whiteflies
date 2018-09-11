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
library("superheat") # download it from github

setwd("~/surfdrive/BleekerLab/00.manuscripts/01.natural_insecticides_from_wild_tomatoes/natural-insecticides-wild-tomatoes-paper/Figure5_heatmaps/")

# result directory (use the working directory and the today's date)
outdir = file.path(getwd(),strsplit(x = format(Sys.time()),split = " ")[[1]][1])
dir.create(outdir,showWarnings = F,recursive = T)

############################################
# load input dataframe and convert to logs  
############################################

terpenoids = read.delim("pheno_terpenoids.txt",header=T,stringsAsFactors = F,check.names = F)
acylsugars = read.delim("pheno_acylsugars.txt",header=T,stringsAsFactors = F,check.names = F)

# add row names
row.names(terpenoids)=terpenoids$sample;terpenoids$sample <- NULL
row.names(acylsugars)=acylsugars$sample;acylsugars$sample <- NULL

#################################
# Heatmap terpenoids and whitefly
#################################  
# I order the accessions based on whitefly median survival time
# add 1 to the zeros
# calculate log2 values
terpenoids.wf.ordered = terpenoids[order(terpenoids$wf),]
mat.terpenoids.wf.ordered = as.matrix(terpenoids.wf.ordered[,3:ncol(terpenoids.wf.ordered)])
mat.terpenoids.wf.ordered[mat.terpenoids.wf.ordered == 0] <- 1
mat.terpenoids.wf.ordered = log10(mat.terpenoids.wf.ordered)


#height = 700,width=900
png(filename = file.path(outdir,"terpenoids.wf.png"),width = 900,height = 700)
superheat(mat.terpenoids.wf.ordered,
          scale=F,
          smooth.heat = F,
          # metabolites
          bottom.label.text.angle = 90,
          bottom.label.text.size = 1.5,
          column.title = "Metabolites",
          # colors
          heat.pal=brewer.pal(9,"OrRd"),
          grid.hline = F,
          grid.vline = F,
          # scatterplot on the side ()
          yr=terpenoids.wf.ordered$wf,
          yr.axis.name = "Whitefly median survival (%)",
          yr.axis.size = 10,
          yr.axis.name.size = 10,
          yr.plot.type = "scatterline",
          yr.obs.col = rep("grey",nrow(mat.terpenoids.wf.ordered)),
          yr.point.size = 5,
          yr.num.ticks = 10,
          # genotype info
          left.label.text.alignment = "left",
          row.title = "Genotype",
          legend = T,
          pretty.order.cols = T
)
dev.off()

#################################
# Heatmap terpenoids and thrips
#################################  
# I order the accessions based on thrips median survival time
# add 1 to the zeros
# calculate log2 values
terpenoids.thripsOrdered = terpenoids[order(terpenoids$thrips),]
mat.terpenoids.thripsOrdered = as.matrix(terpenoids.thripsOrdered[,3:ncol(terpenoids.thripsOrdered)])
mat.terpenoids.thripsOrdered[mat.terpenoids.thripsOrdered == 0] <- 1
mat.terpenoids.thripsOrdered = log10(mat.terpenoids.thripsOrdered)

# plot it
png(file.path(outdir,"terpenoids.thrips.png"),height=700,width=900)
superheat(mat.terpenoids.thripsOrdered,
          scale=F,
          smooth.heat = F,
          # metabolites
          bottom.label.text.angle = 90,
          bottom.label.text.size = 1.5,
          column.title = "Metabolites",
          # colors
          heat.pal=brewer.pal(9,"OrRd"),
          grid.hline = F,
          grid.vline = F,
          # scatterplot on the side ()
          yr=terpenoids.thripsOrdered$thrips,
          yr.axis.name = "Thrips median survival time (days)",
          yr.axis.size = 10,
          yr.axis.name.size = 10,
          yr.plot.type = "scatterline",
          yr.obs.col = rep("grey",nrow(mat.terpenoids.thripsOrdered)),
          yr.point.size = 5,
          yr.num.ticks = 9,
          # genotype info
          left.label.text.alignment = "left",
          row.title = "Genotype",
          legend = T,
          pretty.order.cols = T
)
dev.off()


#################################
# Heatmap acylsugars and whitefly
#################################  
acylsugars.wf.Ordered = acylsugars[order(acylsugars$wf),]
mat.acylsugars.wf.Ordered = as.matrix(acylsugars.wf.Ordered[,3:ncol(acylsugars.wf.Ordered)])
mat.acylsugars.wf.Ordered[mat.acylsugars.wf.Ordered == 0] <- 1
mat.acylsugars.wf.Ordered = log10(mat.acylsugars.wf.Ordered)

png(file.path(outdir,"acylsugars.wf.png"),height=700,width=900)
superheat(mat.acylsugars.wf.Ordered,
          scale=F,
          smooth.heat = F,
          # metabolites
          bottom.label.text.angle = 90,
          bottom.label.text.size = 1.5,
          column.title = "Metabolites",
          # colors
          heat.pal=brewer.pal(9,"BuPu"),
          grid.hline = F,
          grid.vline = F,
          # scatterplot on the side ()
          yr=acylsugars$wf,
          yr.axis.name = "Whitefly median survival rate (%)",
          yr.axis.size = 10,
          yr.axis.name.size = 10,
          yr.plot.type = "scatterline",
          yr.obs.col = rep("grey",nrow(mat.acylsugars.wf.Ordered)),
          yr.point.size = 5,
          yr.num.ticks = 9,
          # genotype info
          left.label.text.alignment = "left",
          row.title = "Genotype",
          legend = T,
          pretty.order.cols = T,
          force.bottom.label = T,
          #legend.breaks = c(0,2,4,6,8)
)
dev.off()

#################################
# Heatmap acylsugars and thrips
#################################  
# if I order the accessions based on thrips median survival time
# add 1 to the zeros
# calculate log2 values
acylsugars.thripsOrdered = acylsugars[order(acylsugars$thrips),]
mat.acylsugars.thripsOrdered = as.matrix(acylsugars.thripsOrdered[,5:ncol(acylsugars.thripsOrdered)])
mat.acylsugars.thripsOrdered[mat.acylsugars.thripsOrdered == 0] <- 1
mat.acylsugars.thripsOrdered = log10(mat.acylsugars.thripsOrdered)

png(file.path(outdir,"acylsugars.thrips.png"),height=700,width=900)
superheat(mat.acylsugars.thripsOrdered,
          scale=F,
          smooth.heat = F,
          # metabolites
          bottom.label.text.angle = 90,
          bottom.label.text.size = 2,
          column.title = "Metabolites",
          # colors
          heat.pal=brewer.pal(9,"BuPu"),
          grid.hline = F,
          grid.vline = F,
          # scatterplot on the side ()
          yr=acylsugars.thripsOrdered$thrips,
          yr.axis.name = "Thrips median survival (days)",
          yr.axis.size = 10,
          yr.axis.name.size = 10,
          yr.plot.type = "scatterline",
          yr.obs.col = rep("grey",nrow(mat.acylsugars.thripsOrdered)),
          yr.point.size = 5,
          yr.num.ticks = 9,
          # genotype info
          left.label.text.alignment = "left",
          row.title = "Genotype",
          legend = T,
          pretty.order.cols = T,
          force.bottom.label = T,
          legend.breaks = c(0,2,4,6,8)
)
dev.off()