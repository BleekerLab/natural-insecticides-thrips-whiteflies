# library
library(tidyverse)
library(RColorBrewer)
library(superheat) # install the dev version >> devtools::install_github("rlbarter/superheat")


# load metabolites data
metabolites = read.delim("Figure8_volatiles_F2/20180625_F2metabolite_table.tsv",header=T,stringsAsFactors = F,check.names = F)
head(metabolites)


# load whitefly data
wf = read.delim("Figure7_F2_whitefly_survival/20170912_wf_bioassay_F2.tsv",header = T,stringsAsFactors = F)


####### plot heatmap ########

# calculate mean WF survival and order
wf = wf %>% dplyr::mutate(percentage = alive/(alive+dead)*100)
wf = wf %>% group_by(line) %>% summarise(survival = mean(percentage)) 
wf =  wf[order(wf$survival),]

# add whitefly survival values to the metabolite table
# remove missing lines
metabolites = dplyr::left_join(metabolites,wf,by="line")
metabolites = na.omit(metabolites)

# reorder metabolites according to increasing wf survival values
metabolites.sorted = metabolites[order(metabolites$survival,decreasing = F),]

# prepare matrix for heatmaps
metabolites.sorted.nosurvival = metabolites.sorted
metabolites.sorted.nosurvival$survival = NULL
mat = as.matrix(metabolites.sorted.nosurvival[,2:ncol(metabolites.sorted.nosurvival)])
row.names(mat)=metabolites.sorted.nosurvival$line
mat[mat == 0] <- 1
mat = log2(mat)

# superheatmap
png(file.path("Figure8_volatiles_F2/","F2_volatiles_with_wf.png"),height=900,width=1600)
superheat(mat,
          scale=F,
          smooth.heat = F,
          # metabolites
          bottom.label.text.angle = 90,
          bottom.label.text.size = 6,
          column.title = "Volatiles",
          # colors
          heat.pal=brewer.pal(9,"OrRd"),
          grid.hline = F,
          grid.vline = F,
          # scatterplot on the side ()
          yr=metabolites.sorted$survival,
          yr.axis.name = "Whitefly survival (%)",
          yr.axis.size = 10,
          yr.axis.name.size = 10,
          yr.plot.type = "scatterline",
          yr.obs.col = rep("grey",nrow(mat)),
          yr.point.size = 10,
          yr.num.ticks = 9,
          # genotype info
          left.label.text.alignment = "left",
          row.title = "Genotype",
          legend = T,
          pretty.order.rows = F,
          pretty.order.cols = F,
          force.bottom.label = T,
          # legend
          legend.breaks = c(-5,0,5,10,15,20),
          legend.text.size = 20
)
dev.off()


# pdf(file.path("Figure8_volatiles_F2/","F2_volatiles_with_wf.pdf"))
# superheat(mat,
#           scale=F,
#           smooth.heat = F,
#           # metabolites
#           bottom.label.text.angle = 90,
#           bottom.label.text.size = 2,
#           column.title = "Volatiles",
#           # colors
#           heat.pal=brewer.pal(9,"OrRd"),
#           grid.hline = F,
#           grid.vline = F,
#           # scatterplot on the side ()
#           yr=metabolites.sorted$survival,
#           yr.axis.name = "Whitefly survival (%)",
#           yr.axis.size = 10,
#           yr.axis.name.size = 10,
#           yr.plot.type = "scatterline",
#           yr.obs.col = rep("grey",nrow(mat)),
#           yr.point.size = 5,
#           yr.num.ticks = 9,
#           # genotype info
#           left.label.text.alignment = "left",
#           row.title = "Genotype",
#           legend = T,
#           pretty.order.rows = F,
#           pretty.order.cols = T,
#           force.bottom.label = T,
#           heat.lim = c(-2,6),
#           legend.breaks = c(-2,0,2,4,6)
# )
# dev.off()
