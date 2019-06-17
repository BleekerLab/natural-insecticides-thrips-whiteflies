################
# load libraries
################
#library(gplots)
library(RColorBrewer)
library(dplyr)
library(ggrepel)

setwd("~/surfdrive/BleekerLab/00.manuscripts/01.natural_insecticides_from_wild_tomatoes/01.Figures/_Figure4_scatterplot/")
outdir = getwd()

##############
## Scatterplot
##############
df.scatterplot = read.delim("data4scatterplot.txt",header=T,stringsAsFactors = F)
g <- ggplot(df.scatterplot) + 
  geom_point(aes(wf,thrips),fill="grey",color="black",shape=21,size=4) + 
  theme_bw() +
  geom_label_repel(aes(x = wf,y=thrips,label=sample,fill=species)) +
  labs(x = "Tomato genotype rank for whitefly survival (low to high survival)",y = "Tomato genotype rank for thrips survival (low to high survival)") +
  scale_x_continuous(breaks=seq(0,19,1)) +
  scale_y_continuous(breaks=seq(0,19,1)) 
#guides(fill=FALSE)
print(g)
ggsave(filename = file.path(outdir,"scatterplot.svg"),plot = g,width = 7,height = 5)
