################
# load libraries
################
#library(gplots)
library(RColorBrewer)
library(dplyr)
library(ggrepel)

##############
## Scatterplot
##############
df = read.delim("Figure3_scatterplot/data4scatterplot.txt",header=T,stringsAsFactors = F)

g <- ggplot(df) +
  geom_point(aes(wf,thrips),fill="grey",color="black",shape=21,size=4) +
  theme_bw() +
  geom_label_repel(aes(x = wf,y=thrips,label=sample,fill=species)) +
  labs(x = "Tomato genotype rank for whitefly survival (low to high survival)",y = "Tomato genotype rank for thrips survival (low to high survival)") +
  scale_x_continuous(breaks=seq(0,19,1)) +
  scale_y_continuous(breaks=seq(0,19,1))

ggsave(filename = file.path("Figure3_scatterplot/","scatterplot.svg"),plot = g,width = 7,height = 5)
ggsave(filename = file.path("Figure3_scatterplot/","scatterplot.png"),plot = g,width = 7,height = 5)
