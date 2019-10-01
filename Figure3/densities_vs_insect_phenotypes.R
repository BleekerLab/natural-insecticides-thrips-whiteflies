library(RColorBrewer)
library(ggrepel)
library(tidyverse)
library(Rmisc)

#############################
# Custom theme for plotting #
#############################

my.theme = theme(axis.text.x = element_text(color = "black", size = 6),
                 axis.text.y = element_text(color = "black", size = 6),
                 axis.title.x = element_text(color = "black", size = 8),
                 axis.title.y = element_text(color = "black", size = 8)
)

##################################################
# Load datasets
##################################################

# trichome density data
density = read.delim("FigureS4/trichome.counts.processed.tsv", header =T)

# whitefly and thrips ranks based on toxicity
pheno = read.delim("Figure2/data4scatterplot.txt", header = T)
pheno = separate(pheno, col = sample, into = c("x", "y", "accession")) %>% #separate "sample" into accession names only so it can be used for fusing the datasets
  select(., c("accession", "wf", "thrips"))

################################################################################################
# Calculate the average trichome densities per genotype (abaxial and adaxial are taken together)
density.avg = summarySE(
  density, 
  measurevar = "density", 
  groupvars = c("genotype", "accession", "species" ,"type", "color")
  )

#####################################
# merge the two datasets to a new one
#####################################
density.pheno  = left_join(density.avg, pheno, by = "accession")
density.pheno$accession = as.factor(density.pheno$accession)

# Label the trichome classes nicely
density.pheno$type = factor(density.pheno$type, levels = c("non.glandular", "typeIandIV", "typeVI"),
                            labels = c("Non glandular", "Type I/IV", "Type VI"))

# Plot whiteflies vs. trichome densities
p.whitelfies = 
density.pheno %>%
ggplot() +
  geom_point(aes(x = wf, y = density),fill="grey",color="black",shape=21,size=2) +
  theme_bw() +
  geom_label_repel(aes(x = wf, y=density,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  facet_wrap(~type, scale = "free", ncol = 1)+
  labs(x = "Tomato genotype rank for whitefly survival (low to high survival)",y = "Trichome density (trichomes / mm2 leaf)") +
  scale_x_continuous(breaks=seq(0,19,1))+
  my.theme

# Plot thrips vs. trichome densities
p.thrips = 
  density.pheno %>%
  ggplot() +
  geom_point(aes(x = thrips, y = density),fill="grey",color="black",shape=21,size=2) +
  theme_bw() +
  geom_label_repel(aes(x = thrips, y=density,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  facet_wrap(~type, scale = "free", ncol = 1)+
  labs(x = "Tomato genotype rank for thrips survival (low to high survival)",y = "Trichome density (trichomes / mm2 leaf)") +
  scale_x_continuous(breaks=seq(0,19,1))+
  my.theme

################
# Save the plots
################
ggsave(filename = file.path("Figure3/whitefly_toxicity_vs_trichome_densities.pdf"),plot = p.whitelfies, width = 9,height = 21, units = "cm")
ggsave(filename = file.path("Figure3/thrips_toxicity_vs_trichome_densities.pdf"),plot = p.thrips, width = 9,height = 21, units = "cm")

#################
# session info
##################
writeLines(capture.output(sessionInfo()), "Figure3/sessionInfo.txt")
