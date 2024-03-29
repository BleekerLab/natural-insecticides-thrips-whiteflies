library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(ggplot2)

#############################
# Custom theme for plotting #
#############################

my.theme = theme(axis.text.x = element_text(color = "black", size = 6),
                 axis.text.y = element_text(color = "black", size = 6),
                 axis.title.x = element_text(color = "black", size = 8),
                 axis.title.y = element_text(color = "black", size = 8)
)

#Load Trichome density data + Whitefly / Thrips ranked phenotypes
density = read.delim("Figure_S3/trichome.counts.processed.tsv", header =T)
pheno = read.delim("Figure_1C/data4scatterplot.txt", header = T)
pheno = separate(pheno, col = sample, into = c("x", "y", "accession")) %>% #separate "sample" into accession names only so it can be used for fusing the datasets
  select(., c("accession", "wf", "thrips"))

#Calculate the average trichome densities per genotype (abaxial and adaxial are taken together)
density.avg = summarySE(density, measurevar = "density", groupvars = c("genotype", "accession", "species" ,"type", "color"))

#Fuse the two datasets to a new one
density.pheno  = left_join(density.avg, pheno, by = "accession")
density.pheno$accession = as.factor(density.pheno$accession)

#Label the trichome classes nicely
density.pheno$type = factor(density.pheno$type, levels = c("non.glandular", "typeIandIV", "typeVI"),
                            labels = c("Non glandular", "Type I/IV", "Type VI"))

# Plot Whiteflies vs. Trichome densities
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

# Plot Thrips vs. Trichome densities
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

##############
# SAVE PLOTS #
##############
ggsave(filename = file.path("Figure_S3/plots/WF_vs_densities.pdf"),plot = p.whitelfies, width = 9,height = 21, units = "cm")
ggsave(filename = file.path("Figure_S3/plots/Thrips_vs_densities.pdf"),plot = p.thrips, width = 9,height = 21, units = "cm")
