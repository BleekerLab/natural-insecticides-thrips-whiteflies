library(RColorBrewer)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(Rmisc)

#Load trichome density data + Whitefly / Thrips ranked phenotypes
density = read.delim("Figure5_densities/trichome.counts.processed.tsv", header =T)
pheno = read.delim("Figure3_scatterplot/data4scatterplot.txt", header = T)
pheno = separate(pheno, col = sample, into = c("x", "y", "accession")) %>% #separate "sample" into accession names only so it can be used for fusing the datasets
  select(., c("accession", "wf", "thrips"))

#Calculate the average trichome densities per genotype (abaxial and adaxial are taken together)
density.avg = summarySE(density, measurevar = "density", groupvars = c("genotype", "accession", "species" ,"type", "color"))

#Fuse the two datasets to a new one
density.pheno  = left_join(density.avg, pheno, by = "accession")
density.pheno$accession = as.factor(density.pheno$accession)


#Load acylsugars + make data tidy
acylsugars = read.delim("Figure6_heatmaps/pheno_acylsugars.tsv", header = T) 
acylsugars = separate(acylsugars, col = sample, into = c("x", "y", "accession"))
acylsugars = acylsugars %>% select(., -c("x", "y", "wf", "thrips")) #remove unwanted columns
acylsugars = gather(acylsugars,
                    key = "acylsugar",
                    value = "acylsugar_value",
                    -accession)
acylsugars.pheno = left_join(acylsugars, density.pheno, by = "accession") %>% select(., c("accession", "acylsugar", "acylsugar_value", "species", "type", "density", "wf", "thrips"))

acylsugars.pheno %>%
  ggplot()+
  geom_point(aes(x = thrips, y = acylsugar_value),fill="grey",color="black",shape=21,size=3) +
  theme_bw() +
  geom_label_repel(aes(x = thrips, y= acylsugar_value,label=accession, fill=species),
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  geom_smooth(aes(x = thrips, y = acylsugar_value),
              method = lm,
              alpha = 0.25,
              linetype = "dotted",
              size = 1)+
  facet_wrap(~type, scale = "free", ncol = 1)+
  labs(x = "Tomato genotype rank for thrips survival (low to high survival)",y = "Trichome density (trichomes / mm2 leaf)") +
  scale_x_continuous(breaks=seq(0,19,1))


#Load volatiles + make it tidy
volatiles = read.delim("figure6_heatmaps/pheno_terpenoids.tsv", header = T)
volatiles = separate(volatiles, col = sample, into = c("x", "y", "accession"))
volatiles = volatiles %>% select(., -c("x", "y", "wf", "thrips")) #remove unwanted columns





#Calculate the average trichome densities per genotype (abaxial and adaxial are taken together)
density.avg = summarySE(density, measurevar = "density", groupvars = c("genotype", "accession", "species" ,"type", "color"))

#Fuse the two datasets to a new one
density.pheno  = left_join(density.avg, pheno, by = "accession")
density.pheno$accession = as.factor(density.pheno$accession)

#Label the trichome classes nicely
density.pheno$type = factor(density.pheno$type, levels = c("non.glandular", "typeIandIV", "typeVI"),
                            labels = c("Non glandular", "Type I/IV", "Type VI"))