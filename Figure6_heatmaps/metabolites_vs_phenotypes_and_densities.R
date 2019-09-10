library(RColorBrewer)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(Rmisc)

#############################
# Custom theme for plotting #
#############################

my.theme = theme(axis.text.x = element_text(color = "black", size = 6),
                 axis.text.y = element_text(color = "black", size = 6),
                 axis.title.x = element_text(color = "black", size = 8),
                 axis.title.y = element_text(color = "black", size = 8)
)
                 

#################################
# Load and shape phenotype data #
#################################

#Load trichome density data + Whitefly / Thrips ranked phenotypes
density = read.delim("Figure5_densities/trichome.counts.processed.tsv", header =T)
pheno = read.delim("Figure3_scatterplot/data4scatterplot.txt", header = T)
pheno = separate(pheno, col = sample, into = c("x", "y", "accession")) %>% #separate "sample" into accession names only so it can be used for fusing the datasets
  select(., c("accession", "wf", "thrips"))

#Calculate the average trichome densities per genotype (abaxial and adaxial are taken together)
density.avg = summarySE(density, measurevar = "density", groupvars = c("genotype", "accession", "species" ,"type", "color"))

#Fuse the two datasets to a new one
density.pheno  = left_join(density.avg, pheno, by = "accession") %>% select(., -c("genotype", "N", "sd", "se", "ci"))
density.pheno$accession = as.factor(density.pheno$accession)
#Label the trichome classes nicely
density.pheno$type = factor(density.pheno$type, levels = c("non.glandular", "typeIandIV", "typeVI"),
                            labels = c("Non glandular", "Type I/IV", "Type VI"))

##############
# Acylsugars #
##############

#Load acylsugars + make data tidy
acylsugars = read.delim("Figure6_heatmaps/acylsugars_normalised.tsv", header = T)
acylsugars = separate(acylsugars, col = accession, into = c("accession", "x"))  %>% select(., -c("x", "species")) #remove unwanted names

acylsugars = gather(acylsugars,
                    key = "metabolite",
                    value = "value",
                    -accession)

#fuse phenotype + acylsugar dataset
acylsugars.pheno = left_join(density.pheno, acylsugars, by = "accession")


## Acylsugars vs phenotypes ##

p.acylsugars.wf =
acylsugars.pheno %>% filter(., type == "Non glandular" & metabolite == "summed_total") %>% #also filter on trichome type because of data redundancy due to the tidy format
  ggplot()+
  geom_point(aes(x = wf, y = value),fill="grey",color="black",shape=21,size=2) +
  theme_bw() +
  geom_label_repel(aes(x = wf, y= value,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  labs(x = "Tomato genotype rank for whitefly survival (low to high survival)",y = "Summed acylsugars (ion counts / mg leaf)") +
  scale_x_continuous(breaks=seq(0,19,1))+
  my.theme

p.acylsugars.thrips =
  acylsugars.pheno %>% filter(., type == "Non glandular" & metabolite == "summed_total") %>% #also filter on trichome type because of data redundancy due to the tidy format
  ggplot()+
  geom_point(aes(x = thrips, y = value),fill="grey",color="black",shape=21,size=2) +
  theme_bw() +
  geom_label_repel(aes(x = thrips, y= value,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  labs(x = "Tomato genotype rank for thrips survival (low to high survival)",y = "Summed acylsugars (ion counts / mg leaf)") +
  scale_x_continuous(breaks=seq(0,19,1))+
  my.theme


## Acylsugars vs trichome type-I/IV densitites ##

p.acylsugars.densities =
acylsugars.pheno %>% filter(metabolite == "summed_total" & type == "Type I/IV") %>% 
  ggplot()+
  geom_point(aes(x = density, y = value),fill="grey",color="black",shape=21,size=2) +
  theme_bw() +
  geom_label_repel(aes(x = density, y= value,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE,
                   alpha = 1,
                   segment.size = 0.2) +
  labs(x = "Type-I/IV trichome density (trichomes / mm2 leaf)",y = "Summed acylsugars (ion counts / mg leaf)")+
  my.theme

#############
# Volatiles #
#############

#Load volatiles + make calculate summed volatiles
volatiles = read.delim("figure6_heatmaps/pheno_terpenoids.tsv", header = T)
volatiles = separate(volatiles, col = sample, into = c("x", "y", "accession")) #turn 'sample names' into accessions
volatiles = volatiles %>% select(., -c("x", "y", "wf", "thrips")) #remove unwanted columns
volatiles$summed_volatiles = rowSums(volatiles[,2:ncol(volatiles)])


#fuse phenotype + acylsugar dataset
volatiles.pheno = left_join(density.pheno, volatiles, by = "accession") %>% select(., "accession", "species", "type", "color", "density", "wf", "thrips", "summed_volatiles")
str(volatiles.pheno)

## Volatiles vs phenotypes
p.volatiles.wf =
  volatiles.pheno %>% filter(type == "Non glandular") %>% # filter on trichome type because of data redundancy due to the tidy format
  ggplot()+
  geom_point(aes(x = wf, y = summed_volatiles),fill="grey",color="black",shape=21,size=2) +
  theme_bw() +
  geom_label_repel(aes(x = wf, y= summed_volatiles,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  labs(x = "Tomato genotype rank for whitefly survival (low to high survival)",y = "Summed volatiles (ion counts / mg trichomes)") +
  scale_x_continuous(breaks=seq(0,19,1))+
  my.theme

p.volatiles.thrips =
  volatiles.pheno %>% filter(type == "Non glandular") %>% #filter on trichome type because of data redundancy due to the tidy format
  ggplot()+
  geom_point(aes(x = thrips, y = summed_volatiles),fill="grey",color="black",shape=21,size=2) +
  theme_bw() +
  geom_label_repel(aes(x = thrips, y= summed_volatiles,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  labs(x = "Tomato genotype rank for thrips survival (low to high survival)",y = "Summed volatiles (ion counts / mg trichomes)") +
  scale_x_continuous(breaks=seq(0,19,1))+
  my.theme

## Volatiles vs trichome type-VI densites

p.volatiles.densities =
  volatiles.pheno %>% filter(., type == "Type VI") %>%
  ggplot()+
  geom_point(aes(x = density, y = summed_volatiles),fill="grey",color="black",shape=21,size=2) +
  theme_bw() +
  geom_label_repel(aes(x = density, y= summed_volatiles,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE,
                   alpha = 1,
                   segment.size = 0.2) +
  labs(x = "Type-VI trichome density (trichomes / mm2 leaf)",y = "Summed volatiles (ion counts / mg trichomes)")+
  my.theme

##############
# SAVE PLOTS #
##############

ggsave(file = "Figure6_heatmaps/plots_metabolites_vs_phenotypes/acylsugars_vs_WF.pdf", plot = p.acylsugars.wf, width = 9, height = 9, units = "cm")
ggsave(file = "Figure6_heatmaps/plots_metabolites_vs_phenotypes/acylsugars_vs_thrips.pdf", plot = p.acylsugars.thrips, width = 9, height = 9, units = "cm")
ggsave(file = "Figure6_heatmaps/plots_metabolites_vs_phenotypes/acylsugars_vs_type_I_IV.pdf", plot = p.acylsugars.densities, width = 9, height = 9, units = "cm")
ggsave(file = "Figure6_heatmaps/plots_metabolites_vs_phenotypes/volatiles_vs_WF.pdf", plot = p.volatiles.wf, width = 9, height = 9, units = "cm")
ggsave(file = "Figure6_heatmaps/plots_metabolites_vs_phenotypes/volatiles_vs_thrips.pdf", plot = p.volatiles.thrips, width = 9, height = 9, units = "cm")
ggsave(file = "Figure6_heatmaps/plots_metabolites_vs_phenotypes/volatiles_vs_type_VI.pdf", plot = p.volatiles.densities, width = 9, height = 9, units = "cm")