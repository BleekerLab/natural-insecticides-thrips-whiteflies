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
density.pheno  = left_join(density.avg, pheno, by = "accession") %>% select(., -c("genotype", "N", "sd", "se", "ci"))
density.pheno$accession = as.factor(density.pheno$accession)


#Load acylsugars + make data tidy
acylsugars = read.delim("Figure6_heatmaps/acylsugars_normalised.tsv", header = T)
acylsugars = separate(acylsugars, col = accession, into = c("accession", "x"))  %>% select(., -c("x", "species")) #remove unwanted names

acylsugars = gather(acylsugars,
                    key = "metabolite",
                    value = "value",
                    -accession)

#fuse phenotype + acylsugar dataset
acylsugars.pheno = left_join(density.pheno, acylsugars, by = "accession")


## Acylsugars VS PHENOTYPES

p.acylsugars.wf =
acylsugars.pheno %>% filter(., type == "non.glandular" & metabolite %in% c("summed_glucose", "summed_sucrose", "summed_total")) %>%
  ggplot()+
  geom_point(aes(x = wf, y = value),fill="grey",color="black",shape=21,size=3) +
  theme_bw() +
  geom_label_repel(aes(x = wf, y= value,label=accession, fill=species),
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  facet_wrap(~metabolite, scale = "free", ncol = 1)+
  labs(x = "Tomato genotype rank for whitefly survival (low to high survival)",y = "Summed acylsugars (ion counts / mg leaf)") +
  scale_x_continuous(breaks=seq(0,19,1))

p.acylsugars.thrips =
  acylsugars.pheno %>% filter(., type == "non.glandular" & metabolite %in% c("summed_glucose", "summed_sucrose", "summed_total")) %>%
  ggplot()+
  geom_point(aes(x = thrips, y = value),fill="grey",color="black",shape=21,size=3) +
  theme_bw() +
  geom_label_repel(aes(x = thrips, y= value,label=accession, fill=species),
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  facet_wrap(~metabolite, scale = "free", ncol = 1)+
  labs(x = "Tomato genotype rank for thrips survival (low to high survival)",y = "Summed acylsugars (ion counts / mg leaf)") +
  scale_x_continuous(breaks=seq(0,19,1))

<<<<<<< HEAD

=======
<<<<<<< HEAD
>>>>>>> c05bfd541c8868926bc90d50d1ba180e1c9ed89d
## Acylsugars VS densitites
acylsugars.pheno %>% filter(metabolite == "summed_total") %>%
  ggplot()+
  geom_point(aes(x = density, y = value),fill="grey",color="black",shape=21,size=3) +
  theme_bw() +
  geom_label_repel(aes(x = density, y= value,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE,
                   alpha = 1,
                   segment.size = 0.2) +
  facet_wrap(~type, scale = "free", ncol = 1)+
  labs(x = "Trichome density (trichomes / mm2 leaf)",y = "Summed acylsugars (ion counts / mg leaf)")


<<<<<<< HEAD
=======

=======
>>>>>>> 7eb9db3633ac8c52b8c3c4f9ca919fa7af987233


>>>>>>> c05bfd541c8868926bc90d50d1ba180e1c9ed89d
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